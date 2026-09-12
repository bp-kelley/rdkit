//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Correctness tests for the fragment conformer generator.
//  These are fast, deterministic and self-contained (no external data files or
//  reference tools).  The timing / measurement harnesses that sweep real
//  molecules -- including the crowded-junction strain sweep, which needs the
//  on-disk torsion library -- live in the standalone fragConfGenBench
//  benchmark, not here.
//
#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <Confgen/FragmentConfGen.h>
#include <Confgen/Search/InterFragScore.h>
#include <Confgen/Sampler/TorsionSampler.h>
#include <Confgen/Embedder/Fraglib.h>
#include <Confgen/Utils/ParamSentinels.h>
#include <Confgen/Utils/ParamsIO.h>
#include <Confgen/Search/RotorDriver.h>
#include <Confgen/Search/TreeSearch.h>
#include <Confgen/Joiner/FragmentJoiner.h>
#include <Confgen/Search/RigidRotorSearch.h>
#include <Confgen/Utils/ParamsIO.h>
#include <Confgen/Utils/TheobaldRmsd.h>
#include <Confgen/Utils/SymmetricRmsd.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <ForceField/ForceField.h>
#include <ForceField/MMFF/Params.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Geometry/Transform3D.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <set>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace RDKit;

TEST_CASE("SmirnoffTorsionSampler: parse offxml + energy-minima angles",
          "[smirnoff]") {
  // Synthetic fixture written by the test -- NOT the CC-BY OpenFF data.  One
  // 3-fold sp3 C-C proper torsion (E = k(1+cos(3*theta)) has minima at
  // 60/180/300 deg).
  const std::string offxml =
      "<SMIRNOFF version=\"0.3\">\n"
      " <ProperTorsions version=\"0.4\" "
      "potential=\"k*(1+cos(periodicity*theta-phase))\" default_idivf=\"auto\">\n"
      "  <Proper smirks=\"[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]\" periodicity1=\"3\" "
      "phase1=\"0.0 * degree\" k1=\"0.3 * mole**-1 * kilocalorie\" idivf1=\"1.0\" "
      "id=\"t1\"></Proper>\n"
      " </ProperTorsions>\n"
      "</SMIRNOFF>\n";
  const std::string path = std::string(std::tmpnam(nullptr)) + ".offxml";
  {
    std::ofstream out(path);
    out << offxml;
  }

  // n-butane, embedded so getAngles has a conformer + a real j-k frame.
  std::unique_ptr<RWMol> m(SmilesToMol("CCCC"));
  REQUIRE(m);
  MolOps::addHs(*m);
  DGeomHelpers::EmbedParameters ep = DGeomHelpers::ETKDGv3;
  ep.randomSeed = 42;
  REQUIRE(DGeomHelpers::EmbedMolecule(*m, ep) >= 0);

  // fallbackStepDeg=30 -> the unmatched fallback is a 12-point grid,
  // distinguishable from the 3 minima of the matched rule.
  SmirnoffTorsionSampler sampler(path, 30.0);
  // torsion over the four backbone carbons 0-1-2-3 (central bond 1-2)
  const auto angles = sampler.getAngles(*m, 0, 1, 2, 3);

  SECTION("periodicity-3 term yields three minima ~120 deg apart") {
    REQUIRE(angles.size() == 3);
    std::vector<double> a = angles;
    for (auto &x : a) x = std::fmod(x + 720.0, 360.0);
    std::sort(a.begin(), a.end());
    // consecutive gaps (circular) ~120 deg
    for (size_t idx = 0; idx < a.size(); ++idx) {
      double gap = a[(idx + 1) % a.size()] - a[idx];
      if (gap < 0) gap += 360.0;
      CHECK(std::abs(gap - 120.0) < 6.0);
    }
  }

  SECTION("no matching SMIRKS -> uniform fallback (not the 3-term minima)") {
    // a C-N bond has no C-C rule; ask for a bogus quartet with a non-C central
    // atom by reusing an aromatic molecule with no sp3 C-C-C-C torsion.
    std::unique_ptr<RWMol> benzene(SmilesToMol("c1ccccc1"));
    MolOps::addHs(*benzene);
    DGeomHelpers::EmbedMolecule(*benzene, ep);
    const auto fb = sampler.getAngles(*benzene, 0, 1, 2, 3);
    CHECK(fb.size() > 3);  // fell back to the coarse uniform grid
  }

  std::remove(path.c_str());
}

TEST_CASE("build round-trips the molecular graph", "[roundtrip]") {
  // The whole pipeline (cut -> embed -> molzip -> assemble) must reconstruct
  // the SAME molecule: every output conformer, reduced to its heavy-atom
  // canonical SMILES, must equal the input's largest-component heavy-atom
  // canonical SMILES (salts are dropped, Hs are cosmetic).  Guards against
  // lost/duplicated atoms, mis-zipped junctions and dropped stereochemistry.
  auto heavyLargest = [](const ROMol &mol) -> std::string {
    RWMol m(mol);
    MolOps::removeHs(m);
    auto comps = MolOps::getMolFrags(m, /*sanitizeFrags=*/true);
    const ROMol *best = &m;
    if (comps.size() > 1) {
      unsigned int bestH = 0;
      for (const auto &c : comps) {
        unsigned int h = 0;
        for (const auto a : c->atoms()) {
          if (a->getAtomicNum() > 1) {
            ++h;
          }
        }
        if (h > bestH) {
          bestH = h;
          best = c.get();
        }
      }
    }
    return MolToSmiles(*best);
  };

  FragmentConfGenParams p;
  p.embedding.numConfsPerFragment = 4;
  p.numOutputConfs = 6;
  p.randomSeed = 0xf00d;

  const std::vector<std::string> smis = {
      "c1ccccc1",                // rigid, one fragment
      "c1ccccc1-c1ccccc1",       // one rotor, two identical fragments
      "CCOC(=O)c1ccc(OC)cc1",    // a few rotors
      "CC(C)C[C@@H](N)C(=O)O",   // defined stereocentre (leucine)
      "CCOC(=O)c1ccc(OCCN)cc1",  // multi-rotor chain
      "Oc1ccc(C[C@@H](C)N)cc1",  // stereo + aromatic OH
      // chiral ROTOR: the stereocentre is defined by two different rotatable
      // exit vectors (benzyl vs pyridylmethyl) -- must survive cut/reassemble
      "OC(=O)[C@@H](Cc1ccccc1)Cc1ccncc1",
      "CCCCN(CCCC)CCCNC1=C2CCCCC2=NC3=C1C=CC=C3.OP(O)(O)=O",  // salt -> drop it
  };

  for (const auto &smi : smis) {
    std::unique_ptr<ROMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    const std::string want = heavyLargest(*mol);
    auto res = FragmentConfGen(p).build(*mol);
    INFO("input=" << smi << "  want=" << want);
    REQUIRE(!res.conformers.empty());
    for (const auto &c : res.conformers) {
      CHECK(heavyLargest(*c) == want);
    }
  }
}

TEST_CASE("FragmentConfGen basic fragmentation and generation",
          "[fragconfgen]") {
  FragmentConfGenParams params;
  params.embedding.numConfsPerFragment = 5;
  params.numOutputConfs = 10;
  params.randomSeed = 0xf00d;

  SECTION("rigid molecule -> one fragment") {
    std::unique_ptr<ROMol> mol(SmilesToMol("c1ccccc1"));
    REQUIRE(mol);
    FragmentConfGen gen(params);
    auto res = gen.build(*mol);
    CHECK(res.numLinkBonds == 0);
    CHECK(res.numFragments == 1);
    CHECK(!res.conformers.empty());
    CHECK(std::isfinite(res.bestEnergy));
  }

  SECTION("one rotatable bond -> two fragments") {
    // biphenyl: single rotatable bond between the two rings
    std::unique_ptr<ROMol> mol(SmilesToMol("c1ccccc1-c1ccccc1"));
    REQUIRE(mol);
    FragmentConfGen gen(params);
    auto res = gen.build(*mol);
    CHECK(res.numLinkBonds == 1);
    CHECK(res.numFragments == 2);
    CHECK(!res.conformers.empty());
    for (const auto &c : res.conformers) {
      CHECK(c->getNumConformers() == 1);
      double e = 0.0;
      CHECK(c->getPropIfPresent<double>(kEnergyProp, e));
      CHECK(std::isfinite(e));
    }
  }

  SECTION("salt / disconnected components are dropped") {
    // Main molecule + phosphoric-acid counterion; only the largest component
    // should be assembled, and it must still produce conformers.
    std::unique_ptr<ROMol> mol(
        SmilesToMol("CCCCN(CCCC)CCCNC1=C2CCCCC2=NC3=C1C=CC=C3.OP(O)(O)=O"));
    REQUIRE(mol);
    FragmentConfGen gen(params);
    auto res = gen.build(*mol);
    CHECK(!res.conformers.empty());
    // the phosphate (5 heavy atoms) must be gone from the output
    for (const auto &c : res.conformers) {
      for (const auto atom : c->atoms()) {
        CHECK(atom->getAtomicNum() != 15);  // no phosphorus
      }
    }
  }

  SECTION("multi-rotor chain -> several fragments") {
    std::unique_ptr<ROMol> mol(SmilesToMol("CCOC(=O)c1ccc(OCCN)cc1"));
    REQUIRE(mol);
    FragmentConfGen gen(params);
    std::vector<unsigned int> links;
    auto frags = gen.fragmentAndEmbed(*mol, &links);
    CHECK(links.size() >= 3);
    CHECK(frags.size() == links.size() + 1);
    // conformer count is governed per-fragment-class by
    // getDefaultFragmentParams() (numConfsPerFragment is no longer a sizing
    // knob); every fragment must embed.
    for (const auto &f : frags) {
      CHECK(f->getNumConformers() >= 1);
    }
  }
}

namespace {
// indices of the exit-vector dummy atoms (atomic number 0), in atom order
std::vector<unsigned int> exitIdxs(const ROMol &m) {
  std::vector<unsigned int> out;
  for (const auto a : m.atoms()) {
    if (a->getAtomicNum() == 0) {
      out.push_back(a->getIdx());
    }
  }
  return out;
}

// Cut `mol` into its rigid fragments exactly as FragmentConfGen does (largest
// component, link-bond cuts with isotope-labelled exit dummies 1000+i).
std::vector<ROMOL_SPTR> cutPieces(const ROMol &mol) {
  auto heavy = boost::make_shared<RWMol>(mol);
  MolOps::removeHs(*heavy);
  auto comps = MolOps::getMolFrags(*heavy, /*sanitizeFrags=*/true);
  if (comps.size() > 1) {
    ROMOL_SPTR largest;
    unsigned int best = 0;
    for (const auto &c : comps) {
      unsigned int h = 0;
      for (const auto a : c->atoms()) {
        if (a->getAtomicNum() > 1) {
          ++h;
        }
      }
      if (h > best) {
        best = h;
        largest = c;
      }
    }
    if (largest) {
      heavy = boost::make_shared<RWMol>(*largest);
    }
  }
  auto links = FragmentConfGen::findLinkBonds(*heavy);
  if (links.empty()) {
    return {heavy};
  }
  std::vector<std::pair<unsigned int, unsigned int>> labels;
  for (size_t i = 0; i < links.size(); ++i) {
    labels.emplace_back(1000u + i, 1000u + i);
  }
  std::unique_ptr<ROMol> fragged(
      MolFragmenter::fragmentOnBonds(*heavy, links, true, &labels));
  std::vector<ROMOL_SPTR> out;
  if (fragged) {
    for (auto &p : MolOps::getMolFrags(*fragged, /*sanitizeFrags=*/true)) {
      out.push_back(p);
    }
  }
  return out;
}

// A same-graph "query" variant of `piece`: relabel every exit vector to a new
// id and scramble the atom order, so it exercises both label- and order-
// independence of generateKey.
boost::shared_ptr<RWMol> relabelAndScramble(const ROMol &piece) {
  auto v = boost::make_shared<RWMol>(piece);
  for (auto a : v->atoms()) {
    if (a->getAtomicNum() == 0 && a->getIsotope() > 0) {
      a->setIsotope(a->getIsotope() + 1000);  // 1000+i -> 2000+i
    }
  }
  std::vector<unsigned int> rev(v->getNumAtoms());
  for (unsigned int i = 0; i < v->getNumAtoms(); ++i) {
    rev[i] = v->getNumAtoms() - 1 - i;  // reverse = a deterministic scramble
  }
  std::unique_ptr<ROMol> scr(MolOps::renumberAtoms(*v, rev));
  return boost::make_shared<RWMol>(*scr);
}
}  // namespace

TEST_CASE("Fraglib::generateKey is label-agnostic and canonicalizes order",
          "[generatekey]") {
  using RDKit::Fraglib;

  SECTION("exit-vector isotope labels do not change the key") {
    auto a = std::unique_ptr<RWMol>(new RWMol(*SmilesToMol("[1000*]c1ccccc1")));
    auto b = std::unique_ptr<RWMol>(new RWMol(*SmilesToMol("[1007*]c1ccccc1")));
    REQUIRE(a);
    REQUIRE(b);
    CHECK(Fraglib::generateKey(*a) == Fraglib::generateKey(*b));
    // ... but a genuinely different fragment differs
    auto c = std::unique_ptr<RWMol>(new RWMol(*SmilesToMol("[1000*]c1ccncc1")));
    REQUIRE(c);
    CHECK(Fraglib::generateKey(*a) != Fraglib::generateKey(*c));
  }

  SECTION("the fragment keeps its own isotopes after keying") {
    auto a =
        std::unique_ptr<RWMol>(new RWMol(*SmilesToMol("[1000*]CCO[1001*]")));
    REQUIRE(a);
    Fraglib::generateKey(*a);
    std::vector<unsigned int> isos;
    for (const auto atom : a->atoms()) {
      if (atom->getAtomicNum() == 0) {
        isos.push_back(atom->getIsotope());
      }
    }
    std::sort(isos.begin(), isos.end());
    CHECK(isos == std::vector<unsigned int>{1000, 1001});
  }

  SECTION("same graph -> identical canonical atom order (basis for remap)") {
    // Asymmetric linker so the two exit vectors are distinguishable, with
    // different labels AND different input atom orderings.
    auto a = std::unique_ptr<RWMol>(
        new RWMol(*SmilesToMol("[1000*]CCO[1001*]")));  // dummy-C-C-O-dummy
    auto b = std::unique_ptr<RWMol>(
        new RWMol(*SmilesToMol("[2001*]OCC[2000*]")));  // dummy-O-C-C-dummy
    REQUIRE(a);
    REQUIRE(b);
    const std::string ka = Fraglib::generateKey(*a);  // canonicalizes *a
    const std::string kb = Fraglib::generateKey(*b);  // canonicalizes *b
    CHECK(ka == kb);
    REQUIRE(a->getNumAtoms() == b->getNumAtoms());
    // atom-for-atom the two are now in the same canonical order
    for (unsigned int i = 0; i < a->getNumAtoms(); ++i) {
      CHECK(a->getAtomWithIdx(i)->getAtomicNum() ==
            b->getAtomWithIdx(i)->getAtomicNum());
    }
    CHECK(exitIdxs(*a) == exitIdxs(*b));

    // The remap: copy b's junction labels onto a's atoms by index; a must then
    // read out as the same labelled molecule as b (labels on the right atoms).
    for (unsigned int i = 0; i < a->getNumAtoms(); ++i) {
      if (a->getAtomWithIdx(i)->getAtomicNum() == 0) {
        a->getAtomWithIdx(i)->setIsotope(b->getAtomWithIdx(i)->getIsotope());
      }
    }
    CHECK(MolToSmiles(*a) == MolToSmiles(*b));
  }

  SECTION("remap=false does not renumber; key is the same either way") {
    auto a = std::unique_ptr<RWMol>(
        new RWMol(*SmilesToMol("[1001*]OCC[1000*]")));  // non-canonical order
    auto b = std::unique_ptr<RWMol>(new RWMol(*a));
    REQUIRE(a);
    const std::string kRemap = Fraglib::generateKey(*a, /*remap=*/true);
    const std::string kNo = Fraglib::generateKey(*b, /*remap=*/false);
    CHECK(kRemap == kNo);  // key is canonical regardless of reordering
  }

  SECTION("fragment with no exit vectors is handled") {
    auto a = std::unique_ptr<RWMol>(new RWMol(*SmilesToMol("c1ccccc1")));
    REQUIRE(a);
    const std::string k = Fraglib::generateKey(*a);
    CHECK(!k.empty());
    // idempotent: keying an already-canonical fragment gives the same key
    CHECK(Fraglib::generateKey(*a) == k);
  }

  SECTION(
      "chiral rotor (exit-defined stereocentre): generateKey now handles it") {
    // A stereocentre whose chirality is fixed by its two exit vectors.  The
    // chiral exits are marked [1*]/[2*] (canonically disambiguated), so the key
    // is now ORDER-STABLE (unlike the old collapse-to-[*] behaviour).
    auto a = std::unique_ptr<RWMol>(
        new RWMol(*SmilesToMol("[1000*][C@H](CC)[1001*]")));
    REQUIRE(a);
    std::vector<unsigned int> rev(a->getNumAtoms());
    for (unsigned int i = 0; i < a->getNumAtoms(); ++i) {
      rev[i] = a->getNumAtoms() - 1 - i;
    }
    std::unique_ptr<ROMol> aRev(MolOps::renumberAtoms(*a, rev));

    auto ag = std::unique_ptr<RWMol>(new RWMol(*a));
    auto agRev = std::unique_ptr<RWMol>(new RWMol(*aRev));
    const std::string k1 = Fraglib::generateKey(*ag);
    const std::string k2 = Fraglib::generateKey(*agRev);
    INFO("chiral-rotor key1=" << k1 << " key2=" << k2);
    CHECK(k1 == k2);                            // order-stable now
    CHECK(k1.find("1*") != std::string::npos);  // chiral exits got small marks
    CHECK(k1.find("2*") != std::string::npos);
  }

  SECTION("chiral rotor: relabel + scramble remaps correctly") {
    // The load-bearing check: a chiral rotor, relabelled AND atom-scrambled,
    // must key the same and remap onto the cached copy by plain index while
    // preserving the exit chirality.
    auto cached = std::unique_ptr<RWMol>(
        new RWMol(*SmilesToMol("[1000*][C@H](CC)[1001*]")));
    REQUIRE(cached);
    auto query = relabelAndScramble(*cached);
    const std::string kc = Fraglib::generateKey(*cached);
    const std::string kq = Fraglib::generateKey(*query);
    INFO("kc=" << kc << " kq=" << kq);
    CHECK(kc == kq);
    REQUIRE(cached->getNumAtoms() == query->getNumAtoms());
    const auto exits = exitIdxs(*cached);
    REQUIRE(exits == exitIdxs(*query));
    REQUIRE(exits.size() == 2);
    for (unsigned int idx : exits) {
      cached->getAtomWithIdx(idx)->setIsotope(
          query->getAtomWithIdx(idx)->getIsotope());
    }
    CHECK(MolToSmiles(*cached) == MolToSmiles(*query));
  }

  SECTION(
      "chiral rotor: one cached geometry serves BOTH enantiomers via remap") {
    // The two enantiomers of an exit-defined chiral rotor share a cache key
    // (one stored geometry).  That is correct only if remapping the query's
    // junction labels onto the cached copy reconstructs the RIGHT enantiomer --
    // because swapping the two exit labels on a stereocentre inverts it.
    // Verify: cache one enantiomer, ask for the other, and the relabelled
    // cached fragment must read out as the other enantiomer.
    auto cached = std::unique_ptr<RWMol>(
        new RWMol(*SmilesToMol("[1000*][C@H](CC)[1001*]")));  // stored
    auto query = std::unique_ptr<RWMol>(
        new RWMol(*SmilesToMol("[1000*][C@@H](CC)[1001*]")));  // requested
    REQUIRE(MolToSmiles(*cached) != MolToSmiles(*query));  // truly different
    const std::string kc = Fraglib::generateKey(*cached);  // canonicalizes
    const std::string kq = Fraglib::generateKey(*query);
    CHECK(kc == kq);  // same cache entry (dedup)
    REQUIRE(cached->getNumAtoms() == query->getNumAtoms());
    const auto exits = exitIdxs(*cached);
    REQUIRE(exits == exitIdxs(*query));
    // remap the query's labels onto the cached geometry, by index
    for (unsigned int idx : exits) {
      cached->getAtomWithIdx(idx)->setIsotope(
          query->getAtomWithIdx(idx)->getIsotope());
    }
    // the cached fragment now reads out as the QUERY enantiomer -- correct
    CHECK(MolToSmiles(*cached) == MolToSmiles(*query));
  }

  SECTION(
      "real big molecules: fragments remap correctly regardless of "
      "labels/order") {
    // Big, complex, stereo-rich molecules pulled from canonSmiles.long.smi.
    // For every rigid fragment they cut into, a relabelled + atom-scrambled
    // copy must (a) key the same and (b) after canonicalization line up
    // atom-for-atom so the query's junction labels transfer onto the cached
    // fragment by plain index -- reproducing exactly the labelled molecule.
    // This is the mapping the label-agnostic cache will rely on.
    const std::vector<std::string> bigs = {
        // steroid glycoside (many rings, many stereocentres, many rotors)
        "CO[CH]1C[CH](O[CH](C)[CH]1O[CH]2O[CH](CO)[CH](O)[CH](O)[CH]2O)O[CH]3CC"
        "[C]4(C=O)[CH]5CC[C]6(C)[CH](CC[C]6(O)[CH]5CC[C]4(O)C3)C7=CC(=O)OC7",
        // rifamycin-like macrocycle
        "CNC(=O)[CH](C)[CH]1N([CH]2CC[CH](O)[CH](C)O2)C(=O)C(=C(O)C=CC(C)=C[CH]"
        "(C)[CH]3O[C]4(C)O[CH](C=C[C]45CO5)[CH]3C)C1=O",
        // triarylmethane dye (asymmetric, sulfonate arms, charged)
        "CCN(CC1=CC(=CC=C1)S(O)(=O)=O)C2=CC=C(C=C2)C(=C3C=CC(C=C3)=[N+](CC)CC4="
        "CC(=CC=C4)S(O)(=O)=O)C5=C(C=CC=C5)S([O-])(=O)=O",
        // bis-quaternary diphenyl ester
        "CC[N+]1(CCCCC[N+]2(CC)CCC[CH](C2)OC(=O)C(C3=CC=CC=C3)C4=CC=CC=C4)CCC"
        "[CH](C1)OC(=O)C(C5=CC=CC=C5)C6=CC=CC=C6",
    };

    unsigned int fragsChecked = 0, multiExitChecked = 0;
    for (const auto &smi : bigs) {
      std::unique_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      auto pieces = cutPieces(*mol);
      REQUIRE(!pieces.empty());
      for (const auto &piece : pieces) {
        auto cached = boost::make_shared<RWMol>(*piece);  // the "stored" one
        auto query = relabelAndScramble(*piece);          // the "lookup" one

        const std::string kc = Fraglib::generateKey(*cached);  // canonicalizes
        const std::string kq = Fraglib::generateKey(*query);   // canonicalizes
        INFO("fragment key(cached)=" << kc << " key(query)=" << kq);
        CHECK(kc == kq);  // label- and order-agnostic

        REQUIRE(cached->getNumAtoms() == query->getNumAtoms());
        // canonicalization put them in the same order: atom-for-atom match
        for (unsigned int i = 0; i < cached->getNumAtoms(); ++i) {
          REQUIRE(cached->getAtomWithIdx(i)->getAtomicNum() ==
                  query->getAtomWithIdx(i)->getAtomicNum());
        }
        REQUIRE(exitIdxs(*cached) == exitIdxs(*query));

        // The remap: transfer the query's junction labels onto the cached
        // fragment by index; it must now read out as the query molecule.
        const auto exits = exitIdxs(*cached);
        for (unsigned int idx : exits) {
          cached->getAtomWithIdx(idx)->setIsotope(
              query->getAtomWithIdx(idx)->getIsotope());
        }
        CHECK(MolToSmiles(*cached) == MolToSmiles(*query));

        ++fragsChecked;
        if (exits.size() > 1) {
          ++multiExitChecked;
        }
      }
    }
    INFO("fragments checked = " << fragsChecked << " (multi-exit = "
                                << multiExitChecked << ")");
    CHECK(fragsChecked > 10);  // the sweep actually exercised real fragments
    CHECK(multiExitChecked >= 1);  // ...including linkers with >1 exit vector
  }
}

TEST_CASE("FragmentConfGenResult::isValid flags MMFF-untypeable molecules",
          "[boron]") {
  // MMFF cannot type boron, so FragmentConfGen cannot score/assemble these
  // benzoxaborole drugs; isValid() must report that up front (build() returns
  // no conformers for them).  Typeable drug-like molecules must report valid.
  FragmentConfGen gen;  // default params

  SECTION("boron molecules are not valid") {
    // tavaborole (rigid) and crisaborole (has a rotatable bond)
    for (const auto *smi :
         {"B1(OCc2c1ccc(c2)F)O", "N#Cc1ccc(Oc2ccc3c(c2)COB3O)cc1"}) {
      std::unique_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      INFO("smiles = " << smi);
      auto res = gen.build(*mol);
      CHECK(!res.isValid());
      // and build() indeed yields nothing for an invalid molecule
      CHECK(res.conformers.empty());
    }
  }

  SECTION("typeable molecules are valid") {
    for (const auto *smi :
         {"c1ccccc1", "CCOC(=O)c1ccc(OC)cc1", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"}) {
      std::unique_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      INFO("smiles = " << smi);
      CHECK(gen.build(*mol).resultType() == FragConfGenResultType::OK);
    }
  }

  SECTION("validity follows the largest component (salt of a boron drug)") {
    // a typeable counterion does not make a boron drug typeable...
    std::unique_ptr<ROMol> saltOfBoron(
        SmilesToMol("B1(OCc2c1ccc(c2)F)O.[Na+].[Cl-]"));
    REQUIRE(saltOfBoron);
    CHECK(gen.build(*saltOfBoron).resultType() ==
          FragConfGenResultType::FF_FAIL);
    // ...and a small untypeable counterion does not spoil a typeable main
    // molecule
    std::unique_ptr<ROMol> saltOfDrug(
        SmilesToMol("CC(C)Cc1ccc(cc1)C(C)C(=O)O.OB(O)O"));
    REQUIRE(saltOfDrug);
    CHECK(gen.build(*saltOfDrug).resultType() == FragConfGenResultType::OK);
  }

  SECTION("Fraglib still embeds boron fragments (DG/ETKDG need no MMFF)") {
    std::unique_ptr<ROMol> mol(SmilesToMol("B1(OCc2c1ccc(c2)F)O"));
    REQUIRE(mol);
    FraglibParams fp;
    fp.numConfsPerFragment = 4;
    fp.randomSeed = 0xf00d;
    fp.minimizeMode =
        FragmentMinimize::Full;  // must be skipped gracefully, not crash
    Fraglib lib(fp);
    ROMOL_SPTR e = lib.get(*mol);
    REQUIRE(e);
    CHECK(e->getNumConformers() >=
          1);  // per-class count; boron embeds via DG/ETKDG
  }
}

TEST_CASE("Fraglib embed cache", "[fraglib]") {
  FraglibParams params;
  params.numConfsPerFragment = 4;
  params.randomSeed = 0xf00d;

  auto mol = [](const std::string &smi) {
    return std::unique_ptr<ROMol>(SmilesToMol(smi));
  };

  auto hasIso = [](const ROMol &m, unsigned int iso) {
    for (const auto a : m.atoms()) {
      if (a->getAtomicNum() == 0 && a->getIsotope() == iso) {
        return true;
      }
    }
    return false;
  };

  SECTION(
      "embeds on first get, caches on second (returns a fresh owned copy)") {
    Fraglib lib(params);
    CHECK(lib.size() == 0);

    auto benzene = mol("c1ccccc1");
    REQUIRE(benzene);
    ROMOL_SPTR m = lib.get(*benzene);
    REQUIRE(m);
    CHECK(m->getNumConformers() >=
          1);  // per-class count from getDefaultFragmentParams
    CHECK(lib.size() == 1);
    CHECK(lib.numFragmentConfs(*benzene) == m->getNumConformers());

    // The audit path must be a true lookup: asking about a miss cannot embed
    // it or change the cache size (benchmark timing relies on this property).
    auto propane = mol("CCC");
    REQUIRE(propane);
    CHECK_FALSE(lib.numFragmentConfs(*propane));
    CHECK(lib.size() == 1);

    // second request is a cache HIT (no new entry) but a distinct owned copy
    ROMOL_SPTR again = lib.get(*benzene);
    REQUIRE(again);
    CHECK(again.get() != m.get());  // fresh copy, not the cached master
    CHECK(again->getNumConformers() == m->getNumConformers());
    CHECK(lib.size() == 1);

    // a non-canonical spelling of the SAME molecule hits the same entry
    auto benzene2 = mol("C1=CC=CC=C1");
    REQUIRE(benzene2);
    REQUIRE(lib.get(*benzene2));
    CHECK(lib.size() == 1);

    // a different fragment adds an entry
    auto etoh = mol("CCO");
    REQUIRE(etoh);
    REQUIRE(lib.get(*etoh));
    CHECK(lib.size() == 2);
  }

  SECTION("classParams injects a per-class recipe override") {
    // cyclohexane classifies as SmallRing.  With the default recipe it keeps a
    // pucker set; an injected override that caps SmallRing to a single kept
    // conformer must win.
    auto ring = mol("C1CCCCC1");
    REQUIRE(ring);

    Fraglib def(params);  // empty classParams -> built-in recipe
    ROMOL_SPTR d = def.get(*ring);
    REQUIRE(d);
    CHECK(d->getNumConformers() >= 1);

    FraglibParams ov = params;
    FragmentParams sr = getDefaultFragmentParams(FragmentClass::SmallRing);
    sr.maxConfs = 1;  // keep exactly one
    ov.setClassParam(FragmentClass::SmallRing,
                     sr);  // helper == classParams[cls] = sr
    Fraglib lib(ov);
    ROMOL_SPTR m = lib.get(*ring);
    REQUIRE(m);
    CHECK(m->getNumConformers() == 1);  // the override took effect
    CHECK(m->getNumConformers() <= d->getNumConformers());

    // the override is keyed by CLASS: an Acyclic override leaves the SmallRing
    // fragment alone.
    FraglibParams other = params;
    FragmentParams ac = getDefaultFragmentParams(FragmentClass::Acyclic);
    ac.maxConfs = 1;
    other.classParams[FragmentClass::Acyclic] = ac;
    Fraglib lib2(other);
    ROMOL_SPTR m2 = lib2.get(*ring);
    REQUIRE(m2);
    CHECK(m2->getNumConformers() == d->getNumConformers());  // unaffected
  }

  SECTION(
      "exit-vector fragments are label-agnostic; copy keeps the caller's id") {
    // isotope-labelled exit vectors, as produced by the joiner's cuts
    Fraglib lib(params);
    auto frag = mol("[1000*]c1ccccc1");
    REQUIRE(frag);
    ROMOL_SPTR m = lib.get(*frag);
    REQUIRE(m);
    CHECK(m->getNumConformers() >=
          1);                 // per-class count from getDefaultFragmentParams
    CHECK(hasIso(*m, 1000));  // returned copy carries the caller's junction id
    CHECK(lib.size() == 1);

    // a DIFFERENT exit-vector isotope is now the SAME fragment (label-agnostic
    // key): no new entry, and the returned copy is relabelled to the new id.
    auto frag2 = mol("[1001*]c1ccccc1");
    REQUIRE(frag2);
    ROMOL_SPTR m2 = lib.get(*frag2);
    REQUIRE(m2);
    CHECK(lib.size() == 1);    // shared entry (was 2 before generateKey)
    CHECK(hasIso(*m2, 1001));  // relabelled to the caller's id
    CHECK(!hasIso(*m2, 1000));
  }

  SECTION("unembeddable fragment returns null") {
    Fraglib lib(params);
    // a lone dummy has nothing to embed once promoted to H
    auto empty = mol("[1000*]");
    REQUIRE(empty);
    ROMOL_SPTR r = lib.get(*empty);
    if (!r) {
      CHECK(lib.size() == 0);
    }
  }

  SECTION("serialize / initFromStream round-trips the cache with conformers") {
    Fraglib lib(params);
    auto benzene = mol("c1ccccc1");
    auto etohFrag = mol("[1000*]CCO");
    REQUIRE(benzene);
    REQUIRE(etohFrag);
    ROMOL_SPTR a = lib.get(*benzene);
    ROMOL_SPTR b = lib.get(*etohFrag);
    REQUIRE(a);
    REQUIRE(b);
    const unsigned int na = a->getNumAtoms(), ca = a->getNumConformers();
    const unsigned int nb = b->getNumAtoms(), cb = b->getNumConformers();

    std::stringstream ss;
    lib.serialize(ss);

    // Construct the target with DIFFERENT params to prove initFromStream
    // restores the serialized embedding params, not the constructor's.
    FraglibParams other;
    other.numConfsPerFragment = 99;
    other.fragmentEmbedMode = FragmentEmbedMode::DG;
    other.randomSeed = 12345;
    Fraglib loaded(other);
    loaded.initFromStream(ss);
    CHECK(loaded.size() == 2);
    CHECK(loaded.params().numConfsPerFragment == params.numConfsPerFragment);
    CHECK(loaded.params().fragmentEmbedMode == params.fragmentEmbedMode);
    CHECK(loaded.params().randomSeed == params.randomSeed);
    CHECK(loaded.numFragmentConfs(*benzene) == ca);
    CHECK(loaded.numFragmentConfs(*etohFrag) == cb);

    // reloaded fragments carry the same atoms and conformers, and re-fetching
    // by the same molecule is a hit (no re-embed), so size is unchanged
    ROMOL_SPTR a2 = loaded.get(*benzene);
    ROMOL_SPTR b2 = loaded.get(*etohFrag);
    REQUIRE(a2);
    REQUIRE(b2);
    CHECK(a2->getNumAtoms() == na);
    CHECK(a2->getNumConformers() == ca);
    CHECK(b2->getNumAtoms() == nb);
    CHECK(b2->getNumConformers() == cb);
    CHECK(loaded.size() == 2);

    // coordinates survived
    const RDGeom::Point3D p = a->getConformer(0).getAtomPos(0);
    const RDGeom::Point3D p2 = a2->getConformer(0).getAtomPos(0);
    CHECK((p - p2).length() < 1e-6);
  }

  SECTION("fragmentAndEmbed uses an injected shared cache") {
    // A shared cache passed through params should be populated by the fragment
    // embedding step and reused; the cut fragments must land in it.
    auto sharedLib = std::make_shared<Fraglib>(params);
    FragmentConfGenParams p;
    p.embedding.numConfsPerFragment = params.numConfsPerFragment;
    p.randomSeed = params.randomSeed;
    p.fraglib = sharedLib;

    auto biphenyl = mol("c1ccccc1-c1ccccc1");
    REQUIRE(biphenyl);
    std::vector<unsigned int> links;
    auto frags = FragmentConfGen(p).fragmentAndEmbed(*biphenyl, &links);
    CHECK(!frags.empty());
    // two phenyl fragments were CUT (chunks) -> they go to the shared cache so
    // a future conf gen reuses them
    CHECK(sharedLib->size() >= 1);
  }

  SECTION("a whole, uncut molecule ALSO goes through the one library") {
    // There is exactly one fragment cache per generator.  A rigid molecule (no rotatable bonds)
    // is embedded as a single unit, and that unit is cached like any other: rigid fragments
    // recur across inputs, so caching them is the point.  This previously used a separate
    // private cache, which meant the same molecule was re-embedded on every call and left one
    // code path silently using different embedding parameters from the rest of the pipeline.
    auto sharedLib = std::make_shared<Fraglib>(params);
    FragmentConfGenParams p;
    p.embedding.numConfsPerFragment = params.numConfsPerFragment;
    p.randomSeed = params.randomSeed;
    p.fraglib = sharedLib;

    auto benzene = mol("c1ccccc1");  // no rotatable bonds -> not fragmented
    REQUIRE(benzene);
    std::vector<unsigned int> links;
    auto frags = FragmentConfGen(p).fragmentAndEmbed(*benzene, &links);
    CHECK(!frags.empty());
    CHECK(links.empty());              // it really is the uncut path
    CHECK(sharedLib->size() == 1);     // and it landed in the shared cache

    // a second pass must HIT the cache rather than grow it
    auto again = FragmentConfGen(p).fragmentAndEmbed(*benzene, &links);
    CHECK(!again.empty());
    CHECK(sharedLib->size() == 1);
  }

  SECTION("a generator with no injected library still gets one") {
    // The constructor guarantees d_params.fraglib is never null, so every path can use it
    // unconditionally -- that invariant is what removed the private-cache branch.
    FragmentConfGenParams p;
    p.embedding.numConfsPerFragment = params.numConfsPerFragment;
    p.randomSeed = params.randomSeed;
    const FragmentConfGen gen(p);          // no p.fraglib set
    auto biphenyl = mol("c1ccccc1-c1ccccc1");
    REQUIRE(biphenyl);
    std::vector<unsigned int> links;
    CHECK(!gen.fragmentAndEmbed(*biphenyl, &links).empty());
  }
}

namespace {
// Embed one deterministic 3D conformer with explicit Hs; returns the H-added
// mol.
std::unique_ptr<ROMol> embed3D(const std::string &smi,
                               unsigned int seed = 0xC0FFEE) {
  std::unique_ptr<RWMol> m(SmilesToMol(smi));
  if (!m) {
    return nullptr;
  }
  MolOps::addHs(*m);
  DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
  ps.randomSeed = static_cast<int>(seed);
  if (DGeomHelpers::EmbedMolecule(*m, ps) < 0) {
    return nullptr;
  }
  return std::unique_ptr<ROMol>(m.release());
}
}  // namespace

TEST_CASE("RotorDriver drives dihedrals on coordinates only", "[rotordriver]") {
  SECTION("drive sets the reported dihedral to match RDKit's own measurement") {
    // phenethylamine-ish: an aromatic-CH2 junction plus two acyclic rotors, so
    // we exercise several rotors with differently-sized moving sides.
    auto m = embed3D("c1ccccc1CCN");
    REQUIRE(m);
    RotorDriver rd(*m);
    CHECK(rd.numRotors() >= 2);

    // A scratch conformer we fill from the driver's buffer to measure with
    // RDKit.
    Conformer probe(rd.numAtoms());

    for (unsigned int r = 0; r < rd.numRotors(); ++r) {
      auto t = rd.torsion(r);
      // j and k are the axis; they must never be in the moving set.
      const auto &moving = rd.movingAtoms(r);
      for (unsigned int a : moving) {
        CHECK(a != t[1]);
        CHECK(a != t[2]);
      }
      // l rotates with the k-side; i (fixed side) does not.
      CHECK(std::find(moving.begin(), moving.end(), t[3]) != moving.end());
      CHECK(std::find(moving.begin(), moving.end(), t[0]) == moving.end());

      for (double target : {-120.0, -30.0, 45.0, 150.0}) {
        rd.drive(r, target);
        // driver's own read-back
        CHECK(rd.dihedralDeg(r) == Catch::Approx(target).margin(1e-4));
        // independent read-back through RDKit on the same coordinates and
        // quartet
        const auto &pos = rd.positions();
        for (unsigned int a = 0; a < rd.numAtoms(); ++a) {
          probe.setAtomPos(
              a, RDGeom::Point3D(pos[3 * a], pos[3 * a + 1], pos[3 * a + 2]));
        }
        double viaRDKit =
            MolTransforms::getDihedralDeg(probe, t[0], t[1], t[2], t[3]);
        CHECK(viaRDKit == Catch::Approx(target).margin(1e-3));
      }
    }
  }

  SECTION("rotation is rigid: bonds within the moving set are preserved") {
    auto m = embed3D("CCCCCO");
    REQUIRE(m);
    RotorDriver rd(*m);
    REQUIRE(rd.numRotors() >= 1);

    // record a bond length between two moving atoms (if the moving set has one)
    for (unsigned int r = 0; r < rd.numRotors(); ++r) {
      const auto &moving = rd.movingAtoms(r);
      if (moving.size() < 2) {
        continue;
      }
      auto distIn = [&](unsigned int a, unsigned int b) {
        const auto &p = rd.positions();
        double dx = p[3 * a] - p[3 * b];
        double dy = p[3 * a + 1] - p[3 * b + 1];
        double dz = p[3 * a + 2] - p[3 * b + 2];
        return std::sqrt(dx * dx + dy * dy + dz * dz);
      };
      unsigned int a = moving[0], b = moving[1];
      double before = distIn(a, b);
      rd.drive(r, 99.0);
      CHECK(distIn(a, b) == Catch::Approx(before).margin(1e-6));
    }
  }

  SECTION(
      "rotor containment nests: parent moving sets superset their children") {
    // A benzene core (largest fragment -> root) with a long acyclic tail; every
    // rotor down the tail is contained in the ones nearer the ring.
    auto m = embed3D("c1ccccc1CCCCO");
    REQUIRE(m);
    RotorDriver rd(*m);
    REQUIRE(rd.numRotors() >= 3);

    auto asSet = [](const std::vector<unsigned int> &v) {
      return std::set<unsigned int>(v.begin(), v.end());
    };

    // Ground truth: s is carried by r iff every atom s moves is also moved by
    // r.
    for (unsigned int r = 0; r < rd.numRotors(); ++r) {
      auto movR = asSet(rd.movingAtoms(r));
      std::set<unsigned int> expected;
      for (unsigned int s = 0; s < rd.numRotors(); ++s) {
        if (s == r) {
          continue;
        }
        auto movS = asSet(rd.movingAtoms(s));
        bool subset =
            std::includes(movR.begin(), movR.end(), movS.begin(), movS.end());
        if (subset && movS.size() < movR.size()) {
          expected.insert(s);
        }
      }
      CHECK(asSet(rd.rotorSubsets(r)) == expected);
      // reported dependents are ordered largest-movement-first
      const auto &dep = rd.rotorSubsets(r);
      for (size_t x = 1; x < dep.size(); ++x) {
        CHECK(rd.movingAtoms(dep[x - 1]).size() >=
              rd.movingAtoms(dep[x]).size());
      }
    }

    // The widest-movement rotor is a ring->tail bond and carries every other
    // rotor.
    const auto &order = rd.numRotorAtoms();
    REQUIRE(order.size() == rd.numRotors());
    for (size_t x = 1; x < order.size(); ++x) {
      CHECK(rd.movingAtoms(order[x - 1]).size() >=
            rd.movingAtoms(order[x]).size());
    }
    CHECK(rd.rotorSubsets(order.front()).size() == rd.numRotors() - 1);
    // root is a benzene carbon (the largest fragment), never in any moving set
    for (unsigned int r = 0; r < rd.numRotors(); ++r) {
      const auto &mv = rd.movingAtoms(r);
      CHECK(std::find(mv.begin(), mv.end(), rd.rootAtom()) == mv.end());
    }
  }

  SECTION("sibling branches do not contain each other") {
    // Two distinct substituents off a benzene core: their rotors are siblings,
    // so neither is contained in the other.
    auto m = embed3D("OCc1ccc(CN)cc1");
    REQUIRE(m);
    RotorDriver rd(*m);
    REQUIRE(rd.numRotors() >= 2);
    // find the two CH2-heteroatom rotors (moving set of size 1 heavy + Hs);
    // they are on opposite ends of the ring and must be mutually
    // non-containing.
    for (unsigned int r = 0; r < rd.numRotors(); ++r) {
      for (unsigned int s : rd.rotorSubsets(r)) {
        // containment must be strict and asymmetric
        const auto &sc = rd.rotorSubsets(s);
        CHECK(std::find(sc.begin(), sc.end(), r) == sc.end());
      }
    }
  }

  SECTION("driving a parent carries a contained child rigidly") {
    auto m = embed3D("c1ccccc1CCCCO");
    REQUIRE(m);
    RotorDriver rd(*m);
    unsigned int parent = rd.numRotorAtoms().front();
    REQUIRE(!rd.rotorSubsets(parent).empty());
    unsigned int child = rd.rotorSubsets(parent).front();

    double childBefore = rd.dihedralDeg(child);
    rd.drive(parent, 73.0);  // big swing on the parent
    // the child's own dihedral is preserved because it moved as a rigid body
    CHECK(rd.dihedralDeg(child) == Catch::Approx(childBefore).margin(1e-6));
  }

  SECTION("scoring delegates to the supplied callables") {
    auto m = embed3D("c1ccccc1CCN");
    REQUIRE(m);
    // fullFF = sum of x; fragmentFF = sum of y.  Distinct so we can tell which
    // ran.
    RotorDriver::ScoreFn sumX = [](const double *p, unsigned int n) {
      double s = 0.0;
      for (unsigned int a = 0; a < n; ++a) {
        s += p[3 * a];
      }
      return s;
    };
    RotorDriver::ScoreFn sumY = [](const double *p, unsigned int n) {
      double s = 0.0;
      for (unsigned int a = 0; a < n; ++a) {
        s += p[3 * a + 1];
      }
      return s;
    };

    // no FF -> NaN
    RotorDriver bare(*m);
    CHECK(std::isnan(bare.score()));

    // fragmentFF present -> score() prefers it over fullFF
    RotorDriver rd(*m, -1, sumX, sumY);
    double drivenScore = rd.drive(0, 60.0);
    CHECK(drivenScore == Catch::Approx(rd.scoreFragment()));
    CHECK(rd.scoreFragment() != Catch::Approx(rd.scoreFull()));
    CHECK(rd.score() == Catch::Approx(rd.scoreFragment()));

    // fullFF only -> score() uses it
    RotorDriver rf(*m, -1, sumX, {});
    rf.drive(0, 60.0);
    CHECK(rf.score() == Catch::Approx(rf.scoreFull()));
  }
}

TEST_CASE("RotorTree searches rotor angles coordinate-only", "[rotortree]") {
  // MMFF energy of a flat coordinate buffer, matching mol's atom order.
  auto mmffE = [](const ROMol &m, const std::vector<double> &pos) {
    auto fn = makeFullFFScoreFn(m, /*electrostatics=*/false);
    REQUIRE(fn);
    return fn(pos.data(), m.getNumAtoms());
  };

  SECTION("finds a lower-energy ensemble than the seed, no molecule rebuild") {
    // butane-diol-ish flexible chain: a few real rotors to search.
    auto m = embed3D("OCCCCO");
    REQUIRE(m);
    auto score = makeFullFFScoreFn(*m, false);
    REQUIRE(score);
    RotorDriver drv(*m, -1, score);
    REQUIRE(drv.numRotors() >= 3);

    double seedE = drv.score();
    unsigned int nAtomsBefore = drv.numAtoms();

    RigidRotorSearchParams p;
    p.tree.beamWidth = 40;
    RotorTree tree(drv, p);
    auto res = tree.search();

    REQUIRE(!res.empty());
    // topology is untouched: still the same atom count / same buffer size
    CHECK(drv.numAtoms() == nAtomsBefore);
    CHECK(res.front().coords.size() == 3u * nAtomsBefore);
    // results are sorted best-first and the best beats the seed geometry
    for (size_t i = 1; i < res.size(); ++i) {
      CHECK(res[i - 1].score <= res[i].score + 1e-6);
    }
    CHECK(res.front().score <= seedE + 1e-6);
    // the reported score is the true MMFF energy of the reported coordinates
    CHECK(mmffE(*m, res.front().coords) ==
          Catch::Approx(res.front().score).margin(1e-4));
  }

  SECTION("produces several distinct rotamers within the energy window") {
    auto m = embed3D("OCCCCO");
    REQUIRE(m);
    RotorDriver drv(*m, -1, makeFullFFScoreFn(*m, false));
    RigidRotorSearchParams p;
    p.tree.beamWidth = 30;
    p.energyWindow = 12.0;
    RotorTree tree(drv, p);
    auto res = tree.search();
    REQUIRE(res.size() >= 2);
    // best-first, all within the window of the best
    for (const auto &r : res) {
      CHECK(r.score <= res.front().score + p.energyWindow + 1e-6);
    }
    // the top two are geometrically different (angle dedup did its job)
    double dmax = 0.0;
    for (unsigned int a = 0; a < drv.numAtoms(); ++a) {
      double dx = res[0].coords[3 * a] - res[1].coords[3 * a];
      double dy = res[0].coords[3 * a + 1] - res[1].coords[3 * a + 1];
      double dz = res[0].coords[3 * a + 2] - res[1].coords[3 * a + 2];
      dmax = std::max(dmax, std::sqrt(dx * dx + dy * dy + dz * dz));
    }
    CHECK(dmax > 0.2);
  }

  SECTION("beam search beats a single-rotor-at-a-time greedy pass") {
    // greedy (beamWidth 1) can get trapped; a wider beam should reach <= its
    // energy.
    auto m = embed3D("OCCCCCCO");
    REQUIRE(m);
    RotorDriver drv(*m, -1, makeFullFFScoreFn(*m, false));

    RigidRotorSearchParams greedyP;
    greedyP.tree.beamWidth = 1;
    double greedy = RotorTree(drv, greedyP).search().front().score;

    // reseed the driver's coordinates for a fair comparison
    RotorDriver drv2(*m, -1, makeFullFFScoreFn(*m, false));
    RigidRotorSearchParams wideP;
    wideP.tree.beamWidth = 60;
    double wide = RotorTree(drv2, wideP).search().front().score;

    CHECK(wide <= greedy + 1e-6);
  }
}

TEST_CASE(
    "inter-fragment energy decomposition is exact under junction rotation",
    "[rotortree][decomp]") {
  // Cut a molecule at ONE rotatable bond -> two fragments.  Driving that
  // junction rigidly rotates one fragment, so the full MMFF energy must change
  // by EXACTLY the change in the inter-fragment field (intra-fragment energy is
  // invariant).  This is the guarantee behind "use the stored fragment energy,
  // only rescore the inter-fragment rotors".
  auto check = [](const std::string &smi) {
    INFO("smiles=" << smi);
    auto m = embed3D(smi);
    REQUIRE(m);
    // the junction = the driver's first rotor; partition atoms by its two sides
    RotorDriver drv(*m);
    REQUIRE(drv.numRotors() >= 1);
    auto t = drv.torsion(0);
    unsigned int j = t[1], k = t[2];
    // fragOf: side reachable from k without crossing j = fragment 1, else 0
    std::vector<int> fragOf(m->getNumAtoms(), 0);
    {
      std::vector<char> seen(m->getNumAtoms(), 0);
      seen[j] = 1;
      std::vector<unsigned int> stk{k};
      seen[k] = 1;
      while (!stk.empty()) {
        unsigned int a = stk.back();
        stk.pop_back();
        fragOf[a] = 1;
        for (const auto nbr : m->atomNeighbors(m->getAtomWithIdx(a))) {
          if (!seen[nbr->getIdx()]) {
            seen[nbr->getIdx()] = 1;
            stk.push_back(nbr->getIdx());
          }
        }
      }
    }
    std::vector<std::pair<unsigned int, unsigned int>> junctions = {{j, k}};

    auto full = makeFullFFScoreFn(*m, false);
    auto inter = makeInterFragmentScoreFn(*m, fragOf, junctions, false);
    REQUIRE(full);
    REQUIRE(inter);

    RotorDriver d(*m, -1, full);
    d.drive(0, -60.0);
    double full0 = full(d.positions().data(), d.numAtoms());
    double int0 = inter(d.positions().data(), d.numAtoms());
    for (double ang : {0.0, 60.0, 120.0, 175.0}) {
      d.drive(0, ang);
      double fullA = full(d.positions().data(), d.numAtoms());
      double intA = inter(d.positions().data(), d.numAtoms());
      // change in full energy == change in inter-fragment energy
      CHECK((fullA - full0) == Catch::Approx(intA - int0).margin(1e-4));
    }
  };
  check("c1ccccc1-c1ccccc1");  // biphenyl: torsion + ortho-H clash across
                               // junction
  check("c1ccccc1CC");         // ethylbenzene: aryl-CH2 junction
  check("CC(=O)OCC");          // ester chain junction
}

TEST_CASE("coordinate-only fragment placement reassembles a cut molecule",
          "[joiner][placement]") {
  // Cut a molecule at one rotatable bond, rigidly DISPLACE the child fragment
  // (as if it were embedded in its own frame), then place it back
  // coordinate-only and drive the junction torsion to the original dihedral.
  // The result must reproduce the original geometry -- proving place + drive
  // reassembles without any molecule.
  auto reassemble = [](const std::string &smi) {
    INFO("smiles=" << smi);
    auto m = embed3D(smi);
    REQUIRE(m);
    RotorDriver drv(*m);
    REQUIRE(drv.numRotors() >= 1);

    const unsigned int n = drv.numAtoms();
    std::vector<double> orig = drv.positions();
    auto P = [](const std::vector<double> &b, unsigned int a) {
      return RDGeom::Point3D(b[3 * a], b[3 * a + 1], b[3 * a + 2]);
    };

    auto t = drv.torsion(0);
    unsigned int j = t[1], k = t[2];  // junction bond j(parent)-k(child)
    double origDih = drv.dihedralDeg(0);

    // child atoms = the k-side moving set plus k itself
    std::vector<char> isChild(n, 0);
    isChild[k] = 1;
    for (unsigned int a : drv.movingAtoms(0)) {
      isChild[a] = 1;
    }

    // the child in its own frame: copy of orig, then a fixed rigid displacement
    std::vector<RDGeom::Point3D> child(n);
    RDGeom::Transform3D disp;
    RDGeom::Point3D axis(0.3, -0.7, 0.5);
    axis.normalize();
    disp.SetRotation(1.1, axis);
    for (unsigned int a = 0; a < n; ++a) {
      RDGeom::Point3D p = P(orig, a);
      disp.TransformPoint(p);
      child[a] = p + RDGeom::Point3D(4.0, -2.0, 1.5);  // translate away too
    }

    // place the child: its bonding atom k joins the parent atom j; its exit
    // marker is j (the parent-side atom).  Use the true bond length so k lands
    // on orig[k].
    double bondLen = (P(orig, k) - P(orig, j)).length();
    placeChildCoords(child, /*childBondAtom=*/k, /*childExit=*/j,
                     /*parentBondAtom=*/P(orig, j), /*parentExit=*/P(orig, k),
                     bondLen);
    // the bonding atom must land on the original k position
    CHECK((child[k] - P(orig, k)).length() < 1e-6);

    // reassemble: parent atoms from orig, child atoms from the placed child
    std::vector<double> buf = orig;
    for (unsigned int a = 0; a < n; ++a) {
      if (isChild[a]) {
        buf[3 * a] = child[a].x;
        buf[3 * a + 1] = child[a].y;
        buf[3 * a + 2] = child[a].z;
      }
    }
    drv.positions() = buf;
    drv.drive(0, origDih);  // set the junction torsion back to the original

    // every atom should be back at its original position
    double maxDev = 0.0;
    for (unsigned int a = 0; a < n; ++a) {
      maxDev = std::max(maxDev, (P(drv.positions(), a) - P(orig, a)).length());
    }
    INFO("maxDev=" << maxDev);
    CHECK(maxDev < 1e-3);
  };
  reassemble("c1ccccc1-c1ccccc1");  // biphenyl
  reassemble("c1ccccc1CC");         // ethylbenzene
  reassemble("c1ccccc1OCC");        // aryl ether chain
}

TEST_CASE("RotorTree diversity-preserving retention keeps a spread ensemble",
          "[rotortree][diversity]") {
  auto coordRms = [](const std::vector<double> &a,
                     const std::vector<double> &b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
      double d = a[i] - b[i];
      s += d * d;
    }
    return std::sqrt(s / (a.size() / 3));
  };

  auto m =
      embed3D("OCCCCCCO");  // several rotors -> many within-window rotamers
  REQUIRE(m);

  // energy-greedy retention (divThr = 0)
  RotorDriver g(*m, -1, makeFullFFScoreFn(*m, false));
  RigidRotorSearchParams gp;
  gp.tree.beamWidth = 30;
  gp.diversityRmsThresh = 0.0;
  auto greedy = RotorTree(g, gp).search();

  // diversity-preserving retention
  RotorDriver d(*m, -1, makeFullFFScoreFn(*m, false));
  RigidRotorSearchParams dp;
  dp.tree.beamWidth = 30;
  dp.energyWindow = 20.0;
  dp.diversityRmsThresh = 0.6;
  auto diverse = RotorTree(d, dp).search();

  REQUIRE(diverse.size() >= 2);

  SECTION(
      "retention invariant: kept conformers are pairwise >= threshold apart") {
    for (size_t i = 0; i < diverse.size(); ++i) {
      for (size_t j = i + 1; j < diverse.size(); ++j) {
        CHECK(coordRms(diverse[i].coords, diverse[j].coords) >= 0.6 - 1e-6);
      }
    }
  }

  SECTION("diversity spreads the ensemble wider than energy-greedy") {
    // max spread from the best (lowest-energy) conformer
    auto maxSpread = [&](const std::vector<SearchResult> &r) {
      double mx = 0.0;
      for (size_t i = 1; i < r.size(); ++i) {
        mx = std::max(mx, coordRms(r[0].coords, r[i].coords));
      }
      return mx;
    };
    CHECK(maxSpread(diverse) >= maxSpread(greedy) - 1e-9);
    // both still report the same lowest energy first (diversity keeps the best
    // too)
    CHECK(diverse.front().score <= greedy.front().score + 1e-6);
  }
}

TEST_CASE("FragmentJoiner greedily reassembles a molecule coordinate-only",
          "[joiner][greedy]") {
  // Cut a real molecule into fragments (one conformer each = its crystal-frame
  // coords), then let the joiner place them and search the junctions.  The
  // assembled ensemble must reach the original geometry -- end-to-end proof of
  // place + drive + diversity search with the decomposed scorer.
  auto run = [](const std::string &smi, double thresh) {
    INFO("smiles=" << smi);
    auto m = embed3D(smi);
    REQUIRE(m);
    const unsigned int n = m->getNumAtoms();
    const Conformer &oc = m->getConformer();

    // rotatable (junction) bonds and the fragment partition (union-find over
    // the NON-rotatable bonds)
    auto links = FragmentConfGen::findLinkBonds(*m);
    REQUIRE(links.size() >= 1);
    std::vector<char> isLink(m->getNumBonds(), 0);
    for (auto b : links) isLink[b] = 1;
    std::vector<int> uf(n);
    for (unsigned int i = 0; i < n; ++i) uf[i] = i;
    std::function<int(int)> find = [&](int x) {
      while (uf[x] != x) {
        uf[x] = uf[uf[x]];
        x = uf[x];
      }
      return x;
    };
    for (const auto b : m->bonds()) {
      if (!isLink[b->getIdx()]) {
        uf[find(b->getBeginAtomIdx())] = find(b->getEndAtomIdx());
      }
    }
    std::map<int, unsigned int> fragId;
    std::vector<int> fragOf(n);
    for (unsigned int a = 0; a < n; ++a) {
      int r = find(a);
      if (!fragId.count(r)) fragId[r] = fragId.size();
      fragOf[a] = fragId[r];
    }
    unsigned int nFrag = fragId.size();

    // full-molecule-indexed coordinate copy (valid for every atom ->
    // owned+partners)
    std::vector<RDGeom::Point3D> pos(n);
    for (unsigned int a = 0; a < n; ++a) pos[a] = oc.getAtomPos(a);

    std::vector<JoinFragment> frags(nFrag);
    for (unsigned int a = 0; a < n; ++a) frags[fragOf[a]].atoms.push_back(a);
    for (auto &f : frags)
      f.confs.push_back({pos, 0.0});  // single conf = crystal frame

    std::vector<JoinJunction> junctions;
    for (auto bi : links) {
      const Bond *b = m->getBondWithIdx(bi);
      unsigned int a = b->getBeginAtomIdx(), c = b->getEndAtomIdx();
      JoinJunction J;
      J.fragA = fragOf[a];
      J.fragB = fragOf[c];
      J.atomA = a;
      J.atomB = c;
      J.bondLen = (oc.getAtomPos(a) - oc.getAtomPos(c)).length();
      junctions.push_back(J);
    }

    RigidRotorSearchParams p;
    p.rootSeeds = 1;
    p.tree.beamWidth = 80;
    p.energyWindow = 30.0;
    p.diversityRmsThresh = 0.4;
    // this checks RAW same-index reproduction of the input geometry, so keep
    // the raw prune (symmetry dedup would legitimately keep a symmetry-image
    // representative that has a large raw -- but ~0 symmetry-aware -- RMSD to
    // the original).
    p.finalSymmetryDedup = false;
    p.defaultAngles.clear();
    for (double a = -180.0; a < 180.0; a += 20.0)
      p.defaultAngles.push_back(a);

    const auto ctx = joinFragments(*m, fragOf, frags, junctions);
    REQUIRE(ctx.isValid());
    auto res = runRigidRotorSearch(ctx, p).results;
    REQUIRE(!res.empty());
    CHECK(res.front().coords.size() == 3u * n);

    // min heavy-atom RMSD to the original (buffer is in the crystal frame: root
    // is placed at its own crystal coords, so no global alignment is needed)
    double best = 1e9;
    for (const auto &r : res) {
      double sq = 0.0;
      unsigned int nh = 0;
      for (unsigned int a = 0; a < n; ++a) {
        if (m->getAtomWithIdx(a)->getAtomicNum() <= 1) continue;
        double dx = r.coords[3 * a] - oc.getAtomPos(a).x;
        double dy = r.coords[3 * a + 1] - oc.getAtomPos(a).y;
        double dz = r.coords[3 * a + 2] - oc.getAtomPos(a).z;
        sq += dx * dx + dy * dy + dz * dz;
        ++nh;
      }
      best = std::min(best, std::sqrt(sq / nh));
    }
    INFO("best heavy RMSD to original = " << best);
    CHECK(best < thresh);
  };
  run("c1ccccc1CCO", 0.8);  // 2 junctions
  run("c1ccccc1CCc1ccccc1",
      0.9);  // bibenzyl, 3 junctions, single-atom CH2 frags
  run("CCOC(=O)c1ccc(OC)cc1", 0.9);  // ester + methoxy, several junctions
}

TEST_CASE("buildFragmentJoinerInput + FragmentJoiner produce a valid conformer",
          "[joiner][integration]") {
  auto run = [](const std::string &smi) {
    INFO("smiles=" << smi);
    std::unique_ptr<ROMol> raw(SmilesToMol(smi));
    REQUIRE(raw);

    FragmentJoinerInput in =
        buildFragmentJoinerInput(*raw, /*nConfs=*/12, /*seed=*/0xC0FFEE);
    REQUIRE(in.mol);
    REQUIRE(in.fragments.size() >= 2);
    // every fragment embedded at least one conformer
    for (const auto &f : in.fragments) {
      REQUIRE(!f.confs.empty());
      // confs are energy-sorted
      for (size_t i = 1; i < f.confs.size(); ++i) {
        CHECK(f.confs[i - 1].energy <= f.confs[i].energy + 1e-6);
      }
    }

    RigidRotorSearchParams p;
    p.rootSeeds = 3;
    p.tree.beamWidth = 60;
    p.diversityRmsThresh = 0.5;
    auto res = runRigidRotorSearch(joinFragments(in), p).results;
    REQUIRE(!res.empty());

    const unsigned int n = in.mol->getNumAtoms();
    CHECK(res.front().coords.size() == 3u * n);
    for (double v : res.front().coords) CHECK(std::isfinite(v));

    // the best conformer must be a geometrically VALID conformer: every bond
    // length is physical (nothing broken or collapsed) -- proves place + drive
    // assembled real geometry, not garbage.
    auto P = [&](unsigned int a) {
      return RDGeom::Point3D(res.front().coords[3 * a],
                             res.front().coords[3 * a + 1],
                             res.front().coords[3 * a + 2]);
    };
    for (const auto b : in.mol->bonds()) {
      double d = (P(b->getBeginAtomIdx()) - P(b->getEndAtomIdx())).length();
      INFO("bond " << b->getBeginAtomIdx() << "-" << b->getEndAtomIdx()
                   << " len=" << d);
      CHECK(d > 0.8);
      CHECK(d < 1.9);
    }
  };
  run("c1ccccc1CCO");
  run("CCOC(=O)c1ccc(OC)cc1");
  run("c1ccccc1CCc1ccccc1");
}

TEST_CASE("every RigidRotorSearchMode runs and returns a valid ensemble",
          "[joiner][search]") {
  // The searches are extracted behind RigidRotorSearch; each mode must place +
  // drive real geometry.  Experimental modes currently delegate to tree, so
  // they must also succeed.
  std::unique_ptr<ROMol> raw(
      SmilesToMol("CCOC(=O)c1ccc(OCCN)cc1"));  // multi-rotor
  REQUIRE(raw);
  FragmentJoinerInput in =
      buildFragmentJoinerInput(*raw, /*nConfs=*/12, /*seed=*/0xC0FFEE);
  REQUIRE(in.mol);
  REQUIRE(in.fragments.size() >= 2);
  const unsigned int n = in.mol->getNumAtoms();

  for (auto mode :
       {RigidRotorSearchMode::Auto, RigidRotorSearchMode::Tree,
        RigidRotorSearchMode::Thompson, RigidRotorSearchMode::Systematic,
        RigidRotorSearchMode::Merged}) {
    INFO("searchMode=" << static_cast<int>(mode));
    RigidRotorSearchParams p;
    p.searchMode = mode;
    p.diversityRmsThresh = 0.5;
    auto res = runRigidRotorSearch(joinFragments(in), p).results;
    REQUIRE(!res.empty());
    CHECK(res.front().coords.size() == 3u * n);
    for (double v : res.front().coords) CHECK(std::isfinite(v));
  }
}

TEST_CASE("build() CoordsOnly assembly mode produces valid conformers",
          "[fragconfgen][coordsonly]") {
  auto heavyCanon = [](const ROMol &m) {
    RWMol h(m);
    MolOps::removeHs(h);
    return MolToSmiles(h);
  };
  for (const std::string &smi :
       {"c1ccccc1", "c1ccccc1CCO", "CCOC(=O)c1ccc(OC)cc1",
        "c1ccccc1CCc1ccccc1"}) {
    INFO("smiles=" << smi);
    std::unique_ptr<ROMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    const std::string want = heavyCanon(*mol);

    FragmentConfGenParams p;
    p.embedding.numConfsPerFragment = 10;
    p.numOutputConfs = 30;
    p.randomSeed = 0xf00d;
    auto res = FragmentConfGen(p).build(*mol);

    REQUIRE(!res.conformers.empty());
    CHECK(res.conformers.size() <= 30u);
    for (const auto &c : res.conformers) {
      // topology preserved (same molecule) and a finite MMFF energy tag present
      CHECK(heavyCanon(*c) == want);
      double e = 0.0;
      CHECK(c->getPropIfPresent<double>(kEnergyProp, e));
      CHECK(std::isfinite(e));
    }
    // best energy first
    for (size_t i = 1; i < res.conformers.size(); ++i) {
      double a = 0, b = 0;
      res.conformers[i - 1]->getPropIfPresent<double>(kEnergyProp, a);
      res.conformers[i]->getPropIfPresent<double>(kEnergyProp, b);
      CHECK(a <= b + 1e-6);
    }
    CHECK(std::isfinite(res.bestEnergy));
  }
}

TEST_CASE("createFragmentConfGen validates params and builds",
          "[fragconfgen][params][factory]") {
  SECTION(
      "default params are valid and the factory builds a working generator") {
    FragmentConfGenParams p;
    p.numOutputConfs = 20;
    p.randomSeed = 0xf00d;
    CHECK(p.validate().empty());
    auto gen = createFragmentConfGen(p);  // must not throw
    std::unique_ptr<ROMol> mol(SmilesToMol("c1ccccc1CCO"));
    REQUIRE(mol);
    auto res = gen.build(*mol);
    CHECK(!res.conformers.empty());
  }
  SECTION("incoherent params are rejected with a message") {
    FragmentConfGenParams bad;
    bad.joiner.ffVariant = "MMFF95";  // not a real variant
    CHECK_FALSE(bad.validate().empty());
    CHECK_THROWS_AS(createFragmentConfGen(bad), std::invalid_argument);

    FragmentConfGenParams bad2;
    bad2.energyWindow = -1.0;  // negative window
    CHECK_FALSE(bad2.validate().empty());

    FragmentConfGenParams bad3;
    bad3.search.thompson.minBudget = 5000;
    bad3.search.thompson.maxBudget = 300;  // min > max
    CHECK_FALSE(bad3.validate().empty());

    FragmentConfGenParams bad4;
    bad4.search.searchMode = RigidRotorSearchMode::Tree;
    bad4.search.tree.beamWidth = 0;  // must be > 0 for the search that runs
    CHECK_FALSE(bad4.validate().empty());

    // ... but under the DEFAULT Auto the same knob is not this run's problem:
    // thompson.autoBudget is on by default, so Auto's non-Systematic arm
    // resolves to Thompson and the tree beam is never read.
    FragmentConfGenParams unreached;
    REQUIRE(unreached.search.searchMode == RigidRotorSearchMode::Auto);
    REQUIRE(unreached.search.thompson.autoBudget);
    unreached.search.tree.beamWidth = 0;
    CHECK(unreached.validate().empty());
  }
  SECTION("Tree ignores thompsonBudget: a warning, not an error") {
    // RotorTree used to answer a non-zero budget with its OWN rotor-only
    // Thompson, so the combination may be set deliberately and now means
    // something else.  It stays VALID -- the run is well defined, the budget
    // is simply unused -- and warns once instead.
    FragmentConfGenParams p;
    p.search.searchMode = RigidRotorSearchMode::Tree;
    p.search.thompsonBudget = 500;
    CHECK(p.validate().empty());
    CHECK_NOTHROW(createFragmentConfGen(p));

    // Under Auto the same budget SELECTS Thompson rather than being ignored,
    // so there is nothing to warn about.
    FragmentConfGenParams autoP;
    REQUIRE(autoP.search.searchMode == RigidRotorSearchMode::Auto);
    autoP.search.thompsonBudget = 500;
    CHECK(autoP.validate().empty());
  }
  SECTION("label variant may be empty (means: use joiner variant)") {
    FragmentConfGenParams p;
    p.labelFFVariant = "";  // empty is allowed
    CHECK(p.validate().empty());
    p.labelFFVariant = "MMFF94s";  // explicit is allowed
    CHECK(p.validate().empty());
  }
  SECTION("the refined-torsion variants are accepted everywhere a variant is") {
    for (const std::string v :
         {"MMFF94", "MMFF94s", "MMFF94_TOR", "MMFF94s_TOR"}) {
      FragmentConfGenParams p;
      p.joiner.ffVariant = v;
      p.embedding.ffVariant = v;
      p.labelFFVariant = v;
      INFO("variant " << v);
      CHECK(p.validate().empty());
    }
    // a near-miss must be REJECTED rather than silently falling back to the
    // stock table -- that would run an experiment as its own control
    FragmentConfGenParams bad;
      bad.embedding.ffVariant = "MMFF94s_RTOR";  // a near-miss spelling
    CHECK_FALSE(bad.validate().empty());
  }

  SECTION("an incompatible fragment library is rejected at assignment") {
    // A library built with a different embedding holds fragments that do not
    // mean what this run expects.  This used to be a per-molecule throw from
    // inside build(); it is a property of the parameters, so it belongs to
    // validate() -- and setFraglib() applies it at the point of assignment.
    FragmentConfGenParams p;
    FraglibParams other = p.embedding;
    other.numConfsPerFragment = p.embedding.numConfsPerFragment + 7;
    auto wrongLib = std::make_shared<Fraglib>(other);

    CHECK_THROWS_AS(p.setFraglib(wrongLib), std::invalid_argument);
    CHECK(p.fraglib == nullptr);   // rejected: parameters left as they were
    CHECK(p.validate().empty());

    // assigning the field directly bypasses setFraglib, so validate() (and
    // therefore createFragmentConfGen) has to catch it too
    p.fraglib = wrongLib;
    CHECK_FALSE(p.validate().empty());
    CHECK_THROWS_AS(createFragmentConfGen(p), std::invalid_argument);

    // ... and a matching library is accepted
    FragmentConfGenParams q;
    auto goodLib = std::make_shared<Fraglib>(q.embedding);
    CHECK_NOTHROW(q.setFraglib(goodLib));
    CHECK(q.fraglib == goodLib);
    CHECK(q.validate().empty());
  }
}

TEST_CASE("joinFragments reports why a join failed", "[assembler][joiner]") {
  // The joiner used to do all of this in a constructor, which left it no way
  // to say what went wrong -- an unjoinable molecule just produced a context
  // whose coordinates were the all-zero placeholder.
  std::unique_ptr<ROMol> m(SmilesToMol("c1ccccc1CCc1ccccc1"));
  REQUIRE(m);
  auto in = buildFragmentJoinerInput(*m, /*nConfs=*/4, /*seed=*/0xf00d);
  REQUIRE(in.mol);
  REQUIRE(in.fragments.size() >= 2);

  SECTION("a well-formed input joins") {
    const auto ctx = joinFragments(in);
    CHECK(ctx.isValid());
    CHECK(ctx.status == FragmentJoinerStatus::Ok);
    CHECK(!ctx.edges.empty());
    CHECK(ctx.mol.getNumConformers() == 1);
    CHECK(ctx.scorer);
  }

  SECTION("no fragments") {
    auto empty = in;
    empty.fragments.clear();
    empty.junctions.clear();
    const auto ctx = joinFragments(empty);
    CHECK_FALSE(ctx.isValid());
    CHECK(ctx.status == FragmentJoinerStatus::NoFragments);
  }

  SECTION("unviable parameters are reported, not acted on") {
    FragmentJoinerParams bad;
    bad.ffVariant = "MMFF95";  // not a real variant
    CHECK_FALSE(bad.isValid());
    CHECK_FALSE(bad.validate().empty());

    const auto ctx = joinFragments(in, bad);
    CHECK_FALSE(ctx.isValid());
    CHECK(ctx.status == FragmentJoinerStatus::BadParams);
    // and it stopped BEFORE doing any of the work
    CHECK(ctx.edges.empty());
  }

  SECTION("a fragment with an empty conformer pool") {
    auto starved = in;
    starved.fragments.front().confs.clear();
    const auto ctx = joinFragments(starved);
    CHECK_FALSE(ctx.isValid());
    CHECK(ctx.status == FragmentJoinerStatus::NoFragmentConformers);
    // the rotor topology is still worked out -- only the placement is not
    CHECK(!ctx.rotorBonds.empty());
    // and a search on it declines rather than reading unplaced coordinates
    RigidRotorSearchParams sp;
    CHECK(runRigidRotorSearch(ctx, sp).results.empty());
  }
}

TEST_CASE("AUTO sentinels: 0 means DISABLED, AUTO means derive",
          "[fragconfgen][params][sentinel]") {
  // Regression guard for a bug that has bitten four times: a parameter whose
  // default 0 meant "derive a default" was read as "feature off".  The worst
  // case silently skipped the final symmetry dedup, so a search emitted
  // near-duplicates, saturated the output cap and scored worst of every mode --
  // while looking like an algorithmic result.
  SECTION("the derive-able thresholds default to AUTO, not to 0") {
    FragmentConfGenParams p;
    CHECK(isAuto(p.search.systematic.finalRms));
    CHECK(isAuto(p.search.systematic.nodeRms));
  }
  SECTION(
      "AUTO resolves to the derived value; 0 stays OFF; a real value is untouched") {
    CHECK(resolveAuto(kAutoD, 0.5) == Catch::Approx(0.5));
    CHECK(resolveAuto(0.0, 0.5) ==
          Catch::Approx(0.0));  // DISABLED must survive
    CHECK(resolveAuto(0.25, 0.5) == Catch::Approx(0.25));
    CHECK(resolveAuto(kAutoI, 1000) == 1000);
    CHECK(resolveAuto(0, 1000) == 0);
  }
  SECTION("AUTO (-1) round-trips through the serializer as a plain number") {
    FragmentConfGenParams p;
    p.search.systematic.finalRms = kAutoD;
    p.search.systematic.nodeRms = 0.0;  // explicitly DISABLED
    p.search.systematic.dedupEnergyBand = 0.25;  // explicit value
    const std::string text = fragmentConfGenParamsToString(p);
    CHECK(text.find("search.systematic.finalRms = -1") != std::string::npos);

    FragmentConfGenParams q;
    q.search.systematic.finalRms =
        9.0;  // must be overwritten by the AUTO token
    REQUIRE(fragmentConfGenParamsFromString(text, q).empty());
    CHECK(isAuto(q.search.systematic.finalRms));
    CHECK(q.search.systematic.nodeRms == Catch::Approx(0.0));
    CHECK(q.search.systematic.dedupEnergyBand == Catch::Approx(0.25));
  }
  SECTION("the resolved dump reports EFFECTIVE values, not the raw struct") {
    FragmentConfGenParams p;
    p.joiner.ffVariant = "MMFF94s";
    p.embedding.ffVariant =
        "MMFF94";  // the pipeline overrides this from the joiner
    p.energyWindow = 10.0;
    const FragmentConfGenParams r = resolvedFragmentConfGenParams(p);
    CHECK(r.embedding.ffVariant == "MMFF94s");  // driven, not as authored
    CHECK(r.search.energyWindow == Catch::Approx(25.0));  // 2.5x, not 10
    CHECK(r.search.systematic.finalRms ==
          Catch::Approx(0.5));  // AUTO -> derived
    CHECK_FALSE(isAuto(r.search.systematic.nodeRms));
    // and an explicitly DISABLED value must NOT be silently turned on
    FragmentConfGenParams off = p;
    off.search.systematic.finalRms = 0.0;
    CHECK(resolvedFragmentConfGenParams(off).search.systematic.finalRms ==
          Catch::Approx(0.0));
  }
}

TEST_CASE(
    "MMFF94(s)_TOR applies the refined dihedral parameters of Wahl et al.",
    "[fragconfgen][params][mmff]") {
  // J. Wahl, J. Freyss, M. von Korff, T. Sander, J. Cheminform. 2019, 11, 53.
  // Two changed entries; everything outside that chemistry must be untouched.
  struct Probe {
    const char *smiles;
    unsigned int typeJ, typeK;  // the bond whose torsion is refined
    double stockV2, refinedV2;
  };
  const std::vector<Probe> probes = {
      {"CC(=O)Nc1ccccc1", 10, 37, 6.0, 2.7},     // N-aryl amide
      {"c1ccc(-n2cccc2)cc1", 37, 39, 6.0, 2.6},  // phenylpyrrole
  };

  for (const auto &probe : probes) {
    INFO("smiles " << probe.smiles);
    std::unique_ptr<RWMol> mol(SmilesToMol(probe.smiles));
    REQUIRE(mol);
    MolOps::addHs(*mol);
    REQUIRE(DGeomHelpers::EmbedMolecule(*mol, 0, 0xf00d) >= 0);

    MMFF::MMFFMolProperties stock(*mol, "MMFF94s");
    MMFF::MMFFMolProperties refined(*mol, "MMFF94s_TOR");
    REQUIRE(stock.isValid());
    REQUIRE(refined.isValid());
    CHECK(stock.getMMFFVariant() == "MMFF94s");
    CHECK(refined.getMMFFVariant() == "MMFF94s_TOR");

    unsigned int nSeen = 0;
    for (const auto bond : mol->bonds()) {
      const unsigned int j = bond->getBeginAtomIdx(), k = bond->getEndAtomIdx();
      const unsigned int tj = stock.getMMFFAtomType(j),
                         tk = stock.getMMFFAtomType(k);
      if (!((tj == probe.typeJ && tk == probe.typeK) ||
            (tj == probe.typeK && tk == probe.typeJ))) {
        continue;
      }
      for (const auto nbI : mol->atomNeighbors(mol->getAtomWithIdx(j))) {
        if (nbI->getIdx() == k) continue;
        for (const auto nbL : mol->atomNeighbors(mol->getAtomWithIdx(k))) {
          if (nbL->getIdx() == j) continue;
          ForceFields::MMFF::MMFFTor a, b;
          unsigned int ta = 0, tb = 0;
          if (!stock.getMMFFTorsionParams(*mol, nbI->getIdx(), j, k,
                                          nbL->getIdx(), ta, a) ||
              !refined.getMMFFTorsionParams(*mol, nbI->getIdx(), j, k,
                                            nbL->getIdx(), tb, b)) {
            continue;
          }
          CHECK(a.V2 == Catch::Approx(probe.stockV2));
          CHECK(b.V2 == Catch::Approx(probe.refinedV2));
          ++nSeen;
        }
      }
    }
    CHECK(nSeen > 0);  // the motif must actually have been found
  }

  SECTION("a molecule without the motif is bit-identical under both") {
    std::unique_ptr<RWMol> mol(SmilesToMol("CCCCCC"));
    REQUIRE(mol);
    MolOps::addHs(*mol);
    REQUIRE(DGeomHelpers::EmbedMolecule(*mol, 0, 0xf00d) >= 0);
    MMFF::MMFFMolProperties stock(*mol, "MMFF94s");
    MMFF::MMFFMolProperties refined(*mol, "MMFF94s_TOR");
    REQUIRE(stock.isValid());
    REQUIRE(refined.isValid());
    std::unique_ptr<ForceFields::ForceField> ffA(
        MMFF::constructForceField(*mol, &stock));
    std::unique_ptr<ForceFields::ForceField> ffB(
        MMFF::constructForceField(*mol, &refined));
    ffA->initialize();
    ffB->initialize();
    CHECK(ffA->calcEnergy() == ffB->calcEnergy());
  }
}

TEST_CASE("FragmentConfGenParams text round-trips through the serializer",
          "[fragconfgen][params][paramsio]") {
  SECTION("mutated non-default params survive a write -> read cycle") {
    FragmentConfGenParams p;
    // touch a field of every serialized type / sub-struct, all off-default
    p.numOutputConfs = 37;
    p.energyWindow = 12.5;
    p.randomSeed = 99;
    p.embedding.fragmentEmbedMode = FragmentEmbedMode::DG;
    p.embedding.minimizeMode = FragmentMinimize::ShrugScore;
    p.sampleTrivialRotors = true;
    p.embedding.numConfsPerFragment = 250;
    p.search.searchMode = RigidRotorSearchMode::Systematic;
    p.joiner.ffVariant = "MMFF94s";
    p.joiner.fragShrugDisplacement = 0.15;
    p.outputRanking = OutputRanking::InterFragScore;
    p.rankByBasinEnergy = true;
    p.labelFFVariant = "MMFF94s";
    p.search.thompson.priorStrength = 4.5;
    p.search.thompson.maxConfs = 25;
    p.search.thompson.outMode = OutputSelection::Diverse;
    p.search.systematic.eWindow = 18.0;
    p.search.systematic.maxPoolConfs = 9000;
    p.search.systematic.maxFragmentConfs = 24;
    p.diagnostics.ASM_PROFILE = true;
    p.diagnostics.SYS_VALIDATE = true;

    const std::string text = fragmentConfGenParamsToString(p);

    FragmentConfGenParams q;  // fresh defaults
    const std::string err = fragmentConfGenParamsFromString(text, q);
    CHECK(err.empty());

    CHECK(q.numOutputConfs == 37);
    CHECK(q.energyWindow == Catch::Approx(12.5));
    CHECK(q.randomSeed == 99);
    CHECK(q.embedding.fragmentEmbedMode == FragmentEmbedMode::DG);
    CHECK(q.embedding.minimizeMode == FragmentMinimize::ShrugScore);
    CHECK(q.sampleTrivialRotors == true);
    CHECK(q.embedding.numConfsPerFragment == 250u);
    CHECK(q.search.searchMode == RigidRotorSearchMode::Systematic);
    CHECK(q.joiner.ffVariant == "MMFF94s");
    CHECK(q.joiner.fragShrugDisplacement == Catch::Approx(0.15));
    CHECK(q.outputRanking == OutputRanking::InterFragScore);
    CHECK(q.rankByBasinEnergy == true);
    CHECK(q.labelFFVariant == "MMFF94s");
    CHECK(q.search.thompson.priorStrength == Catch::Approx(4.5));
    CHECK(q.search.thompson.maxConfs == 25u);
    CHECK(q.search.thompson.outMode == OutputSelection::Diverse);
    CHECK(q.search.systematic.eWindow == Catch::Approx(18.0));
    CHECK(q.search.systematic.maxPoolConfs == 9000);
    CHECK(q.search.systematic.maxFragmentConfs == 24);
    CHECK(q.diagnostics.ASM_PROFILE == true);
    CHECK(q.diagnostics.SYS_VALIDATE == true);
  }
  SECTION("a partial file overrides only the keys it names") {
    FragmentConfGenParams p;  // defaults
    const int origOut = p.numOutputConfs;
    const std::string err = fragmentConfGenParamsFromString(
        "# just one knob\njoiner.ffVariant = MMFF94s\n", p);
    CHECK(err.empty());
    CHECK(p.joiner.ffVariant == "MMFF94s");
    CHECK(p.numOutputConfs == origOut);  // untouched
  }
  SECTION("unknown key and bad value are reported, not silently ignored") {
    FragmentConfGenParams p;
    CHECK_FALSE(fragmentConfGenParamsFromString("noSuchKnob = 3\n", p).empty());
    CHECK_FALSE(fragmentConfGenParamsFromString("numOutputConfs = banana\n", p)
                    .empty());
    CHECK_FALSE(fragmentConfGenParamsFromString("search.searchMode = Nope\n", p)
                    .empty());
  }
}

TEST_CASE("assembled junction angles are physically sane",
          "[fragconfgen][angles]") {
  // Reassembling cut fragments by exit-vector alignment must not produce
  // collapsed or splayed junctions: every heavy-atom bond angle should sit near
  // its ideal (sp3 ~109.5, aromatic ~120).  This guards the assembly geometry
  // across the embedding / torsion / search sweeps -- a mis-aligned exit vector
  // shows up as an out-of-band angle here.  benzene+CCN exercises an aromatic-C
  // junction; an all-sp3 chain+N exercises aliphatic junctions.
  auto angleAt = [](const Conformer &conf, unsigned int i, unsigned int j,
                    unsigned int k) {
    RDGeom::Point3D a = conf.getAtomPos(i) - conf.getAtomPos(j);
    RDGeom::Point3D b = conf.getAtomPos(k) - conf.getAtomPos(j);
    double denom = a.length() * b.length();
    double c = denom > 0 ? a.dotProduct(b) / denom : 0.0;
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c) * 180.0 / M_PI;
  };
  for (const std::string &smi : {"c1ccccc1CCN", "CCCCCCCCN"}) {
    INFO("smiles=" << smi);
    std::unique_ptr<RWMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    FragmentConfGenParams p;
    p.embedding.numConfsPerFragment = 10;
    p.numOutputConfs = 20;
    p.randomSeed = 0xf00d;
    auto res = FragmentConfGen(p).build(*mol);
    REQUIRE(!res.conformers.empty());
    const ROMol &best = *res.conformers.front();
    const Conformer &conf = best.getConformer();
    for (const auto atom : best.atoms()) {
      if (atom->getDegree() < 2) {
        continue;
      }
      std::vector<unsigned int> nbrs;
      for (const auto nb : best.atomNeighbors(atom)) {
        nbrs.push_back(nb->getIdx());
      }
      for (size_t x = 0; x + 1 < nbrs.size(); ++x) {
        for (size_t y = x + 1; y < nbrs.size(); ++y) {
          double ang = angleAt(conf, nbrs[x], atom->getIdx(), nbrs[y]);
          INFO("angle at atom " << atom->getIdx() << " = " << ang);
          CHECK(ang > 100.0);  // no collapsed junction
          CHECK(ang < 126.0);  // no splayed/near-linear junction
        }
      }
    }
  }
}

TEST_CASE("Theobald QCP RMSD matches optimal-superposition alignment",
          "[qcprmsd]") {
  using namespace RDKit;

  // A real 3D conformer to work from.
  std::unique_ptr<RWMol> mol(SmilesToMol("CC(=O)Nc1ccc(O)cc1CCN"));
  REQUIRE(mol);
  MolOps::addHs(*mol);
  DGeomHelpers::EmbedParameters ps(DGeomHelpers::ETKDGv3);
  ps.randomSeed = 0xC0FFEE;
  REQUIRE(DGeomHelpers::EmbedMolecule(*mol, ps) == 0);

  const unsigned int n = mol->getNumAtoms();
  const Conformer &conf = mol->getConformer();
  std::vector<RDGeom::Point3D> a(n);
  for (unsigned int i = 0; i < n; ++i) a[i] = conf.getAtomPos(i);
  std::vector<unsigned int> idx(n);
  std::iota(idx.begin(), idx.end(), 0u);

  SECTION("identical structures give zero") {
    CHECK(qcpRmsd(a, a, idx) == Catch::Approx(0.0).margin(1e-5));
  }

  SECTION("invariant under rotation + translation") {
    RDGeom::Transform3D t;
    RDGeom::Point3D axis(0.3, -0.8, 0.5);
    axis.normalize();  // SetRotation needs a unit axis for a pure rotation
    t.SetRotation(0.7, axis);
    std::vector<RDGeom::Point3D> b(n);
    for (unsigned int i = 0; i < n; ++i) {
      RDGeom::Point3D p = a[i];
      t.TransformPoint(p);
      b[i] = p + RDGeom::Point3D(4.0, -2.0, 7.0);  // rigid translation
    }
    CHECK(qcpRmsd(a, b, idx) == Catch::Approx(0.0).margin(1e-4));
  }

  SECTION("matches RDKit MolAlign optimal RMSD on a perturbed structure") {
    // second conformer = a with a deterministic per-atom perturbation
    auto *c2 = new Conformer(n);
    std::vector<RDGeom::Point3D> b(n);
    for (unsigned int i = 0; i < n; ++i) {
      RDGeom::Point3D p =
          a[i] + RDGeom::Point3D(0.11 * std::sin(0.7 * i),
                                 0.09 * std::cos(1.3 * i),
                                 0.07 * std::sin(2.1 * i + 1.0));
      b[i] = p;
      c2->setAtomPos(i, p);
    }
    const unsigned int cid2 = mol->addConformer(c2, true);

    // ground truth: RDKit's Kabsch/SVD optimal-superposition RMSD (prb=cid2,
    // ref=0)
    RDGeom::Transform3D trans;
    double ref = MolAlign::getAlignmentTransform(*mol, *mol, trans, cid2, 0);

    double got = qcpRmsd(a, b, idx);
    CHECK(got == Catch::Approx(ref).epsilon(1e-3));
    CHECK(got > 0.05);  // sanity: the perturbation is real, not ~0
  }

  SECTION("subset of atoms (heavy only) also matches alignment") {
    std::vector<unsigned int> heavy;
    for (unsigned int i = 0; i < n; ++i) {
      if (mol->getAtomWithIdx(i)->getAtomicNum() > 1) heavy.push_back(i);
    }
    // perturb, build conformer, compare over the heavy-atom map
    auto *c2 = new Conformer(n);
    std::vector<RDGeom::Point3D> b(n);
    for (unsigned int i = 0; i < n; ++i) {
      RDGeom::Point3D p = a[i] + RDGeom::Point3D(0.05 * std::cos(0.9 * i), 0.0,
                                                 0.06 * std::sin(1.7 * i));
      b[i] = p;
      c2->setAtomPos(i, p);
    }
    const unsigned int cid2 = mol->addConformer(c2, true);
    MatchVectType heavyMap;
    for (unsigned int h : heavy) heavyMap.emplace_back(h, h);
    RDGeom::Transform3D trans;
    double ref =
        MolAlign::getAlignmentTransform(*mol, *mol, trans, cid2, 0, &heavyMap);
    double got = qcpRmsd(a, b, heavy);
    CHECK(got == Catch::Approx(ref).epsilon(1e-3));
  }
}

TEST_CASE("symmetry-aware RMSD dedup pool", "[symrmsd]") {
  using namespace RDKit;

  SECTION("automorphism counts match the graph symmetry") {
    // asymmetric drug-like molecule -> only the identity
    std::unique_ptr<RWMol> asym(SmilesToMol("CC(=O)Nc1ccc(O)cc1CCN"));
    std::vector<unsigned int> h;
    CHECK(heavyAtomAutomorphisms(*asym, h).size() == 1u);

    // benzene: heavy graph is a 6-cycle -> dihedral group D6, order 12
    std::unique_ptr<RWMol> bz(SmilesToMol("c1ccccc1"));
    auto bzAutos = heavyAtomAutomorphisms(*bz, h);
    CHECK(h.size() == 6u);
    CHECK(bzAutos.size() == 12u);

    // p-dichlorobenzene: 2-fold + ring flip
    std::unique_ptr<RWMol> px(SmilesToMol("Clc1ccc(Cl)cc1"));
    CHECK(heavyAtomAutomorphisms(*px, h).size() > 1u);
  }

  SECTION("symmetry image is a duplicate; plain RMSD is not") {
    std::unique_ptr<RWMol> mol(SmilesToMol("O=C(c1ccccc1)c1ccccc1"));  // 2-fold
    REQUIRE(mol);
    MolOps::addHs(*mol);
    DGeomHelpers::EmbedParameters ps(DGeomHelpers::ETKDGv3);
    ps.randomSeed = 0xBEEF;
    REQUIRE(DGeomHelpers::EmbedMolecule(*mol, ps) == 0);

    std::vector<unsigned int> heavy;
    auto autos = heavyAtomAutomorphisms(*mol, heavy);
    REQUIRE(autos.size() >= 2u);  // has real topological symmetry

    const unsigned int n = mol->getNumAtoms();
    const Conformer &conf = mol->getConformer();
    std::vector<double> A(3 * n);
    for (unsigned int i = 0; i < n; ++i) {
      const auto p = conf.getAtomPos(i);
      A[3 * i] = p.x;
      A[3 * i + 1] = p.y;
      A[3 * i + 2] = p.z;
    }
    // B = A with a NON-identity automorphism applied to the heavy atoms
    const auto &sigma = autos[1];
    std::vector<double> B = A;
    for (unsigned int p = 0; p < heavy.size(); ++p) {
      const unsigned int dst = heavy[sigma[p]], src = heavy[p];
      B[3 * dst] = A[3 * src];
      B[3 * dst + 1] = A[3 * src + 1];
      B[3 * dst + 2] = A[3 * src + 2];
    }

    // plain (identity-only) heavy QCP sees a real displacement...
    const double plain = qcpRmsd(A, B, heavy, heavy);
    CHECK(plain > 0.1);
    // ...but the symmetry-aware minimum recognises B as the same structure
    RMSDPruner acc(*mol, 0.25);
    REQUIRE(acc.add(A));  // first conformer kept
    CHECK(acc.RMSD(B) == Catch::Approx(0.0).margin(1e-4));
    CHECK_FALSE(acc.add(B));  // rejected: symmetry duplicate
    CHECK(acc.size() == 1u);

    // a genuinely different geometry IS kept
    std::vector<double> C = A;
    for (unsigned int hIdx : heavy) {
      C[3 * hIdx] += 1.5 * std::sin(0.9 * hIdx);
      C[3 * hIdx + 2] += 1.2 * std::cos(1.4 * hIdx);
    }
    CHECK(acc.add(C));
    CHECK(acc.size() == 2u);
  }
}

TEST_CASE("search parameters are validated by the search that will run",
          "[assembler][search][params]") {
  // Each search owns the rules for its own knobs, so FragmentConfGenParams
  // does not have to restate them and let the two drift.
  SECTION("the base rejects what every search reads") {
    RigidRotorSearchParams sp;
    sp.defaultAngles.clear();  // no fallback arms -> no rotor has candidates
    CHECK_FALSE(makeRigidRotorSearch(RigidRotorSearchMode::Tree)
                    ->isValid(sp, "MMFF94"));

    RigidRotorSearchParams neg;
    neg.energyWindow = -1.0;
    CHECK_FALSE(makeRigidRotorSearch(RigidRotorSearchMode::Thompson)
                    ->isValid(neg, "MMFF94"));
  }

  SECTION("a search only answers for its own knobs") {
    RigidRotorSearchParams sp;
    sp.tree.beamWidth = 0;  // Tree's problem, nobody else's
    CHECK_FALSE(makeRigidRotorSearch(RigidRotorSearchMode::Tree)
                    ->isValid(sp, "MMFF94"));
    CHECK(makeRigidRotorSearch(RigidRotorSearchMode::Thompson)
              ->isValid(sp, "MMFF94"));

    RigidRotorSearchParams bp;
    bp.thompson.minBudget = bp.thompson.maxBudget + 1;
    CHECK_FALSE(makeRigidRotorSearch(RigidRotorSearchMode::Thompson)
                    ->isValid(bp, "MMFF94"));
    CHECK(makeRigidRotorSearch(RigidRotorSearchMode::Tree)
              ->isValid(bp, "MMFF94"));
  }

  SECTION("the searches that build their own MMFF terms declare it") {
    // Inert today -- every variant isValidFF() accepts is an MMFF variant --
    // but this is the hook that starts biting when a UFF family arrives.
    CHECK(makeRigidRotorSearch(RigidRotorSearchMode::Systematic)
              ->checkFF("MMFF94s"));
    CHECK_FALSE(makeRigidRotorSearch(RigidRotorSearchMode::Tree)
                    ->checkFF("MMFF94s"));
    CHECK_FALSE(makeRigidRotorSearch(RigidRotorSearchMode::Thompson)
                    ->checkFF("MMFF94s"));
    // an MMFF-only search rejects a force field it cannot build
    RigidRotorSearchParams sp;
    CHECK_FALSE(makeRigidRotorSearch(RigidRotorSearchMode::Systematic)
                    ->isValid(sp, "UFF"));
    CHECK(makeRigidRotorSearch(RigidRotorSearchMode::Systematic)
              ->isValid(sp, "MMFF94s"));
  }

  SECTION("Auto is checked against only the searches it can still pick") {
    // The rotor count decides Systematic vs the Thompson/Tree arm, but WHICH of
    // those two the arm takes is fixed by the budget parameters -- so Auto is
    // never all four, and a knob for an unreachable search must not fail a run.
    RigidRotorSearchParams sp;
    sp.searchMode = RigidRotorSearchMode::Auto;
    sp.thompson.autoBudget = true;  // the arm resolves to Thompson, not Tree
    sp.autoSystematicMinRotors = 12;
    auto modes = reachableSearchModes(sp);
    CHECK(modes.size() == 2);
    CHECK(std::find(modes.begin(), modes.end(),
                    RigidRotorSearchMode::Systematic) != modes.end());
    CHECK(std::find(modes.begin(), modes.end(),
                    RigidRotorSearchMode::Thompson) != modes.end());
    CHECK(std::find(modes.begin(), modes.end(), RigidRotorSearchMode::Tree) ==
          modes.end());
    CHECK(validateSearchParams(sp, "MMFF94").empty());

    // Tree cannot be reached, so its broken beam width is not this run's problem
    sp.tree.beamWidth = 0;
    CHECK(validateSearchParams(sp, "MMFF94").empty());

    // Systematic CAN be reached, so its broken pool is
    sp.systematic.maxPoolConfs = 0;
    CHECK_FALSE(validateSearchParams(sp, "MMFF94").empty());

    // switch Systematic off and it stops being checked
    sp.autoSystematicMinRotors = 0;
    CHECK(reachableSearchModes(sp).size() == 1);
    CHECK(validateSearchParams(sp, "MMFF94").empty());

    // ... and naming one search checks exactly that search
    sp.searchMode = RigidRotorSearchMode::Tree;
    CHECK(reachableSearchModes(sp) ==
          std::vector<RigidRotorSearchMode>{RigidRotorSearchMode::Tree});
    CHECK_FALSE(validateSearchParams(sp, "MMFF94").empty());  // beamWidth == 0
  }

  SECTION("FragmentConfGenParams delegates to both") {
    FragmentConfGenParams p;
    CHECK(p.validate().empty());

    FragmentConfGenParams badJoiner;
    badJoiner.joiner.interFragVdwCutoff = -1.0;
    CHECK_FALSE(badJoiner.validate().empty());
    CHECK_THROWS_AS(createFragmentConfGen(badJoiner), std::invalid_argument);

    FragmentConfGenParams badSearch;
    badSearch.search.defaultAngles.clear();
    CHECK_FALSE(badSearch.validate().empty());
    CHECK_THROWS_AS(createFragmentConfGen(badSearch), std::invalid_argument);
  }
}
