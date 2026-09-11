//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Embedder/Fraglib.h"

#include <GraphMol/MolOps.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <ForceField/ForceField.h>
#include <ForceField/PositionConstraint.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include "Utils/SymmetricRmsd.h"
#include "Joiner/JoinerProfiling.h"
#include <Geometry/Transform3D.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/types.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <map>
#include <chrono>
#include <cstdio>
#include <memory>
#include <set>
#include <stdexcept>
#include <atomic>
#include <vector>

namespace RDKit {

FragmentParams getDefaultFragmentParams(FragmentClass cls) {
  switch (cls) {
    // Rigid = aromatic / fused / <=4-ring: essentially ONE conformer, so the
    // pool is small. See docs/fragment-embedding.md for the pool-size response
    // of each class.
    case FragmentClass::Rigid:
      return {30, 30, 100, 8., 0.1, FragmentSamplesOperator::MAXIMUM, 0};
    case FragmentClass::SmallRing:
      return {30, 800, 800, 8., 0.1, FragmentSamplesOperator::MULTIPLY, 6};
    case FragmentClass::LargeRing:
      return {100, 1600, 1000, 24., 0.1, FragmentSamplesOperator::EXPONENTIAL,
              10};
    case FragmentClass::Acyclic:
      return {5, 50, 1, 8., 0.1, FragmentSamplesOperator::EXPONENTIAL, 10};
    // Non-per-class fallbacks: a single generic sizing when perClassEmbedding
    // is off.
    case FragmentClass::Exhaustive:
      return {30, 1000, 100, 12., 0.1, FragmentSamplesOperator::MAXIMUM, 0};
    case FragmentClass::Fast:
      return {10, 100, 20, 8., 0.1, FragmentSamplesOperator::MAXIMUM, 0};
  };
  return {30, 800, 800, 8., 0.1, FragmentSamplesOperator::MAXIMUM,
          0};  // unreachable
}

std::map<FragmentClass, FragmentParams> getDefaultFragmentParams() {
  std::map<FragmentClass, FragmentParams> m;
  for (auto c : {FragmentClass::Rigid, FragmentClass::SmallRing,
                 FragmentClass::LargeRing, FragmentClass::Acyclic,
                 FragmentClass::Exhaustive, FragmentClass::Fast}) {
    m[c] = getDefaultFragmentParams(c);
  }
  return m;
}

namespace {

//! Classify by ring content.
//!  This is designed to sample puckers/n-inversions etc.
//!  Acyclic = No Rings
//!  Rigid: Aromatic or fused/bridge/spiro
//!  LargeRing > 10 atoms
//!  SmallRing otherwise
FragmentClass classifyFragment(RWMol &frag) {
  if (!frag.getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(frag);
  }
  const RingInfo *ri = frag.getRingInfo();
  if (ri->numRings() == 0) {
    return FragmentClass::Acyclic;
  }
  std::size_t maxRing = 0;
  bool allAromatic = true;
  bool fused = false;
  for (const auto &ring : ri->atomRings()) {
    maxRing = std::max(maxRing, ring.size());
    for (int idx : ring) {
      const auto u = static_cast<unsigned int>(idx);
      if (!frag.getAtomWithIdx(u)->getIsAromatic()) {
        allAromatic = false;
      }
      if (ri->numAtomRings(u) >= 2) {
        fused = true;  // fused / bridged / spiro
      }
    }
  }
  if (allAromatic || maxRing <= 4 || fused) {
    return FragmentClass::Rigid;
  }
  return maxRing >= 10 ? FragmentClass::LargeRing : FragmentClass::SmallRing;
}

const char *fragmentClassName(FragmentClass cls) {
  switch (cls) {
    case FragmentClass::Rigid:
      return "Rigid";
    case FragmentClass::SmallRing:
      return "SmallRing";
    case FragmentClass::LargeRing:
      return "LargeRing";
    case FragmentClass::Exhaustive:
      return "Exhaustive";
    case FragmentClass::Fast:
      return "Fast";
    case FragmentClass::Acyclic:
    default:
      return "Acyclic";
  }
}

//! Embed a fragment
bool embedFragmentInPlace(RWMol &frag, const FraglibParams &params) {
  RDLog::LogStateSetter blocker;
  const long long tEmbed0 = detail::profiling() ? detail::nowNs() : 0;

  if (!frag.getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(frag);
  }
  // if we are a trivial rotor (CF3) do nothing.
  const RingInfo *ri = frag.getRingInfo();
  auto nonHDeg = [&](const Atom *a) {
    unsigned int d = 0;
    for (const auto nbr : frag.atomNeighbors(a)) {
      if (nbr->getAtomicNum() != 1) {
        ++d;
      }
    }
    return d;
  };
  unsigned int rot = 0;
  for (const auto b : frag.bonds()) {
    if (b->getBondType() != Bond::SINGLE || ri->numBondRings(b->getIdx()) > 0) {
      continue;
    }
    const Atom *a1 = b->getBeginAtom(), *a2 = b->getEndAtom();
    if (a1->getAtomicNum() == 0 || a2->getAtomicNum() == 0) {
      continue;  // junctoin bond, ignore
    }
    if (nonHDeg(a1) < 2 || nonHDeg(a2) < 2) {
      continue;  // another trivial rotor
    }
    ++rot;
  }

  const FragmentClass fc = params.perClassEmbedding ? classifyFragment(frag)
                                                    : FragmentClass::Exhaustive;
  const FragmentParams fp = params.classParams.at(fc);
  unsigned int keepN = static_cast<unsigned int>(fp.maxConfs);
  double poolRms = fp.rmsd;
  double eWindow = fp.eWindow;
  unsigned int poolSize =
      static_cast<unsigned int>(fp.computeMaxSamples(static_cast<int>(rot)));
  // Do we minimize or shrug off non MMFF energies?  The class recipe may
  // override the global embed/minimise settings (e.g. floppy Acyclic chains via
  // DG+MMFF while rings keep ETKDG).
  const FragmentMinimize minMode =
      fp.minimizeMode.value_or(params.minimizeMode);
  const unsigned int minIters =
      resolveAuto(fp.minimizeMaxIters, params.minimizeMaxIters);
  bool fullMin = minMode == FragmentMinimize::Full;
  bool scoreKeep = minMode == FragmentMinimize::Score ||
                   minMode == FragmentMinimize::ShrugScore;
  double scoreShrug =
      minMode == FragmentMinimize::ShrugScore ? params.shrugDisplacement : 0.0;

  FragmentEmbedMode embedMode = fp.embedMode.value_or(params.fragmentEmbedMode);

  // Junction bonds are capped with dummies, so for embedding we replace them
  //  with carbon to get decent output geometries.  It is up to the Joiner
  //  to convert them to ideal bond lengths when the full atom types are known
  const int CARBON = 6;
  std::vector<unsigned int> dummies;
  for (auto atom : frag.atoms()) {
    if (atom->getAtomicNum() == 0) {
      dummies.push_back(atom->getIdx());
      atom->setAtomicNum(CARBON);
      atom->setNoImplicit(CARBON == 1);
      atom->updatePropertyCache(false);
    }
  }

  MolOps::addHs(frag);

  // if we added Hs to the dummies we will need to remove them later
  std::vector<unsigned int> capHs;
  if (!dummies.empty()) {
    std::set<unsigned int> dset(dummies.begin(), dummies.end());
    for (const auto atom : frag.atoms()) {
      if (atom->getAtomicNum() != 1) {
        continue;
      }
      for (const auto nbr : frag.atomNeighbors(atom)) {
        if (dset.count(nbr->getIdx())) {
          capHs.push_back(atom->getIdx());
          break;
        }
      }
    }
  }
  DGeomHelpers::EmbedParameters ps(embedMode == FragmentEmbedMode::ETKDG
                                       ? DGeomHelpers::ETKDGv3
                                       : DGeomHelpers::EmbedParameters());
  ps.randomSeed = params.randomSeed;

  // Embed the full sample pool, then reduce to keepN by RMSD-diverse,
  // lowest-energy selection.
  const unsigned int nEmbed = poolSize;
  INT_VECT cids;
  try {
    DGeomHelpers::EmbedMultipleConfs(frag, cids, nEmbed, ps);
  } catch (...) {
    return false;
  }
  if (cids.empty()) {
    return false;  // embedding failed
  }

  // Compute final energies for sorting fragment conformations
  //  Full:        full minimise, MUTATE coords (DG; raw geometry must be
  //  relaxed) Score:       minimise for the energy, RESTORE embed coords
  //  (ETKDG) ShrugScore:  heavy flat-bottom ("shrug") minimise for the energy,
  //  restore coords None:        single-point at the embed coords
  std::vector<std::pair<double, int>> byEnergy;  // (energy, cid)

  // XXX FIX ME: this only honors MMFF94 variants right now.
  //   UFF and others will need to be a fallback or other variant or factory
  if (!cids.empty() && (fullMin || scoreKeep || cids.size() > 1)) {
    MMFF::MMFFMolProperties props(frag, params.ffVariant);
    if (props.isValid()) {
      props.setMMFFEleTerm(false);
      const unsigned int nA = frag.getNumAtoms();
      if (fullMin || scoreKeep) {
        // per-conformer force field (minimisation is conf-bound)
        for (const auto cid : cids) {
          std::unique_ptr<ForceFields::ForceField> ff(
              MMFF::constructForceField(frag, &props, 1.0e8, cid));
          if (!ff) continue;
          Conformer &c = frag.getConformer(static_cast<unsigned int>(cid));
          std::vector<RDGeom::Point3D> saved;
          if (scoreKeep) {  // remember embed coords; restore after scoring
            saved.resize(nA);
            for (unsigned int i = 0; i < nA; ++i) saved[i] = c.getAtomPos(i);
            if (scoreShrug > 0.0) {
              for (unsigned int i = 0; i < nA; ++i)
                if (frag.getAtomWithIdx(i)->getAtomicNum() != 1) {
                  ff->contribs().push_back(ForceFields::ContribPtr(
                      new ForceFields::PositionConstraintContrib(
                          ff.get(), i, scoreShrug, 100.0)));
                }
            }
          }
          ff->initialize();
          const double ePre = params.logFragmentEnergies
                                  ? ff->calcEnergy()
                                  : std::numeric_limits<double>::quiet_NaN();
          const int nIt = ff->minimize(minIters, params.minimizeGradTol);
          double e = ff->calcEnergy();
          if (std::isnan(e)) e = std::numeric_limits<double>::max();
          if (params.logFragmentEnergies) {
            // nIt != 0 from RDKit's minimize() means it hit maxIters WITHOUT
            // converging. Under ShrugScore the restraints live in
            // ff->contribs(), so calcEnergy() (and hence the ranking score) =
            // MMFF + flat-bottom PENALTY.  Re-score with a restraint-free field
            // at the same coords to expose how much of `post` is penalty, not
            // real strain.
            double pureE = e;
            if (scoreShrug > 0.0) {
              std::unique_ptr<ForceFields::ForceField> ff2(
                  MMFF::constructForceField(frag, &props, 1.0e8, cid));
              if (ff2) {
                ff2->initialize();
                pureE = ff2->calcEnergy();
              }
            }
            BOOST_LOG(rdWarningLog)
                << "[fragE] class=" << fragmentClassName(fc) << " rot=" << rot
                << " nAt=" << nA << " pre=" << ePre << " post=" << e
                << " pure=" << pureE << " penalty=" << (e - pureE)
                << " delta=" << (ePre - e) << " iters=" << minIters
                << " converged=" << (nIt == 0) << "\n";
          }
          byEnergy.emplace_back(e, cid);
          if (scoreKeep)
            for (unsigned int i = 0; i < nA; ++i)
              c.setAtomPos(i, saved[i]);  // keep embed coords
        }
      } else {
        // single-point: type once, score every conformer at its own coords.
        std::unique_ptr<ForceFields::ForceField> ff(
            MMFF::constructForceField(frag, &props, 1.0e8, cids.front()));
        if (ff) {
          ff->initialize();
          std::vector<double> pos(3 * static_cast<size_t>(nA));
          for (const auto cid : cids) {
            const Conformer &c =
                frag.getConformer(static_cast<unsigned int>(cid));
            for (unsigned int i = 0; i < nA; ++i) {
              const RDGeom::Point3D &p = c.getAtomPos(i);
              pos[3 * i] = p.x;
              pos[3 * i + 1] = p.y;
              pos[3 * i + 2] = p.z;
            }
            double e = ff->calcEnergy(pos.data());
            if (std::isnan(e)) e = std::numeric_limits<double>::max();
            byEnergy.emplace_back(e, cid);
          }
        }
      }
    }
  }

  // Quick Filter by energy window
  if (eWindow > 0.0 && byEnergy.size() > 1) {
    double best = std::numeric_limits<double>::max();
    for (const auto &e : byEnergy) best = std::min(best, e.first);
    std::vector<std::pair<double, int>> within;
    within.reserve(byEnergy.size());
    for (const auto &e : byEnergy) {
      if (e.first <= best + eWindow) within.push_back(e);
    }
    if (!within.empty()) byEnergy.swap(within);
  }

  // Reduce to keepN: from lowest energy up, keep only RMSD-distinct geometries
  // (the pool + diverse-select recipe -- rigid -> ~1, flexible -> its real
  // number of puckers, up to maxConfs).
  {
    std::vector<int> order;
    if (!byEnergy.empty()) {
      std::sort(byEnergy.begin(), byEnergy.end());
      for (const auto &e : byEnergy) order.push_back(e.second);
    } else {
      order.assign(cids.begin(), cids.end());
    }
    // Prune the confs
    RMSDPruner acc(frag, poolRms);
    const unsigned int nAt = frag.getNumAtoms();
    std::vector<double> buf(3 * static_cast<size_t>(nAt));
    std::vector<int> kept;
    for (int cid : order) {
      if (kept.size() >= keepN) break;
      const Conformer &c = frag.getConformer(static_cast<unsigned int>(cid));
      for (unsigned int i = 0; i < nAt; ++i) {
        const RDGeom::Point3D &p = c.getAtomPos(i);
        buf[3 * i] = p.x;
        buf[3 * i + 1] = p.y;
        buf[3 * i + 2] = p.z;
      }
      if (acc.add(buf)) kept.push_back(cid);
    }

    // Tag each survivor with the MMFF energy that ranked it, so the assembler
    // can sum real energies rather than recomputing or approximating them.  See
    // kFragConfEnergy.
    std::map<int, double> eByCid;
    for (const auto &e : byEnergy) eByCid[e.second] = e.first;

    std::vector<Conformer *> keepConfs;
    keepConfs.reserve(kept.size());
    for (int cid : kept) {
      auto *nc =
          new Conformer(frag.getConformer(static_cast<unsigned int>(cid)));
      const auto it = eByCid.find(cid);
      if (it != eByCid.end()) nc->setProp(kFragConfEnergy, it->second);
      keepConfs.push_back(nc);
    }
    frag.clearConformers();
    unsigned int nid = 0;
    for (auto *c : keepConfs) {
      c->setId(nid++);
      frag.addConformer(c, false);
    }
  }

  // (The reduce-to-keepN block above already emitted the survivors
  // lowest-energy first, which is what the joiner consumes.)

  //  XXX FIX ME -> Is this actually ever used?
  // Tag each surviving conformer with a COARSE geometric state id (cluster by
  // heavy+exit RMSD): groups minor variations of the same ring pucker /
  // N-inversion state under one tag so the assembly's state-aware retention
  // isn't fooled by fine per-conformer differences (conformer indices
  // over-discriminate).
  {
    std::set<unsigned int> exitSet(dummies.begin(), dummies.end());
    MatchVectType stMap;
    for (const auto atom : frag.atoms()) {
      if (atom->getAtomicNum() > 1 || exitSet.count(atom->getIdx())) {
        stMap.emplace_back(static_cast<int>(atom->getIdx()),
                           static_cast<int>(atom->getIdx()));
      }
    }
    const double coarse =
        0.5;  // A; distinguishes puckers/N-faces, groups the rest
    std::vector<int> repCids;
    for (auto ci = frag.beginConformers(); ci != frag.endConformers(); ++ci) {
      const int cid = (*ci)->getId();
      int state = -1;
      for (size_t r = 0; r < repCids.size(); ++r) {
        double rms = std::numeric_limits<double>::max();
        try {
          RDGeom::Transform3D t;
          rms = MolAlign::getAlignmentTransform(frag, frag, t, cid, repCids[r],
                                                &stMap);
        } catch (...) {
        }
        if (rms < coarse) {
          state = static_cast<int>(r);
          break;
        }
      }
      if (state < 0) {
        state = static_cast<int>(repCids.size());
        repCids.push_back(cid);
      }
      (*ci)->setProp<std::string>("_fragCombo", std::to_string(state));
    }
  }

  // XXX FIX ME -> this might be vestigial
  // Un-cap: exit C/H -> atomic-number-0 dummy (its embedded position is the
  // exit vector), and remove the carbon-cap placeholder methyl H's.
  for (auto idx : dummies) {
    Atom *a = frag.getAtomWithIdx(idx);
    a->setAtomicNum(0);
    a->setNoImplicit(true);
    a->updatePropertyCache(false);
  }
  if (!capHs.empty()) {
    std::sort(capHs.rbegin(),
              capHs.rend());  // descending: removal keeps indices valid
    for (unsigned int hi : capHs) {
      frag.removeAtom(hi);
    }
  }
  frag.updatePropertyCache(false);
  if (detail::profiling()) {  // per-class embedding cost (speed breakdown vs
                              // )
    const int ci = static_cast<int>(fc);
    if (ci >= 0 && ci < detail::JoinProf::kNClass) {
      detail::prof().tEmbedClass[ci] += detail::nowNs() - tEmbed0;
      detail::prof().nEmbedClass[ci] += 1;
      detail::prof().nConfClass[ci] += frag.getNumConformers();
    }
  }
  return true;
}

}  // namespace

std::string Fraglib::cacheKey(const ROMol &frag) {
  // XXX FIX ME -> this is currently unused.  It might or might not be important
  //  at a later date to check enhanced stereo, but a single conf should
  //  not be a mixture so...
  // Canonical isomeric SMILES plus ONLY the enhanced-stereo CX tags: stable
  // (no coordinates / atom labels / other CX noise) but still distinguishes
  // fragments that differ solely in enhanced (AND/OR/relative) stereo.
  try {
    SmilesWriteParams ps;  // canonical + isomeric by default
    return MolToCXSmiles(frag, ps,
                         SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO);
  } catch (...) {
    return std::string();
  }
}

std::string Fraglib::generateKey(RWMol &frag, bool remap,
                                 std::vector<unsigned int> *outOrder) {
  static const char *kExitIsoProp = "_fraglibExitIso";

  // Exit-vector isotopes are per-cut bookkeeping for zipping, not fragment
  // identity: stash them so the key ignores which junction-index labels a
  // fragment carries.  They ride along through renumberAtoms and are restored
  // at the end.
  std::vector<unsigned int> exitDummies;
  for (auto atom : frag.atoms()) {
    if (atom->getAtomicNum() == 0 && atom->getIsotope() > 0) {
      atom->setProp<unsigned int>(kExitIsoProp, atom->getIsotope());
      exitDummies.push_back(atom->getIdx());
    }
  }

  // note: A dummy on a STEREOCENTRE that carries >=2 exit dummies is a "chiral
  // exit": we need any dummy or carbon to be labled differently for the RDKit
  // to preserve stereo erase the stereocentre.  Instead we give such dummies
  // DISTINCT small isotopes
  // ([1*],[2*],...), which keeps the stereocentre valid so the key stays
  // order-stable AND chirality-aware.  All other exits are zeroed (fully
  // label-agnostic).
  std::vector<unsigned int> chiralExits;
  {
    std::map<unsigned int, std::vector<unsigned int>> byCentre;
    for (unsigned int d : exitDummies) {
      const Atom *da = frag.getAtomWithIdx(d);
      for (const auto nbr : frag.atomNeighbors(da)) {
        if (nbr->getChiralTag() != Atom::CHI_UNSPECIFIED) {
          byCentre[nbr->getIdx()].push_back(d);
        }
      }
    }
    for (auto &kv : byCentre) {
      if (kv.second.size() >= 2) {
        for (unsigned int d : kv.second) {
          chiralExits.push_back(d);
        }
      }
    }
  }
  const std::set<unsigned int> chiralSet(chiralExits.begin(),
                                         chiralExits.end());
  for (unsigned int d : exitDummies) {
    if (!chiralSet.count(d)) {
      frag.getAtomWithIdx(d)->setIsotope(0);
    }
  }

  std::string key;
  try {
    if (chiralExits.empty() || chiralExits.size() > 6) {
      if (chiralExits.size() > 6) {  // pathological; give up on marking
        for (unsigned int d : chiralExits) {
          frag.getAtomWithIdx(d)->setIsotope(0);
        }
      }
      key = MolToSmiles(frag);
    } else {
      // Canonically DISAMBIGUATE the chiral exits: try every assignment of the
      // small isotopes {1..k} and keep the one with the
      // lexicographically-smallest canonical SMILES.  Deterministic and
      // structure-based, so the same fragment always resolves to the same marks
      // regardless of its original junction labels or atom order.  k is tiny
      // (usually 2).
      std::vector<unsigned int> perm(chiralExits.size());
      for (unsigned int i = 0; i < perm.size(); ++i) {
        perm[i] = i + 1;
      }
      std::vector<unsigned int> best = perm;
      do {
        for (unsigned int i = 0; i < chiralExits.size(); ++i) {
          frag.getAtomWithIdx(chiralExits[i])->setIsotope(perm[i]);
        }
        const std::string s = MolToSmiles(frag);
        if (key.empty() || s < key) {
          key = s;
          best = perm;
        }
      } while (std::next_permutation(perm.begin(), perm.end()));
      for (unsigned int i = 0; i < chiralExits.size(); ++i) {
        frag.getAtomWithIdx(chiralExits[i])->setIsotope(best[i]);
      }
    }

    if (remap || outOrder) {
      // Canonical atom order of the chosen assignment (recompute so
      // _smilesAtomOutputOrder matches `key`).  With explicit Hs present on
      // `frag` this order covers EVERY atom (heavy AND hydrogen), which lets
      // callers map coordinates atom-for-atom, Hs included.
      MolToSmiles(frag);
      std::vector<unsigned int> order;
      frag.getProp(common_properties::_smilesAtomOutputOrder, order);
      if (outOrder) {
        *outOrder = order;
      }
      if (remap && order.size() == frag.getNumAtoms()) {
        // renumberAtoms returns a fresh mol (really an RWMol); move it into
        // frag instead of copying.
        std::unique_ptr<RWMol> renum(
            static_cast<RWMol *>(MolOps::renumberAtoms(frag, order)));
        frag = std::move(*renum);
      }
    }
  } catch (...) {
    key.clear();
    if (outOrder) {
      outOrder->clear();
    }
  }

  for (auto atom : frag.atoms()) {  // restore the exit-vector isotopes
    unsigned int iso = 0;
    if (atom->getPropIfPresent(kExitIsoProp, iso)) {
      atom->setIsotope(iso);
      atom->clearProp(kExitIsoProp);
    }
  }
  return key;
}

Fraglib::~Fraglib() {
  for (auto &kv : d_fraglib) {
    delete kv.second;
  }
}

const RWMol *Fraglib::lookupOrEmbed(const std::string &key,
                                    const std::function<RWMol *()> &make,
                                    bool cache,
                                    std::unique_ptr<RWMol> &owned) const {
  {
    std::lock_guard<std::mutex> lock(d_mutex);
    auto it = d_fraglib.find(key);
    if (it != d_fraglib.end()) {
      return it->second;
    }
  }
  RWMol *embedded = nullptr;
  try {
    embedded = make();
  } catch (...) {
    // we can't embed this
    embedded = nullptr;
  }
  if (!embedded) {
    // If we can't embed record the failure as a nullptr
    //  so we don't try again.  Some of these failures can
    //  take eons before they fail
    if (cache) {
      std::lock_guard<std::mutex> lock(d_mutex);
      d_fraglib.emplace(key, nullptr);
    }
    return nullptr;
  }
  std::lock_guard<std::mutex> lock(d_mutex);
  auto it = d_fraglib.find(key);
  if (it != d_fraglib.end()) {  // another thread won the race
    delete embedded;
    return it->second;
  }
  if (cache) {
    d_fraglib[key] = embedded;  // the map owns it from here
    return embedded;
  }
  owned.reset(embedded);  // caller owns it; freed when `owned` leaves scope
  return embedded;
}

ROMOL_SPTR Fraglib::get(const ROMol &frag, bool cache) const {
  // We need to generate a cache key here for lookup.  To do this we always add
  //  hydrogens to generate the cache key.  Then, we look it up
  //  if it exists, we're golden, otherwise we look at the fraglib
  //  params and generate (and store) a new one.
  auto q = boost::make_shared<RWMol>(frag);
  MolOps::addHs(
      *q);  // key/order (and thus the stored master) are FULL-ATOM canonical
            // -- Hs included -- so getConformersOnto maps coords atom-for-atom
  const std::string key = generateKey(*q, /*remap=*/true);
  if (key.empty()) {
    return ROMOL_SPTR();  // could not canonicalize
  }

  std::unique_ptr<RWMol> owned;
  const RWMol *cachedFrag = lookupOrEmbed(
      key,
      [&]() -> RWMol * {
        auto *embedded = new RWMol(*q);
        if (!embedFragmentInPlace(*embedded, d_params)) {
          delete embedded;
          return nullptr;
        }
        return embedded;
      },
      cache, owned);
  if (!cachedFrag) {
    return ROMOL_SPTR();
  }

  // Return a new copy with the original isotope labels if necessary.
  //  note that the original ones were changed in order to preserve
  //  dummy atom stereo during the cache lookup
  auto out = boost::make_shared<RWMol>(*cachedFrag);
  for (unsigned int i = 0; i < q->getNumAtoms(); ++i) {
    if (q->getAtomWithIdx(i)->getAtomicNum() == 0) {
      out->getAtomWithIdx(i)->setIsotope(q->getAtomWithIdx(i)->getIsotope());
    }
  }
  return out;
}

bool Fraglib::getConformerCoords(RWMol &frag, unsigned int nMolAtoms,
                                 const std::string &molIdxProp,
                                 std::vector<std::vector<RDGeom::Point3D>> &out,
                                 std::vector<double> *energiesOut,
                                 bool cache) const {
  std::vector<unsigned int> order;  // order[canonicalPosition] = frag atom index
  const std::string key = generateKey(frag, /*remap=*/false, &order);
  if (key.empty() || order.size() != frag.getNumAtoms()) {
    return false;
  }

  std::unique_ptr<RWMol> owned;
  const RWMol *cachedFrag = lookupOrEmbed(
      key,
      [&]() -> RWMol * {
        // Embed in the full-atom canonical order (renumber THIS fragment to it).
        auto *embedded = static_cast<RWMol *>(MolOps::renumberAtoms(frag, order));
        if (!embedFragmentInPlace(*embedded, d_params)) {
          delete embedded;
          return nullptr;
        }
        return embedded;
      },
      cache, owned);
  if (!cachedFrag) {
    return false;
  }
  const unsigned int nConf = cachedFrag->getNumConformers();
  if (nConf == 0 || cachedFrag->getNumAtoms() != frag.getNumAtoms()) {
    if (d_params.FRAGLIB_TRACE) {
      static std::atomic<int> nAtomMiss{0};
      std::string mkey;
      try {
        mkey = MolToSmiles(*cachedFrag);
      } catch (...) {
      }
      std::string fkey;
      try {
        fkey = MolToSmiles(frag);
      } catch (...) {
      }
      // A cached fragment whose atom count disagrees with the query means the
      // cache key collided -- serving it would corrupt coordinates.
      BOOST_LOG(rdErrorLog)
          << "[fraglibTRACE] atom-count mismatch #" << ++nAtomMiss
          << " cachedFrag=" << cachedFrag->getNumAtoms() << "(" << mkey
          << ") frag=" << frag.getNumAtoms() << "(" << fkey << ") key=[" << key
          << "]\n";
    }
    return false;
  }
  if (d_params.FRAGLIB_TRACE) {
    static std::atomic<int> nServed{0};
    if ((++nServed % 50) == 1)
      BOOST_LOG(rdWarningLog) << "[fraglibTRACE] served coords from cachedFrag #"
                              << nServed.load() << "\n";
  }

  std::vector<unsigned int> canonPos(frag.getNumAtoms());
  for (unsigned int p = 0; p < order.size(); ++p) {
    canonPos[order[p]] = p;
  }
  std::vector<int> molIdx(frag.getNumAtoms(), -1);
  for (const auto a : frag.atoms()) {
    int mi = -1;
    if (a->getPropIfPresent(molIdxProp, mi)) molIdx[a->getIdx()] = mi;
  }

  // Add the confs to the given frag
  out.assign(nConf, std::vector<RDGeom::Point3D>(nMolAtoms));
  if (energiesOut)
    energiesOut->assign(nConf, std::numeric_limits<double>::quiet_NaN());
  for (unsigned int ci = 0; ci < nConf; ++ci) {
    const Conformer &mc = cachedFrag->getConformer(ci);
    if (energiesOut) {
      double e = std::numeric_limits<double>::quiet_NaN();
      if (mc.getPropIfPresent(kFragConfEnergy, e)) (*energiesOut)[ci] = e;
    }
    auto &buf = out[ci];
    for (unsigned int a = 0; a < frag.getNumAtoms(); ++a) {
      const int mi = molIdx[a];
      if (mi >= 0 && mi < static_cast<int>(nMolAtoms)) {
        buf[mi] = mc.getAtomPos(canonPos[a]);
      }
    }
  }
  return true;
}

size_t Fraglib::size() const {
  std::lock_guard<std::mutex> lock(d_mutex);
  size_t n = 0;
  for (const auto &kv : d_fraglib) {
    if (kv.second) {  // skip tombstones for fragments that cannot be embedded
      ++n;
    }
  }
  return n;
}

size_t Fraglib::numUnembeddable() const {
  std::lock_guard<std::mutex> lock(d_mutex);
  size_t n = 0;
  for (const auto &kv : d_fraglib) {
    if (!kv.second) {
      ++n;
    }
  }
  return n;
}

bool Fraglib::markUnembeddable(const ROMol &frag) {
  RWMol withHs(frag);
  MolOps::addHs(withHs);
  const std::string key = generateKey(withHs, /*remap=*/true);
  if (key.empty()) {
    return false;
  }
  std::lock_guard<std::mutex> lock(d_mutex);
  return d_fraglib.emplace(key, nullptr).second;
}

std::optional<unsigned int> Fraglib::numFragmentConfs(
    const ROMol &frag) const {
  RWMol withHs(frag);
  MolOps::addHs(withHs);
  const std::string key = generateKey(withHs, /*remap=*/true);
  if (key.empty()) {
    return std::nullopt;
  }
  std::lock_guard<std::mutex> lock(d_mutex);
  const auto found = d_fraglib.find(key);
  if (found == d_fraglib.end()) {
    return std::nullopt;  // never attempted: the caller must embed it
  }
  if (!found->second) {
    return 0u;  // tombstone: attempted and known unembeddable, do not retry
  }
  return found->second->getNumConformers();
}

// Simple IO class for the fraglib
namespace {
// Fraglib format version.
//  XXX FIX ME Should we use boost serialize?

constexpr char kFraglibMagic[8] = {'F', 'R', 'A', 'G', 'L', 'I', 'B', '1'};
template <typename T>
void wRaw(std::ostream &os, const T &v) {
  os.write(reinterpret_cast<const char *>(&v), sizeof(T));
}
template <typename T>
void rRaw(std::istream &is, T &v) {
  is.read(reinterpret_cast<char *>(&v), sizeof(T));
}
void wBlob(std::ostream &os, const std::string &s) {
  wRaw(os, static_cast<std::uint64_t>(s.size()));
  os.write(s.data(), static_cast<std::streamsize>(s.size()));
}
std::string rBlob(std::istream &is) {
  std::uint64_t n = 0;
  rRaw(is, n);
  std::string s(static_cast<size_t>(n), '\0');
  if (n) is.read(&s[0], static_cast<std::streamsize>(n));
  return s;
}
}  // namespace

void Fraglib::serialize(std::ostream &os) const {
  std::lock_guard<std::mutex> lock(d_mutex);
  os.write(kFraglibMagic, sizeof(kFraglibMagic));

  // Serialize the params
  wRaw(os, d_params.numConfsPerFragment);
  wRaw(os, static_cast<int>(d_params.fragmentEmbedMode));
  wRaw(os, d_params.randomSeed);
  wRaw(os, static_cast<int>(d_params.minimizeMode));
  wRaw(os, d_params.shrugDisplacement);
  wRaw(os, d_params.minimizeGradTol);
  wRaw(os, d_params.minimizeMaxIters);
  wRaw(os, d_params.perClassEmbedding);
  wRaw(os, d_params.energyWindow);

  // n.b. we need to pickle failures as well (empty mols/pickles)
  //  as these are sentinels for failed embeddings
  wRaw(os, static_cast<std::uint64_t>(d_fraglib.size()));
  for (const auto &kv : d_fraglib) {
    wBlob(os, kv.first);
    std::string pkl;
    if (kv.second) {
      // Save annotations on molecules and atoms
      MolPickler::pickleMol(*kv.second, pkl,
                            PicklerOps::MolProps | PicklerOps::AtomProps);
    }
    wBlob(os, pkl);  // empty blob == marked as unembeddable
  }
}

void Fraglib::initFromStream(std::istream &is) {
  std::lock_guard<std::mutex> lock(d_mutex);
  for (auto &kv : d_fraglib) {
    delete kv.second;
  }
  d_fraglib.clear();

  char magic[8] = {0};
  is.read(magic, sizeof(magic));
  if (std::memcmp(magic, kFraglibMagic, sizeof(magic)) != 0) {
    // Report the magic we WANT, from the constant -- a hardcoded version in
    // this message goes stale the first time the format is bumped.
    throw std::runtime_error(
        "Fraglib::initFromStream: bad magic (not a " +
        std::string(kFraglibMagic, sizeof(kFraglibMagic)) + " file), got \"" +
        std::string(magic, sizeof(magic)) + "\"");
  }

  int embedMode = 0, minMode = 0;
  rRaw(is, d_params.numConfsPerFragment);
  rRaw(is, embedMode);
  d_params.fragmentEmbedMode = static_cast<FragmentEmbedMode>(embedMode);
  rRaw(is, d_params.randomSeed);
  rRaw(is, minMode);
  d_params.minimizeMode = static_cast<FragmentMinimize>(minMode);
  rRaw(is, d_params.shrugDisplacement);
  rRaw(is, d_params.minimizeGradTol);
  rRaw(is, d_params.minimizeMaxIters);
  rRaw(is, d_params.perClassEmbedding);
  rRaw(is, d_params.energyWindow);

  std::uint64_t n = 0;
  rRaw(is, n);
  for (std::uint64_t i = 0; i < n; ++i) {
    std::string key = rBlob(is);
    std::string pkl = rBlob(is);
    if (pkl.empty()) {
      d_fraglib[key] = nullptr;  // tombstone: known-unembeddable
      continue;
    }
    auto *m = new RWMol();
    MolPickler::molFromPickle(pkl, m);
    d_fraglib[key] = m;
  }
}

void Fraglib::writeSDF(std::ostream &os) const {
  std::lock_guard<std::mutex> lock(d_mutex);
  for (const auto &kv : d_fraglib) {
    RWMol *m = kv.second;  
    if (!m || m->getNumConformers() == 0) {
      continue;
    }

    const unsigned int nA = m->getNumAtoms();
    std::vector<double> energies(m->getNumConformers(),
                                 std::numeric_limits<double>::quiet_NaN());
    try {
      // XXX FIX ME -> I think this is dead code
      // The stored cachedFrag keeps exit vectors as atomic-number-0 dummies, which
      // MMFF cannot type.  Mirror embedding: cap the dummies (carbon or H per
      // params) on a throwaway copy so the single-point energy matches the
      // geometry we generated.
      RWMol capped(*m);
      for (auto *at : capped.atoms()) {
        if (at->getAtomicNum() == 0) {
          at->setAtomicNum(6);  // dummies always carbon-capped
          at->setNoImplicit(false);
          at->setIsotope(0);
        }
      }
      MolOps::sanitizeMol(capped);
      MMFF::MMFFMolProperties props(capped, d_params.ffVariant);
      if (props.isValid()) {
        props.setMMFFEleTerm(false);
        std::unique_ptr<ForceFields::ForceField> ff(MMFF::constructForceField(
            capped, &props, 1.0e8, (*capped.beginConformers())->getId()));
        if (ff) {
          ff->initialize();
          std::vector<double> pos(3 * static_cast<size_t>(nA));
          size_t k = 0;
          for (auto ci = capped.beginConformers(); ci != capped.endConformers();
               ++ci, ++k) {
            const Conformer &c = **ci;
            for (unsigned int i = 0; i < nA; ++i) {
              const RDGeom::Point3D &p = c.getAtomPos(i);
              pos[3 * i] = p.x;
              pos[3 * i + 1] = p.y;
              pos[3 * i + 2] = p.z;
            }
            energies[k] = ff->calcEnergy(pos.data());
          }
        }
      }
    } catch (...) {
      // untypeable even when capped -> leave energies NaN, still dump the
      // geometry
    }
    
    RWMol tagged(*m);
    // Mark exits as RGroups for round tripping
    unsigned int rlabel = 0;
    for (auto *at : tagged.atoms()) {
      if (at->getAtomicNum() == 0) {
        at->setProp<unsigned int>(common_properties::_MolFileRLabel, ++rlabel);
      }
    }
    tagged.setProp("fragment_smiles", kv.first);
    tagged.setProp("fragment_class", std::string(fragmentClassName(classifyFragment(tagged))));
    size_t k = 0;
    for (auto ci = m->beginConformers(); ci != m->endConformers(); ++ci, ++k) {
      const int cid = (*ci)->getId();
      tagged.setProp("conf_id", cid);
      tagged.setProp("mmff_energy", energies[k]);
      os << SDWriter::getText(tagged, cid);
    }
  }
}

std::string Fraglib::add(RWMol &frag) const {
  // XXX FIX ME -> this may eventually be dead code
  // Take an external frag and place it in the fraglib
  std::vector<unsigned int> order;
  const std::string key = generateKey(frag, /*remap=*/true, &order);
  if (key.empty()) {
    return std::string();
  }
  std::lock_guard<std::mutex> lock(d_mutex);
  if (d_fraglib.find(key) == d_fraglib.end()) {
    d_fraglib[key] = new RWMol(frag);
  }
  return key;
}

}  // namespace RDKit
