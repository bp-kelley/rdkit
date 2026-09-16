//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "FragmentConfGen.h"
#include "Embedder/Embedder.h"
#include "Joiner/FragmentJoiner.h"
#include "Joiner/JoinerProfiling.h"
#include "Search/InterFragScore.h"
#include "Search/RigidRotorSearch.h"
#include "Utils/ParamsIO.h"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Geometry/Transform3D.h>
#include <ForceField/ForceField.h>
#include <RDGeneral/RDLog.h>

#include <algorithm>
#include <atomic>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <mutex>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace RDKit {

namespace {

//! get the largest heavy component of a molecule
boost::shared_ptr<RWMol> largestHeavyComponent(const ROMol &mol) {
  auto heavy = boost::make_shared<RWMol>(mol);
  MolOps::removeHs(*heavy);
  RWMOL_SPTR_VECT comps = getRWMolFrags(*heavy, /*sanitizeFrags=*/true);
  if (comps.size() > 1) {
    RWMOL_SPTR largest;
    unsigned int bestHeavy = 0;
    for (const auto &c : comps) {
      unsigned int h = 0;
      for (const auto atom : c->atoms()) {
        if (atom->getAtomicNum() > 1) {
          ++h;
        }
      }
      if (h > bestHeavy) {
        bestHeavy = h;
        largest = c;
      }
    }
    if (largest) {
      heavy = largest;  // fresh fragment; take it, don't copy
    }
  }
  return heavy;
}

//! Label exit rotors
constexpr unsigned int EXIT_LABEL_BASE = 1000;

//! Construct fragment library params from the conformer generator params
EmbedderParams getEmbedderParams(const FragmentConfGenParams &p, int seed) {
  EmbedderParams flp = p.embedding;
  flp.randomSeed = seed;
  flp.EMBEDDER_TRACE = p.diagnostics.EMBEDDER_TRACE;
  flp.ffVariant = p.joiner.ffVariant;
  return flp;
}

//! Symmetric spinner
//!  true if we are spinning symmetric atoms.  I.e. CF3, etc.
bool isSymmetricSpinner(const ROMol &mol, const Atom *y, unsigned int fromIdx) {
  std::vector<const Atom *> others;
  for (const auto nbr : mol.atomNeighbors(y)) {
    if (nbr->getIdx() != fromIdx) {
      others.push_back(nbr);
    }
  }
  if (others.size() != 3) {
    return false;
  }
  const int z = others.front()->getAtomicNum();
  for (const auto *o : others) {
    if (o->getAtomicNum() != z || o->getDegree() != 1) {
      return false;
    }
  }
  return true;
}

//! AUTO output-conformer count: scale the ensemble with the molecule's
//! flexibility.
/*!
  A flat cap truncates floppy molecules hardest -- their accessible space grows
  combinatorially with rotor count while the cap does not.  Measured on Platinum
  rot>=11, the searches generate a sub-1A pose for ~2/3 of molecules but a
  250-conformer cap keeps one for under half: the poses are produced and then
  discarded by selection.

  The bands come from the saturation sweep in docs/output-selection.md -- the
  knee of %<1 vs ensemble size per rotor band, not the highest value reached
  (the metric is best-of-ensemble, so it can only flatten, never turn down).
*/
int autoOutputConfs(unsigned int nRotors) {
  if (nRotors <= 4) return 50;
  if (nRotors <= 6) return 100;
  if (nRotors <= 8) return 250;
  if (nRotors <= 10) return 500;
  return 1000;
}

}  // namespace

std::string FragmentConfGenParams::validate() const {
  // The joiner and the search own the rules for their own parameters -- ask
  // them rather than restating the checks here and letting the two drift.
  if (const std::string err = joiner.validate(); !err.empty()) {
    return err;
  }
  if (const std::string err = validateSearchParams(search, joiner.ffVariant);
      !err.empty()) {
    return err;
  }

  // What is left is this layer's own.
  if (!labelFFVariant.empty() && !isValidFF(labelFFVariant)) {
    return std::string("labelFFVariant must be empty, ") +
           " (got \"" + labelFFVariant + "\")";
  }
  if (!isValidFF(embedding.ffVariant)) {
    std::string s;
    const size_t sz = sizeof(kFFVariants)/sizeof(kFFVariants[0]);
    for(size_t i=0; i<sz; i++) {
      s += kFFVariants[0];
      if(i<sz-1)
	s+=",";	 
    }

    return std::string("embedding.ffVariant must be one of") + s +
           " (got \"" + embedding.ffVariant + "\")";
  }
  if (energyWindow < 0.0) {
    return "energyWindow must be >= 0";
  }
  if (embedding.minimizeGradTol < 0.0) {
    return "embedding.minimizeGradTol must be >= 0";
  }
  if (embedding.shrugDisplacement < 0.0) {
    return "embedding.shrugDisplacement must be >= 0";
  }

  // Valide the embedder against the current params
  if (embedder && !sameEmbeddingType(embedder->params(),
                                    getEmbedderParams(*this, randomSeed))) {
    return "embedder was built with a different embedding than the current "
           "FragmentConfGenParams";
  }
  return {};  // coherent
}

void FragmentConfGenParams::setEmbedder(std::shared_ptr<Embedder> lib) {
  auto prev = std::move(embedder);
  embedder = std::move(lib);
  const std::string err = validate();
  if (!err.empty()) {
    embedder = std::move(prev);  // leave the params as we found them
    throw std::invalid_argument("FragmentConfGen Parameter error: " + err);
  }
}

FragmentConfGen createFragmentConfGen(FragmentConfGenParams params) {
  const std::string err = params.validate();
  if (!err.empty()) {
    throw std::invalid_argument("FragmentConfGen Parameter error: " + err);
  }
  return FragmentConfGen(std::move(params));
}

FragmentConfGen::FragmentConfGen(FragmentConfGenParams params)
    : d_params(std::move(params)) {
  if (!d_params.embedder) {
    EmbedderParams flp = getEmbedderParams(d_params, d_params.randomSeed);
    d_params.embedder = std::make_shared<Embedder>(flp);
  }
  // Set once, here, where the parameters are fixed.  The joiner profiler is a
  // PROCESS-GLOBAL flag, so doing this per-build let concurrent generators
  // flip each other's setting mid-run.
  detail::setJoinerProfiling(d_params.diagnostics.ASM_PROFILE);
}

RotatableBonds FragmentConfGen::findRotatableBonds(const ROMol &mol,
                                                   bool sampleTrivial,
                                                   bool wholeAcyclicFragments) {
  RotatableBonds res;
  // Rotatable-bond definition from calcNumRotatableBonds (NonStrict): acyclic
  // single bond between two non-terminal, non-triple-bonded atoms.
  static const std::string rotSmarts = "[!$(*#*)&!D1]-!@[!$(*#*)&!D1]";
  std::unique_ptr<ROMol> patt(SmartsToMol(rotSmarts));
  if (!patt) {
    return res;
  }

  // XXX FIX ME -> this feels fragile and needs to be reconsidered
  // NEVER cut the conjugated C(=X)-Y bond of an
  // amide/urea/carbamate/ester/thioamide (the carbonyl carbon directly bonded
  // to N/O/S).  That bond is planar (partial double bond); cutting it separates
  // the C=O from its heteroatom, and the isolated 1-2 heavy-atom pieces cannot
  // reproduce the planar sp2 geometry This borks the exit geometry.
  //  n.b. We only protect the C(=X)-Y bond
  std::set<unsigned int> amideBonds;
  {
    static const std::string amideSmarts = "[CD3]([#7,#8,#16])=[N,O,S]";
    std::unique_ptr<ROMol> amide(SmartsToMol(amideSmarts));
    if (amide) {
      std::vector<MatchVectType> am;
      SubstructMatch(mol, *amide, am, /*uniquify=*/true);
      for (const auto &m :
           am) {  // m[0]=carbonyl C, m[1]=the single-bonded heteroatom
        const Bond *b = mol.getBondBetweenAtoms(m[0].second, m[1].second);
        if (b) amideBonds.insert(b->getIdx());
      }
    }
  }

  // We have two potential driving sets
  //  inter - bonds between fragments
  //  intra - bonds IN fragments if using wholeAcyclicFragments
  // Determine both
  std::vector<MatchVectType> matches;
  SubstructMatch(mol, *patt, matches, /*uniquify=*/true);
  std::set<unsigned int> bonds;
  std::map<unsigned int, IntraRotorType> intra;
  for (const auto &m : matches) {
    const Bond *b = mol.getBondBetweenAtoms(m[0].second, m[1].second);
    if (!b) {
      continue;
    }
    if (!sampleTrivial &&
        (isSymmetricSpinner(mol, b->getBeginAtom(), b->getEndAtomIdx()) ||
         isSymmetricSpinner(mol, b->getEndAtom(), b->getBeginAtomIdx()))) {
      continue;  // trivial rotor: neither cut nor driven
    }
    if (b->getStereo() == Bond::STEREOATROPCW ||
        b->getStereo() == Bond::STEREOATROPCCW) {
      // We need to treat atropisomers differently as we can't change
      //  the stereochemistry during driving.
      intra[b->getIdx()] = IntraRotorType::Atropisomer;
      continue;
    }
    if (amideBonds.count(b->getIdx())) {
      // AMIDE Bond -> sample at 0/180  C(=X)-Y:
      intra[b->getIdx()] = IntraRotorType::PlanarAmide;
      continue;
    }
    bonds.insert(b->getIdx());
  }

  // if wholeAcyclicFragments: we oversample the entire fragment
  //  and only drive the junction bond.
  //  NOTE -> this was the initial driver for synthon searching
  if (wholeAcyclicFragments && !bonds.empty()) {
    const RingInfo *ri = mol.getRingInfo();
    if (ri && ri->isInitialized()) {
      std::set<unsigned int> ringAdjacent;
      for (unsigned int bi : bonds) {
        const Bond *b = mol.getBondWithIdx(bi);
        if (ri->numAtomRings(b->getBeginAtomIdx()) > 0 ||
            ri->numAtomRings(b->getEndAtomIdx()) > 0) {
          ringAdjacent.insert(bi);
        } else {
          intra[bi] = IntraRotorType::Free;  // uncut chain-chain bond
        }
      }
      bonds.swap(ringAdjacent);
    }
  }

  res.inter.assign(bonds.begin(), bonds.end());
  for (const auto &kv : intra) {
    res.intra.push_back({kv.first, kv.second});
  }
  return res;
}

std::vector<ROMOL_SPTR> FragmentConfGen::fragmentAndEmbed(
    const ROMol &mol, std::vector<unsigned int> *linkBondsOut,
    const std::vector<unsigned int> *linkBonds) const {
  RDLog::LogStateSetter blocker;

  auto heavy = largestHeavyComponent(mol);

  // XXX Fix me -> we find link bonds twice, let's consolidate this
  if (!heavy->getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(*heavy);
  }
  auto links = linkBonds ? *linkBonds
                        : findLinkBonds(*heavy, d_params.sampleTrivialRotors,
                                        d_params.wholeAcyclicFragments);

  if (linkBondsOut) {
    *linkBondsOut = links;
  }

  std::vector<ROMOL_SPTR> fragments;
  if (links.empty()) {
    auto whole = d_params.embedder->get(*heavy);
    if (whole) {
      fragments.push_back(whole);
    }
    return fragments;
  }

  std::vector<std::pair<unsigned int, unsigned int>> dummyLabels;
  dummyLabels.reserve(links.size());
  for (size_t i = 0; i < links.size(); ++i) {
    const unsigned int label = EXIT_LABEL_BASE + static_cast<unsigned int>(i);
    dummyLabels.emplace_back(label, label);
  }
  std::unique_ptr<ROMol> fragged(
      MolFragmenter::fragmentOnBonds(*heavy, links, true, &dummyLabels));
  if (!fragged) {
    return fragments;
  }

  const bool sanitizeFrags = true;
  std::vector<ROMOL_SPTR> pieces = MolOps::getMolFrags(*fragged, sanitizeFrags);
  for (auto &piece : pieces) {
    auto frag = d_params.embedder->get(*piece);
    if (!frag) {
      return {};  // a fragment failed to embed -> no assembly is possible
    }
    fragments.push_back(frag);
  }
  return fragments;
}

FragmentConfGenResult FragmentConfGen::build(
    const ROMol &mol, const std::vector<unsigned int> *linkBonds) const {
  FragmentConfGenResult result;

  if (d_params.maxRotatableBonds > 0) {
    const RotatableBonds rb = findRotatableBonds(
        mol, d_params.sampleTrivialRotors, d_params.wholeAcyclicFragments);
    const size_t nRot = rb.inter.size() + rb.intra.size();
    if (nRot > static_cast<size_t>(d_params.maxRotatableBonds)) {
      result.status = FragConfGenResultType::TOO_MANY_ROTORS;
      result.numLinkBonds = static_cast<unsigned int>(rb.inter.size());
      return result;
    }
  }

  // Generate the (hopefully energy minimized) conformer ensemble for this
  // molecule.
  buildEnsemble(mol, result, linkBonds);

  // Report the best energy.
  //  XXX FIX ME -> the ensemble should already be sorted from lowest to highest
  //  energy so this full pass might be redundant
  if (!result.conformers.empty()) {
    double best = std::numeric_limits<double>::max();
    for (const auto &c : result.conformers) {
      double e = 0.0;
      if (c->getPropIfPresent<double>(kEnergyProp, e)) {
        best = std::min(best, e);
      }
    }
    if (std::isfinite(best)) {
      result.bestEnergy = best;
    }
  }
  return result;
}

// The workhorse: joins fragments into reasonable linkage
// geometries, then rotorsearches the rotatable bonds
void FragmentConfGen::buildEnsemble(
    const ROMol &mol, FragmentConfGenResult &result,
    const std::vector<unsigned int> *linkBonds) const {
  RDLog::LogStateSetter blocker;
  const int seed = d_params.randomSeed >= 0 ? d_params.randomSeed : 0xf00d;

  // Embed frags
  const EmbedderParams want = getEmbedderParams(d_params, seed);
  std::shared_ptr<Embedder> privateLib;
  // Compatibility with a supplied library is a parameter-level invariant,
  // established once by FragmentConfGenParams::validate()/setEmbedder().
  const Embedder *lib = d_params.embedder.get();
  if (!lib) {
    privateLib = std::make_shared<Embedder>(want);
    lib = privateLib.get();
  }

  FragmentJoinerParams gp = d_params.joiner;
  gp.diagnostics = d_params.diagnostics;
  gp.sampleTrivialRotors = d_params.sampleTrivialRotors;
  gp.wholeAcyclicFragments = d_params.wholeAcyclicFragments;

  RigidRotorSearchParams sp = d_params.search;
  sp.diagnostics = d_params.diagnostics;
  sp.randomSeed = static_cast<unsigned int>(seed);
  sp.energyWindow =
      d_params.energyWindow > 0 ? d_params.energyWindow * 2.5 : 25.0;
  sp.timeBudgetMs = d_params.timeBudgetMs;

  FragmentJoinerInput in = buildFragmentJoinerInput(
      mol, d_params.embedding.numConfsPerFragment, seed,
      d_params.joiner.ffVariant,
      d_params.embedding.fragmentEmbedMode == FragmentEmbedMode::DG,
      d_params.embedding.minimizeMode != FragmentMinimize::None, lib, &gp,
      linkBonds);

  if (!in.mol) {
    return;
  }
  result.numLinkBonds = static_cast<unsigned int>(in.junctions.size());
  result.numFragments = static_cast<unsigned int>(in.fragments.size());
  // total flexibility = cut junctions + rotatable bonds left inside fragments;
  // the AUTO output cap scales with this, not with the junction count alone
  const size_t nIntraRotors = in.intraRotorBonds.size();

  const std::string labelVariant = d_params.labelFFVariant.empty()
                                       ? d_params.joiner.ffVariant
                                       : d_params.labelFFVariant;
  auto fullE = makeFullFFScoreFn(*in.mol, /*electrostatics=*/false, labelVariant);
  if (!fullE) {
    return;
  }

  //! Materialize the buffer into an RDKit Conformer
  auto materializeOne = [&](const std::vector<double> &coords) {
    auto m = boost::make_shared<RWMol>(*in.mol);
    auto *conf = new Conformer(m->getNumAtoms());
    for (unsigned int a = 0; a < m->getNumAtoms(); ++a) {
      conf->setAtomPos(a, RDGeom::Point3D(coords[3 * a], coords[3 * a + 1],
                                          coords[3 * a + 2]));
    }
    m->clearConformers();
    m->addConformer(conf, true);
    if (fullE) {
      m->setProp(kEnergyProp, fullE(coords.data(), m->getNumAtoms()));
      m->setProp(kForceFieldProp, labelVariant);
    }
    result.conformers.push_back(m);
    return m;
  };

  // scratch space for coords
  auto materialize = [&](const std::vector<std::vector<double>> &buffers) {
    for (const auto &coords : buffers) {
      materializeOne(coords);
    }
  };

  if (in.fragments.size() <= 1) {
    // no rotatable bonds: just return lowest energy embedded frag,
    std::vector<std::vector<double>> buffers;
    const unsigned int n = in.mol->getNumAtoms();
    if (!in.fragments.empty()) {
      for (const auto &c : in.fragments[0].confs) {
        std::vector<double> flat(3 * static_cast<size_t>(n));
        for (unsigned int a = 0; a < n; ++a) {
          flat[3 * a] = c.pos[a].x;
          flat[3 * a + 1] = c.pos[a].y;
          flat[3 * a + 2] = c.pos[a].z;
        }
        buffers.push_back(std::move(flat));
      }
    }
    materialize(buffers);
  } else {
    if (d_params.diagnostics.FRAGCG_DUMP_PARAMS) {
      // Dump params used for conf gen
      size_t totalFragConfs = 0;
      for (const auto &f : in.fragments) totalFragConfs += f.confs.size();
      const size_t nRot = in.junctions.size();
      RigidRotorSearchMode eff = sp.searchMode;
      if (eff == RigidRotorSearchMode::Auto) {
        eff = RigidRotorSearchParams::autoModeForRotors(nRot);
        if (eff == RigidRotorSearchMode::Thompson &&
            !(sp.thompsonBudget > 0 || sp.thompson.autoBudget)) {
          eff = RigidRotorSearchMode::Tree;
        }
      }

      {
        const std::string cfg = fragmentConfGenParamsToString(
            resolvedFragmentConfGenParams(d_params));
        static std::mutex dumpMtx;
        static std::set<size_t> dumped;
        const size_t key = std::hash<std::string>{}(cfg);
        bool isNew = false;
        {
          std::lock_guard<std::mutex> lk(dumpMtx);
          isNew = dumped.insert(key).second;
        }
        if (isNew) {
          BOOST_LOG(rdWarningLog)
              << "[fragcg] ===== RESOLVED PARAMETERS (effective values, AUTO "
                 "already derived) =====\n"
              << cfg << "[fragcg] ===== end =====\n";
        }
      }
      BOOST_LOG(rdWarningLog)
          << "[fragcg] per-molecule: rotors=" << nRot
          << " fragments=" << in.fragments.size()
          << " fragConfs=" << totalFragConfs
          << " searchMode=" << searchModeName(eff)
          << (sp.searchMode == RigidRotorSearchMode::Auto ? " (via Auto)" : "")
          << " sampler=" << (sp.torsionSampler ? "set" : "null") << "\n";
    }
    //  1. join the fragments, generating new ones when necessary
    //  2. run the rotor search finding the torsions that minimize energy
    //  3. optionally minimize ensemble (this shouldn't be necessary if 1+2 did
    //  a good job)
    //  4. sort by energy, prune by desired # confs and return
    const auto ctx = joinFragments(in, gp);
    if (!ctx.isValid()) {
      // An empty ensemble already surfaces as FF_FAIL; the joiner's status is
      // the reason WHY, which FragConfGenResultType cannot yet carry.
      BOOST_LOG(rdWarningLog) << "[fragcg] join failed: "
                              << fragmentJoinerStatusMessage(ctx.status)
                              << "\n";
    } else {
      const auto _searchStart = std::chrono::steady_clock::now();
      auto search = runRigidRotorSearch(ctx, sp);
      if (search.timedOut) result.status = FragConfGenResultType::TIMED_OUT;
      result.joinerAssemblyNs =
          std::chrono::duration_cast<std::chrono::nanoseconds>(
              std::chrono::steady_clock::now() - _searchStart)
              .count();
      result.joinerBudget = search.budget;
      result.conformers.reserve(result.conformers.size() +
                                search.results.size());
      for (const auto &r : search.results) {
        materializeOne(r.coords)->setProp(kInterFragScoreProp, r.score);
      }
    }
  }

  // rerank output confs if requested
  //  This could be expensive
  if (d_params.rankByBasinEnergy) {
    for (auto &m : result.conformers) {
      if (!m || m->getNumConformers() == 0) continue;
      try {
        RWMol tmp(*m);
        MMFF::MMFFMolProperties props(tmp, labelVariant);
        if (!props.isValid()) continue;
        props.setMMFFEleTerm(false);
        const int cid = tmp.getConformer().getId();
        std::unique_ptr<ForceFields::ForceField> ff(
            MMFF::constructForceField(tmp, &props, 1.0e8, cid));
        if (ff) {
          ff->initialize();
          ff->minimize(400);
          m->setProp(kEnergyProp,
                     ff->calcEnergy());  // basin E; m's coords unchanged
          m->setProp(kForceFieldProp, labelVariant);
        }
      } catch (...) {
	//  We should have been already FF typed by now, so this should be
	//   unreachable
	// In any case, this is really for diagnostics only so keep the non
	//  basin energy.
      }
    }
  }

  // rank by energy.
  const char *rankProp = d_params.outputRanking == OutputRanking::InterFragScore
                             ? kInterFragScoreProp
                             : kEnergyProp;
  std::sort(result.conformers.begin(), result.conformers.end(),
            [&](const ROMOL_SPTR &a, const ROMOL_SPTR &b) {
              double ea = std::numeric_limits<double>::max();
              double eb = std::numeric_limits<double>::max();
              a->getPropIfPresent<double>(rankProp, ea);
              b->getPropIfPresent<double>(rankProp, eb);
              return ea < eb;
            });
  // -1 AUTO -> derived from the rotor count; 0 DISABLED -> no cap; >0 explicit.
  const int cap =
      isAuto(d_params.numOutputConfs)
          ? autoOutputConfs(result.numLinkBonds +
                            static_cast<unsigned int>(nIntraRotors))
          : d_params.numOutputConfs;
  if (cap > 0 && result.conformers.size() > static_cast<size_t>(cap)) {
    result.conformers.resize(static_cast<size_t>(cap));
  }
}

}  // namespace RDKit
