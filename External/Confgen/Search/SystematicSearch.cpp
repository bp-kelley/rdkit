//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Systematic search
//  generate pools, trim by energy, combine pools at the end
//  prune by rms
#include "Search/SystematicSearch.h"

#include "Search/RotorDriver.h"
#include "Search/InterFragScore.h"
#include "Search/SearchResult.h"
#include "Search/RotorRefine.h"
#include "Utils/NonbondedLookup.h"
#include "Utils/SymmetricRmsd.h"
#include "Joiner/JoinerProfiling.h"
#include "Sampler/TorsionSampler.h"

#include <GraphMol/MolOps.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <ForceField/MMFF/TorsionAngle.h> 
#include <ForceField/MMFF/AngleBend.h>
#include <ForceField/MMFF/StretchBend.h>
#include <ForceField/MMFF/Params.h>
#include <Geometry/point.h>
#include <RDGeneral/RDLog.h>
#include <set>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <atomic>
#include <functional>
#include <limits>
#include <vector>

namespace RDKit {

namespace {

using Vdw = ForceFields::MMFF::InterFragVdWContrib;

inline double confE(const JoinFragment &f, unsigned int c) {
  const double e = f.confs[c].energy;
  return std::isfinite(e) ? e : 0.0;
}

//! Non-superposed coordinate RMSD (over the given flat indices) < thr?
//! for quick rejection
inline bool withinRmsSubset(const std::vector<double> &a,
                            const std::vector<double> &b,
                            const std::vector<unsigned int> &cmpIdx,
                            double thr) {
  const size_t n = cmpIdx.size() / 3;
  const double limit = thr * thr * static_cast<double>(n ? n : 1);
  double s = 0.0;
  for (unsigned int i : cmpIdx) {
    const double d = a[i] - b[i];
    s += d * d;
    if (s >= limit) return false;
  }
  return true;
}

//! One subtree conformer
struct Conf {
  std::vector<double> coords;
  double energy; // FF intergy (can be just junction bonds)
  double fragmentEnergies = 0.0; // energies of joined fragments
};

//! MMFF torsion quartet (idx1-j-k-idx4) and force constants
//! MMFF angle-bend terms
struct AngTerm {
  unsigned int a, b, c;  //!< the angle a-b-c (b central)
  double theta0, ka;     //!< angle-bend
  double sbIJK, sbKJI;   //!< stretch-bend force constants (0 => skip)
  double r0ab, r0cb;     //!< reference bond lengths for the stretch-bend deltas
};

}  // namespace

std::string SystematicSearch::validateParams(
    const RigidRotorSearchParams &sp, const std::string &ffVariant) const {
  const std::string err = RigidRotorSearch::validateParams(sp, ffVariant);
  if (!err.empty()) {
    return err;
  }
  if (sp.systematic.eWindow < 0.0) {
    return "search.systematic.eWindow must be >= 0";
  }
  if (sp.systematic.maxPoolConfs <= 0) {
    return "search.systematic.maxPoolConfs must be > 0";
  }
  if (sp.systematic.combinePoolMultiplier <= 0) {
    return "search.systematic.combinePoolMultiplier must be > 0";
  }
  return {};
}

std::vector<SearchResult> SystematicSearch::search(
    const FragmentJoinerContext &ctx, const RigidRotorSearchParams &sp) {
  const auto &params = sp;
  const size_t nf = ctx.frags.size();
  const unsigned int nAtoms = ctx.mol.getNumAtoms();
  std::vector<SearchResult> out;
  if (nf == 0) return out;
  for (const auto &f : ctx.frags)
    if (f.confs.empty()) return out;

  const SystematicParams &syp = params.systematic;
  const double eWindow = syp.eWindow;
  const size_t maxPoolConfs = static_cast<size_t>(syp.maxPoolConfs);
  // Set the max working memory before we get to the maxPoolConfs
  //   XXX FIX ME -> this is not great and doesn't make much sense
  //                 to have so many of these
  const size_t workingPoolSize =
      std::max<size_t>(maxPoolConfs * std::max(1, syp.combinePoolMultiplier), 2048);
  const size_t maxPool = workingPoolSize * 2;
  const unsigned int maxFragmentConfs = static_cast<unsigned int>(syp.maxFragmentConfs);
  const double finalRms =
      resolveAuto(syp.finalRms, params.diversityRmsThresh > 0.0
                                    ? params.diversityRmsThresh
                                    : 0.5);
  const bool validate = params.diagnostics.SYS_VALIDATE;

  // extract the subtree from the joiner context
  std::vector<std::vector<size_t>> childEdges(nf);
  for (size_t e = 0; e < ctx.edges.size(); ++e)
    childEdges[ctx.edges[e].parentFrag].push_back(e);

  std::vector<std::vector<char>> subtreeAtoms(nf, std::vector<char>(nAtoms, 0));
  
  // Order edges from parent fragment (fixed) to child fragments
  //  search rotations inside out assuming inside edges have larger
  //  effects (i.e. more atoms)
  std::vector<unsigned int> edgeOrder;
  {
    std::function<void(unsigned int)> dfs = [&](unsigned int f) {
      for (size_t e : childEdges[f]) dfs(ctx.edges[e].childFrag);
      edgeOrder.push_back(f);
    };
    dfs(ctx.root);
  }
  for (unsigned int f : edgeOrder) {
    for (unsigned int a : ctx.frags[f].atoms) subtreeAtoms[f][a] = 1;
    for (size_t e : childEdges[f]) {
      const auto &cs = subtreeAtoms[ctx.edges[e].childFrag];
      for (unsigned int a = 0; a < nAtoms; ++a)
        if (cs[a]) subtreeAtoms[f][a] = 1;
    }
  }

  // Compute rigid energy
  //  Get VdW cross pairs
  Vdw *vdw = ctx.scoreHandles.vdw;
  const size_t nPairs = vdw ? vdw->numPairs() : 0;
  std::vector<std::vector<unsigned int>> crossPairs(ctx.edges.size());
  {
    // only compute cross pairs
    for (unsigned int f : edgeOrder) {
      std::vector<char> determined(nAtoms, 0);
      for (unsigned int a : ctx.frags[f].atoms) determined[a] = 1;
      for (size_t e : childEdges[f]) {
        const auto &sc = subtreeAtoms[ctx.edges[e].childFrag];
        for (unsigned int p = 0; p < nPairs; ++p) {
          const int a = vdw->at1(p), b = vdw->at2(p);
          if ((sc[a] && determined[b]) || (sc[b] && determined[a]))
            crossPairs[e].push_back(p);
        }
        for (unsigned int a = 0; a < nAtoms; ++a)
          if (sc[a]) determined[a] = 1;
      }
    }
  }

  // Make the interfrag MMFF score torsions and angles
  std::vector<std::vector<JunctionTorsionTerm>> torTerms(
      ctx.edges.size());
  std::vector<std::vector<AngTerm>> angTerms(ctx.edges.size());
  // NOTE: this search is currently MMFF-ONLY,
  //  it will be dramatically slower for a non-mmff variant so
  //  we don't allow this
  {
    RWMol owned(ctx.mol);
    RDLog::LogStateSetter blk;
    MMFF::MMFFMolProperties props(owned, ctx.ffVariant);
    if (props.isValid()) {
      std::vector<std::pair<unsigned int, unsigned int>> junctionBonds;
      junctionBonds.reserve(ctx.edges.size());
      for (const auto &ed : ctx.edges) {
        junctionBonds.push_back({ed.parentAtom, ed.childAtom});
      }
      torTerms = junctionTorsionTerms(owned, props, junctionBonds);
      // Add angle bend if necessary
      //  XXX FIX ME -> it's nescessary :)
      if (syp.junctionAngleTerms) {
        auto addTriple = [&](size_t e, unsigned int a, unsigned int b,
                             unsigned int c) {
          unsigned int at = 0;
          ForceFields::MMFF::MMFFAngle ap;
          if (!props.getMMFFAngleBendParams(owned, a, b, c, at, ap)) return;
          AngTerm t{a, b, c, ap.theta0, ap.ka, 0.0, 0.0, 0.0, 0.0};
          unsigned int st = 0;
          ForceFields::MMFF::MMFFStbn sp2;
          ForceFields::MMFF::MMFFBond bp[2];
          ForceFields::MMFF::MMFFAngle ap2;
          if (props.getMMFFStretchBendParams(owned, a, b, c, st, sp2, bp,
                                             ap2)) {
            const auto fc =
                ForceFields::MMFF::Utils::calcStbnForceConstants(&sp2);
            t.sbIJK = fc.first;
            t.sbKJI = fc.second;
            t.r0ab = bp[0].r0;
            t.r0cb = bp[1].r0;
          }
          angTerms[e].push_back(t);
        };
        for (size_t e = 0; e < ctx.edges.size(); ++e) {
          const unsigned int j = ctx.edges[e].parentAtom,
                             k = ctx.edges[e].childAtom;
          for (const auto nb : owned.atomNeighbors(owned.getAtomWithIdx(j)))
            if (nb->getIdx() != k) addTriple(e, nb->getIdx(), j, k);
          for (const auto nb : owned.atomNeighbors(owned.getAtomWithIdx(k)))
            if (nb->getIdx() != j) addTriple(e, j, k, nb->getIdx());
        }
      }
    }
  }
  // Junction angle-bend + stretch-bend energy. 
  auto junctionAngles = [&](size_t e, const std::vector<double> &pos) {
    double s = 0.0;
    for (const auto &t : angTerms[e]) {
      auto P = [&](unsigned int a) {
        return RDGeom::Point3D(pos[3 * a], pos[3 * a + 1], pos[3 * a + 2]);
      };
      const RDGeom::Point3D pa = P(t.a), pb = P(t.b), pc = P(t.c);
      RDGeom::Point3D v1 = pa - pb, v2 = pc - pb;
      const double d1 = v1.length(), d2 = v2.length();
      if (d1 < 1e-8 || d2 < 1e-8) continue;
      const double cosT =
          std::max(-1.0, std::min(1.0, v1.dotProduct(v2) / (d1 * d2)));
      s += ForceFields::MMFF::Utils::calcAngleBendEnergy(t.theta0, t.ka, false,
                                                         cosT);
      if (t.sbIJK != 0.0 || t.sbKJI != 0.0) {
        const double dTheta = std::acos(cosT) * 180.0 / M_PI - t.theta0;
        s += ForceFields::MMFF::Utils::calcStretchBendEnergy(
                 d1 - t.r0ab, d2 - t.r0cb, dTheta,
                 std::make_pair(t.sbIJK, t.sbKJI))
                 .first;
      }
    }
    return s;
  };
  auto junctionTorsion = [&](size_t e, const std::vector<double> &pos) {
    double s = 0.0;
    for (const auto &t : torTerms[e]) {
      auto P = [&](unsigned int a) {
        return RDGeom::Point3D(pos[3 * a], pos[3 * a + 1], pos[3 * a + 2]);
      };
      s += ForceFields::MMFF::Utils::calcTorsionEnergy(
          t.V1, t.V2, t.V3,
          ForceFields::MMFF::Utils::calcTorsionCosPhi(P(t.a1), P(t.j), P(t.k),
                                                      P(t.a4)));
    }
    return s;
  };

  // ---- rotor drive
  RotorDriver drv(ctx.mol, ctx.rotorBonds, -1, ctx.scorer);
  std::vector<std::vector<double>> angAt(ctx.edges.size());
  for (size_t e = 0; e < ctx.edges.size(); ++e) {
    std::vector<double> a;
    if (params.torsionSampler) {
      auto t = drv.torsion(static_cast<unsigned int>(e));
      try {
        a = params.torsionSampler->getAngles(
            ctx.mol, t[0], t[1], t[2], t[3],
            useBasinAnglesForRotor(params.junctionBasinAngles, ctx.mol, t[1], t[2]));
      } catch (...) {
      }
    }
    // Set default angles only if we don't find anything
    if (a.empty()) a = params.defaultAngles;
    angAt[e] = std::move(a);
  }

  std::vector<std::vector<Conf>> pool(nf);
  for (unsigned int f = 0; f < nf; ++f) {
    for (size_t c = 0; c < ctx.frags[f].confs.size(); ++c) {
      const double e = confE(ctx.frags[f], static_cast<unsigned int>(c));
      Conf cf;
      cf.coords.resize(3u * nAtoms, 0.0);
      const auto &p = ctx.frags[f].confs[c].pos;
      for (unsigned int a = 0; a < nAtoms; ++a) {
        cf.coords[3 * a] = p[a].x;
        cf.coords[3 * a + 1] = p[a].y;
        cf.coords[3 * a + 2] = p[a].z;
      }
      cf.energy = e;
      cf.fragmentEnergies = e;
      pool[f].push_back(std::move(cf));
      if (pool[f].size() >= maxFragmentConfs) break;
    }
    if (pool[f].empty()) {
      Conf cf;
      cf.coords.assign(3u * nAtoms, 0.0);
      cf.energy = 0.0;
      pool[f].push_back(std::move(cf));
    }
  }

  // prune the bool from lowest energy
  const bool energyOnly = syp.energyOnly;
  const double nodeRms = resolveAuto(syp.nodeRms, finalRms);
  const double dedupBand = syp.dedupEnergyBand;
  auto prune = [&](std::vector<Conf> &v,
                   const std::vector<unsigned int> &cmpIdx) {
    if (v.empty()) return;
    std::sort(v.begin(), v.end(),
              [](const Conf &a, const Conf &b) { return a.energy < b.energy; });
    const double lim = v.front().energy + eWindow;
    std::vector<Conf> keep;
    for (auto &c : v) {
      if (std::isfinite(c.energy) && c.energy > lim) {
        if (detail::profiling())
          detail::prof().nRejWindow +=
              static_cast<long long>(&v.back() - &c) + 1;
        break;
      }
      if (!energyOnly && nodeRms > 0.0) {
	// Scan from the back generally finds the dupes sooner.
	//  this is limited by a minimum energy to stop at.
	// n.b. this isn't really proven equivalent, but runs a lot faster and
	//  doesn't seem to change the result
        bool dup = false;
        for (auto it = keep.rbegin(); it != keep.rend(); ++it) {
          // Only scan in the energy band for speed
          if (dedupBand > 0.0 && std::isfinite(c.energy) &&
              std::isfinite(it->energy) &&
              (c.energy - it->energy) > dedupBand) {
            break;
          }
          if (withinRmsSubset(c.coords, it->coords, cmpIdx, nodeRms)) {
            dup = true;
            break;
          }
        }
        if (dup) {
          if (detail::profiling()) detail::prof().nRejRms += 1;
          continue;
        }
      }
      keep.push_back(std::move(c));
      if (keep.size() >= maxPoolConfs) {
        if (detail::profiling()) detail::prof().nRejCap += 1;
        break;
      }
    }
    v.swap(keep);
  };

  // Search each pool now and combine (most expensive step)
  //  XXX FIX ME -> why doesn't this use rotortree?
  const auto deadline = timeOut(params);
  for (unsigned int f : edgeOrder) {
    if (timedOut(deadline)) break;
    std::vector<Conf> cur = std::move(pool[f]);

    std::vector<char> coords_changed(nAtoms, 0);
    for (unsigned int a : ctx.frags[f].atoms) coords_changed[a] = 1;
    for (size_t e : childEdges[f]) {
      if (timedOut(deadline)) break;
      const unsigned int cf = ctx.edges[e].childFrag;
      const unsigned int pa = ctx.edges[e].parentAtom,
                         ca = ctx.edges[e].childAtom;
      const double bondLen = ctx.edges[e].bondLen;
      const auto &sc = subtreeAtoms[cf];
      std::vector<Conf> next;
      double runningMin = std::numeric_limits<double>::infinity();
      // Only drive things that move, eveything else is constant
      std::vector<unsigned int> scIdx;
      scIdx.reserve(nAtoms);
      for (unsigned int a = 0; a < nAtoms; ++a) {
        if (sc[a]) scIdx.push_back(a);
      }
      for (const auto &pc : cur) {
        // parent coords never change when rotating attached frags
        const RDGeom::Point3D parentBondAtom(
            pc.coords[3 * pa], pc.coords[3 * pa + 1], pc.coords[3 * pa + 2]);
        const RDGeom::Point3D parentExit(
            pc.coords[3 * ca], pc.coords[3 * ca + 1], pc.coords[3 * ca + 2]);
        drv.positions() = pc.coords;  // parent atoms are constant across this
                                      // parent's children
        for (const auto &cc : pool[cf]) {
          // join child subtree onto the parent at the junction (rigid; dihedral
          // arbitrary)
          const ChildPlacement pl = computeChildPlacement(
              RDGeom::Point3D(cc.coords[3 * ca], cc.coords[3 * ca + 1],
                              cc.coords[3 * ca + 2]),
              RDGeom::Point3D(cc.coords[3 * pa], cc.coords[3 * pa + 1],
                              cc.coords[3 * pa + 2]),
              parentBondAtom, parentExit, bondLen);
          for (unsigned int a : scIdx) {
            const RDGeom::Point3D q = pl.apply(RDGeom::Point3D(
                cc.coords[3 * a], cc.coords[3 * a + 1], cc.coords[3 * a + 2]));
            drv.positions()[3 * a] = q.x;
            drv.positions()[3 * a + 1] = q.y;
            drv.positions()[3 * a + 2] = q.z;
          }
          const double baseE = pc.energy + cc.energy;
          for (const double ang : angAt[e]) {
            drv.setDihedral(static_cast<unsigned int>(e),
                            ang);  // rotates the child subtree
            const auto &buf = drv.positions();
            const double eXtra =
                (vdw ? vdw->energyOfPairs(buf.data(), crossPairs[e]) : 0.0) +
                junctionTorsion(e, buf) +
                (syp.junctionAngleTerms ? junctionAngles(e, buf) : 0.0);
            const double energy = baseE + eXtra;
            if (detail::profiling()) detail::prof().nScoreInterFrag += 1;
            if (energy < runningMin) runningMin = energy;
            // Remove conformers above our ewindow (plus a little slack
	    //  before combining)
	    // XXX FIX ME -> I don't think we need the slack...
            if (energy > runningMin + eWindow + syp.upperEnergyWindow) {
              if (detail::profiling()) detail::prof().nRejQuick += 1;
              continue;
            }
            next.push_back({buf, energy, pc.fragmentEnergies + cc.fragmentEnergies});
            if (detail::profiling()) detail::prof().nPoolMade += 1;
          }
          // If we've reached our upper size bound, prune
          if (next.size() > maxPool) {
            std::nth_element(next.begin(), next.begin() + workingPoolSize, next.end(),
                             [](const Conf &a, const Conf &b) {
                               return a.energy < b.energy;
                             });
            if (detail::profiling())
              detail::prof().nRejTrim += next.size() - workingPoolSize;
            next.resize(workingPoolSize);
            if (params.diagnostics.SYS_VALIDATE) { 
              static std::atomic<long> nTrim{0};
              BOOST_LOG(rdWarningLog)
                  << "[sysTRIM] fired #" << ++nTrim
                  << " (workingPoolSize=" << workingPoolSize
                  << " maxPoolConfs=" << maxPoolConfs << ")\n";
            }
          }
        }
      }
      for (unsigned int a = 0; a < nAtoms; ++a)
        if (sc[a]) coords_changed[a] = 1;
      std::vector<unsigned int> cmpIdx;
      cmpIdx.reserve(3u * nAtoms);
      for (unsigned int a = 0; a < nAtoms; ++a)
        if (coords_changed[a]) {
          cmpIdx.push_back(3 * a);
          cmpIdx.push_back(3 * a + 1);
          cmpIdx.push_back(3 * a + 2);
        }
      {
        const long long tp0 = detail::profiling() ? detail::nowNs() : 0;
        prune(next, cmpIdx);
        if (detail::profiling()) {
          detail::prof().tPrune += detail::nowNs() - tp0;
          detail::prof().nPoolKept += static_cast<long long>(next.size());
        }
      }
      cur.swap(next);
    }
    pool[f] = std::move(cur);
  }

  std::vector<Conf> &rootPool = pool[ctx.root];

  // If we want to twiddle torsions in the fragments themselves,
  //   here is were we do it.
  // XXX FIX ME -> this code is complicated and possibly unecessary
  if (params.driveIntraFragmentTorsions && !ctx.intraRotorBonds.empty() &&
      !rootPool.empty()) {
    RotorDriver idrv(ctx.mol, ctx.intraRotorBonds, -1, ctx.scorer);
    const size_t nIntra = ctx.intraRotorBonds.size();
    // candidate angles per torsions between fragments
    std::vector<std::vector<double>> iang(nIntra);
    for (size_t r = 0; r < nIntra; ++r) {
      std::vector<double> a;
      if (params.torsionSampler) {
        try {
          const auto t = idrv.torsion(static_cast<unsigned int>(r));
          a = params.torsionSampler->getAngles(ctx.mol, t[0], t[1], t[2], t[3],
                                               true);
        } catch (...) {
        }
      }
      if (a.empty()) a = params.defaultAngles;
      iang[r] = std::move(a);
    }
    std::sort(rootPool.begin(), rootPool.end(),
              [](const Conf &a, const Conf &b) { return a.energy < b.energy; });
    const size_t nSeed = std::min<size_t>(rootPool.size(), maxPoolConfs);
    std::vector<Conf> extra;
    for (size_t s = 0; s < nSeed; ++s) {
      const Conf &seed = rootPool[s];
      for (size_t r = 0; r < nIntra; ++r) {
        idrv.positions() = seed.coords;
        const double cur = idrv.dihedralDeg(static_cast<unsigned int>(r));
        for (const double ang : iang[r]) {
          if (std::fabs(std::fmod(ang - cur + 540.0, 360.0) - 180.0) < 15.0)
            continue;  // ~= seed
          idrv.setDihedral(static_cast<unsigned int>(r), ang);
          const double e = ctx.scorer
                               ? ctx.scorer(idrv.positions().data(), nAtoms)
                               : seed.energy;
          if (detail::profiling()) detail::prof().nScoreFullMMFF += 1;
          extra.push_back({idrv.positions(), e + seed.fragmentEnergies, seed.fragmentEnergies});
          idrv.positions() =
              seed.coords;  // re-seed: vary ONE intra bond at a time
        }
      }
      if (extra.size() > maxPoolConfs) break;
    }
    
    if (!extra.empty()) {
      rootPool.insert(rootPool.end(), std::make_move_iterator(extra.begin()),
                      std::make_move_iterator(extra.end()));
      // whole-molecule comparison set: by the root every atom is placed
      std::vector<unsigned int> rootCmp;
      rootCmp.reserve(3u * nAtoms);
      for (unsigned int a = 0; a < nAtoms; ++a) {
        rootCmp.push_back(3 * a);
        rootCmp.push_back(3 * a + 1);
        rootCmp.push_back(3 * a + 2);
      }
      prune(rootPool, rootCmp);
    }
  }
  // - end of fragment twiddling

  // Diagnose our energies versus a full MMFF
  // XXX Fix me -> is this USED?  Might be worth removing.
  if (validate && !rootPool.empty() && ctx.scorer) {
    const Conf *best = &rootPool.front();
    for (const auto &c : rootPool)
      if (c.energy < best->energy) best = &c;
    const double full = ctx.scorer(best->coords.data(), nAtoms);
    if (detail::profiling()) detail::prof().nScoreFullMMFF += 1;
    //  Hope we have no residuals
    const double bottomUpInter = best->energy - best->fragmentEnergies;
    const double resid = bottomUpInter - full;
    static std::atomic<long> nChk{0}, nBad{0};
    const long i = ++nChk;
    const bool bad = std::fabs(resid) > 1e-3 * std::max(1.0, std::fabs(full));
    if (bad) ++nBad;
    // A partition mismatch means the decomposed energy disagrees with the
    // full scorer -- a bug, not a diagnostic, hence rdErrorLog.
    auto &log = bad ? rdErrorLog : rdWarningLog;
    BOOST_LOG(log) << "[SYS_VALIDATE] #" << i << " bottomUpInter="
                   << bottomUpInter << " fullInter=" << full
                   << " resid=" << resid
                   << (bad ? "  <== PARTITION MISMATCH" : "")
                   << " (intra=" << best->fragmentEnergies
                   << " total=" << best->energy << ") badSoFar=" << nBad.load()
                   << "\n";
    std::fflush(stderr);
  }

  // Final deduplication and sorting
  std::vector<SearchResult> res;
  res.reserve(rootPool.size());
  for (auto &c : rootPool) res.push_back({std::move(c.coords), c.energy});
  std::sort(res.begin(), res.end(),
            [](const SearchResult &a, const SearchResult &b) {
              if (std::isnan(a.score)) return false;
              if (std::isnan(b.score)) return true;
              return a.score < b.score;
            });
  // Local coordinate descent around whatever basins were chosen
  if (params.thompson.refineSteps > 0) {
    refineRotorsInPlace(drv, res, params.thompson.refineSteps,
                        params.thompson.refineStepDeg,
                        params.thompson.refinePasses);
  }
  if (finalRms > 0.0) ctx.symmetryDedupInPlace(res, finalRms);
  return res;
}

}  // namespace RDKit
