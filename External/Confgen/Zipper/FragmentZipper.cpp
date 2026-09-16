//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Zipper/FragmentZipper.h"
#include "FragmentConfGen.h"
#include "Utils/TheobaldRmsd.h"
#include "Zipper/ZipperProfiling.h"

#include <Geometry/Transform3D.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include "Embedder/Embedder.h"
#include "Utils/SymmetricRmsd.h"
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <ForceField/ForceField.h>
#include <ForceField/PositionConstraint.h>
#include <RDGeneral/types.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <utility>

namespace {

//! Lowest finite energy (relative to minimum) in a fragment pool, or 0.0 when the pool carries no
//! energies.
inline double poolEnergyBase(const std::vector<double> &energies) {
  double base = std::numeric_limits<double>::max();
  bool any = false;
  for (double e : energies) {
    if (std::isfinite(e)) {
      base = std::min(base, e);
      any = true;
    }
  }
  return any ? base : 0.0;
}

//! Energy of pool entry 'ci' to base energy
inline double relativePoolEnergy(const std::vector<double> &energies, size_t ci,
                                 double eBase) {
  if (ci < energies.size() && std::isfinite(energies[ci]))
    return energies[ci] - eBase;
  return 0.0;
}
}  // namespace

namespace RDKit {

namespace {
//! printf-formatted line onto the RDKit warning log.
/*!
  XXX FIX ME - we should probably turn disable -Wno-format here explicitly
                RDKit builds enable it
*/
#if defined(__GNUC__) || defined(__clang__)
#define CONFGEN_PRINTF_FMT __attribute__((format(printf, 1, 2)))
#else
#define CONFGEN_PRINTF_FMT
#endif

// The attribute must TRAIL the declarator; written before the return type it is
// accepted and silently ignored, which costs the very checking it is here for.
void logfWarn(const char *fmt, ...) CONFGEN_PRINTF_FMT;

void logfWarn(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  va_list ap2;
  va_copy(ap2, ap);
  const int n = std::vsnprintf(nullptr, 0, fmt, ap);
  va_end(ap);
  if (n < 0) {
    va_end(ap2);
    return;
  }
  std::string buf(static_cast<size_t>(n) + 1, '\0');
  std::vsnprintf(buf.data(), buf.size(), fmt, ap2);
  va_end(ap2);
  buf.resize(static_cast<size_t>(n));
  BOOST_LOG(rdWarningLog) << buf;
}
}  // namespace


namespace {
// Enabled by setJoinerProfiling() (default off);
bool &profilingFlag() {
  static bool b = false;
  return b;
}
}  // namespace

namespace detail {
JoinProf &prof() {
  static JoinProf p;
  return p;
}
bool profiling() { return profilingFlag(); }
void setJoinerProfiling(bool on) { profilingFlag() = on; }
}  // namespace detail

using detail::nowNs;
using detail::prof;
using detail::profiling;

namespace {
//! Place every fragment and get the resulting coords
//! PRECONDITION:
//! all frags have confs and confChoice indexes into each fragment's pool.
std::vector<double> placeFragments(const ROMol &mol,
                                   const std::vector<ZipFragment> &frags,
                                   const std::vector<ZipEdge> &edges,
                                   unsigned int root,
                                   const std::vector<unsigned int> &confChoice);
}  // namespace

void printJoinerProfile(const char *label, size_t nMols) {
  auto &p = prof();
  auto ms = [](long long ns) { return ns / 1.0e6; };
  double total = ms(p.tEmbed + p.tPlace + p.tDrive + p.tScore + p.tDiv);
  // Score the appropriate search timing detail for the given search
  if (total <= 0.0 && p.nScoreInterFrag.load() == 0 && p.nPoolMade.load() == 0)
    return;
  if (total <= 0.0)
    total = 1e-9;  // avoid divide-by-zero in the percentage columns
  double warm = ms(p.tPlace + p.tDrive + p.tScore + p.tDiv);  // no embed
  auto pct = [&](long long ns) { return 100.0 * ms(ns) / total; };
  logfWarn("[asm-profile %s] total=%.0fms  embed=%.0f(%.0f%%) "
               "place=%.0f(%.0f%%) drive=%.0f(%.0f%%) score=%.0f(%.0f%%) "
               "div=%.0f(%.0f%%)\n",
               label, total, ms(p.tEmbed), pct(p.tEmbed), ms(p.tPlace),
               pct(p.tPlace), ms(p.tDrive), pct(p.tDrive), ms(p.tScore),
               pct(p.tScore), ms(p.tDiv), pct(p.tDiv));
  logfWarn("               embeds=%lld samples=%lld drives=%lld scores=%lld "
      "kept=%lld out=%lld\n",
      p.nEmbeds.load(), p.nSamples.load(), p.nDrives.load(), p.nScores.load(),
      p.nKept.load(), p.nOut.load());
  if (nMols > 0) {
    logfWarn("               cold(with embed)=%.0f ms/mol   WARM(assembly "
                 "only, cached frags)=%.0f ms/mol\n",
                 total / nMols, warm / nMols);
  }

  // break down timings of individual class embeddings
  static const char *kClassName[detail::JoinProf::kNClass] = {
      "Rigid", "SmallRing", "LargeRing", "Acyclic", "Exhaustive", "Fast"};
  bool anyClass = false;
  long long classNs = 0;
  for (int i = 0; i < detail::JoinProf::kNClass; ++i) {
    anyClass |= p.nEmbedClass[i].load() > 0;
    classNs += p.tEmbedClass[i].load();
  }
  const double denom = std::max(total, ms(classNs));
  if (anyClass) {
    logfWarn("               EMBED by class (ms, n frags, n confs kept):\n");
    for (int i = 0; i < detail::JoinProf::kNClass; ++i) {
      const long long n = p.nEmbedClass[i].load();
      if (!n) continue;
      const double cms = ms(p.tEmbedClass[i].load());
      logfWarn("                 %-11s %8.0f ms (%4.1f%% of total)  frags=%-7lld "
          "confs=%-8lld  %.2f ms/frag  %.1f confs/frag\n",
          kClassName[i], cms, 100.0 * cms / denom, n, p.nConfClass[i].load(),
          cms / n, static_cast<double>(p.nConfClass[i].load()) / n);
    }
  }
  const long long si = p.nScoreInterFrag.load(), sf = p.nScoreFullMMFF.load();
  const long long pm = p.nPoolMade.load(), pk = p.nPoolKept.load();
  if (si || sf || pm) {
    logfWarn("               SCORES: interFrag(vdW+junctionTorsion)=%lld  fullMMFF=%lld"
        "%s\n",
        si, sf, nMols ? "" : "");
    if (nMols > 0 && si) {
      logfWarn("                       = %.0f interFrag scores/mol\n",
                   static_cast<double>(si) / nMols);
    }
    logfWarn("               POOL: materialised=%lld kept=%lld (%.1f%% survive)  "
        "prune=%.0f ms  dedup(div)=%.0f ms\n",
        pm, pk, pm ? 100.0 * pk / pm : 0.0, ms(p.tPrune.load()), ms(p.tDiv));
    const long long rq = p.nRejQuick.load(), rw = p.nRejWindow.load(),
                    rr = p.nRejRms.load(), rc = p.nRejCap.load(),
                    rt = p.nRejTrim.load();
    const long long rtot = rq + rw + rr + rc + rt;
    if (rtot > 0) {
      auto rp = [&](long long v) { return 100.0 * v / rtot; };
      logfWarn("               WHY REJECTED (of %lld): eWindow-quick=%lld(%.1f%%) "
          "eWindow-final=%lld(%.1f%%) rmsDup=%lld(%.1f%%) poolCap=%lld(%.1f%%) "
          "memTrim=%lld(%.1f%%)\n",
          rtot, rq, rp(rq), rw, rp(rw), rr, rp(rr), rc, rp(rc), rt, rp(rt));
    }
  }

  p.tEmbed = p.tPlace = p.tDrive = p.tScore = p.tDiv = 0;
  p.nEmbeds = p.nSamples = p.nDrives = p.nScores = p.nKept = p.nOut = 0;
  p.nScoreInterFrag = p.nScoreFullMMFF = p.nPoolMade = p.nPoolKept = 0;
  p.nRejQuick = p.nRejWindow = p.nRejRms = p.nRejCap = p.nRejTrim = 0;
  p.tCombine = p.tPrune = 0;
  for (int i = 0; i < detail::JoinProf::kNClass; ++i) {
    p.tEmbedClass[i] = 0;
    p.nEmbedClass[i] = 0;
    p.nConfClass[i] = 0;
  }
}

void setJoinerProfiling(bool on) { profilingFlag() = on; }

double takeJoinerWarmMsPerMol(size_t nMols) {
  auto &p = prof();
  // WARM = assume frag is in cache and don't time.
  double warm = (p.tPlace + p.tDrive + p.tScore + p.tDiv) / 1.0e6;
  p.tEmbed = p.tPlace = p.tDrive = p.tScore = p.tDiv = 0;
  p.nEmbeds = p.nSamples = p.nDrives = p.nScores = p.nKept = p.nOut = 0;
  return nMols ? warm / static_cast<double>(nMols) : warm;
}

ChildPlacement computeChildPlacement(const RDGeom::Point3D &childBondPos,
                                     const RDGeom::Point3D &childExitPos,
                                     const RDGeom::Point3D &parentBondAtom,
                                     const RDGeom::Point3D &parentExit,
                                     double bondLen) {
  RDGeom::Point3D uParent = parentExit - parentBondAtom;
  uParent.normalize();

  const RDGeom::Point3D childNbr = childBondPos;
  RDGeom::Point3D vChild = childExitPos - childNbr;
  vChild.normalize();

  // the child's exit points to the parent
  RDGeom::Point3D target = uParent * -1.0;

  RDGeom::Transform3D rot;  // defaults to I
  double cosT = std::max(-1.0, std::min(1.0, vChild.dotProduct(target)));
  if (cosT < 1.0 - 1e-6) {
    RDGeom::Point3D axis;
    if (cosT < -1.0 + 1e-6) {
      // antiparallel: any axis perpendicular to the child works for the 180 flip
      axis = vChild.crossProduct(RDGeom::Point3D(1.0, 0.0, 0.0));
      if (axis.lengthSq() < 1e-6) {
        axis = vChild.crossProduct(RDGeom::Point3D(0.0, 1.0, 0.0));
      }
      axis.normalize();
      rot.SetRotation(M_PI, axis);
    } else {
      axis = vChild.crossProduct(target);
      axis.normalize();
      rot.SetRotation(std::acos(cosT), axis);
    }
  }

  // rotate the child about its bonding atom, then translate that atom onto the
  // new junction site at the MMFF-ideal bond length along the parent exit
  // direction
  const RDGeom::Point3D site = parentBondAtom + uParent * bondLen;
  return ChildPlacement{rot, childNbr, site};
}

void placeChildCoords(std::vector<RDGeom::Point3D> &childCoords,
                      unsigned int childBondAtom, unsigned int childExit,
                      const RDGeom::Point3D &parentBondAtom,
                      const RDGeom::Point3D &parentExit, double bondLen) {
  const ChildPlacement pl =
      computeChildPlacement(childCoords[childBondAtom], childCoords[childExit],
                            parentBondAtom, parentExit, bondLen);
  for (auto &p : childCoords) {
    p = pl.apply(p);
  }
}

namespace {
double covBondLen(unsigned int z1, unsigned int z2) {
  const auto *pt = PeriodicTable::getTable();
  double l = pt->getRcovalent(z1) + pt->getRcovalent(z2);
  return l > 0.1 ? l : 1.5;
}

// embed nConfs
//  note:  when not using MMFF minimization, shrugDisplacement
//   can be run to try to shrug off non mmff ideal bond lengths
//   and energies.
std::vector<double> embedMin(RWMol &m, unsigned int nConfs, int seed,
                             const std::string &variant, bool useDG = false,
                             bool minimizeMMFF = true, double shrugDispl = 0.0,
                             double shrugK = 100.0) {
  long long t0 = profiling() ? nowNs() : 0;

  MolOps::addHs(m);

  DGeomHelpers::EmbedParameters ps;
  if (!useDG) ps = DGeomHelpers::ETKDGv3;
  ps.randomSeed = seed;
  INT_VECT cids;
  DGeomHelpers::EmbedMultipleConfs(m, cids, nConfs, ps);
  std::vector<double> es(m.getNumConformers(),
                         std::numeric_limits<double>::quiet_NaN());
  MMFF::MMFFMolProperties props(m, variant);
  if (props.isValid()) {
    props.setMMFFEleTerm(false);
    // shrugDispl > 0: use flat bottom wells to help shrug off energies
    for (int cid : cids) {
      try {
        std::unique_ptr<ForceFields::ForceField> ff(
            MMFF::constructForceField(m, &props, 1.0e8, cid));
        if (ff) {
          if (shrugDispl > 0.0) {
            for (unsigned int a = 0; a < m.getNumAtoms(); ++a) {
              if (m.getAtomWithIdx(a)->getAtomicNum() != 1) {
                ff->contribs().push_back(ForceFields::ContribPtr(
                    new ForceFields::PositionConstraintContrib(
                        ff.get(), a, shrugDispl, shrugK)));
              }
            }
          }
          ff->initialize();
          if (minimizeMMFF || shrugDispl > 0.0)
            ff->minimize(200);         // constrained = shrug
          es[cid] = ff->calcEnergy();  // single-point energy for sorting
        }
      } catch (...) {
        // degenerate embed?
      }
    }
  }
  if (profiling()) {
    prof().tEmbed += nowNs() - t0;
    prof().nEmbeds += 1;
  }
  return es;
}
}  // namespace

FragmentZipperInput buildFragmentJoinerInput(
    const ROMol &input, unsigned int nConfs, int seed,
    const std::string &ffVariant, bool fragUseDG, bool fragMinimizeFF,
    const Embedder *lib, const FragmentZipperParams *asmParams,
    const std::vector<unsigned int> *linkBondsOverride) {
  FragmentZipperInput out;

  const double shrugDispl = asmParams ? asmParams->fragShrugDisplacement : 0.0;
  const double shrugK = asmParams ? asmParams->fragShrugForceConst : 100.0;
  const bool exactGeom = asmParams && asmParams->diagnostics.ASM_EXACT_GEOM;

  // largest component, with explicit Hs.
  RWMOL_SPTR_VECT comps = getRWMolFrags(input, /*sanitizeFrags=*/true);
  RWMOL_SPTR big;
  unsigned int bestH = 0;
  for (const auto &c : comps) {
    unsigned int h = 0;
    for (const auto a : c->atoms()) {
      if (a->getAtomicNum() > 1) ++h;
    }
    if (h >= bestH) {
      bestH = h;
      big = c;
    }
  }
  RWMOL_SPTR m = big ? big : boost::make_shared<RWMol>(input);

  MolOps::addHs(*m);
  out.mol = m;
  const unsigned int n = m->getNumAtoms();

  if (!m->getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(*m);
  }
  const auto rb = FragmentConfGen::findRotatableBonds(
      *m, asmParams && asmParams->sampleTrivialRotors,
      asmParams && asmParams->wholeAcyclicFragments);
  auto links = rb.inter;
  out.intraRotorBonds = rb.intra;  // uncut bonds, each carrying its own rule
  if (linkBondsOverride) {
    // COARSE mode: the caller knows where the real seams are (e.g. the synthon
    // junctions a molzip just made) and wants the pieces between them kept
    // WHOLE.  Every rotatable bond we are told not to cut becomes an
    // intra-fragment rotor instead of disappearing, so the flexibility is
    // still accounted for -- it is just not searched by the junction driver.
    const std::set<unsigned int> keep(linkBondsOverride->begin(),
                                      linkBondsOverride->end());
    for (unsigned int b : links) {
      if (!keep.count(b)) {
        out.intraRotorBonds.push_back({b, IntraRotorType::Free});
      }
    }
    links.assign(keep.begin(), keep.end());
  }

  // fragment partition (union-find over the non-link bonds)
  std::vector<char> isLink(m->getNumBonds(), 0);
  for (auto b : links) isLink[b] = 1;
  std::vector<int> uf(n);
  std::iota(uf.begin(), uf.end(), 0);
  std::function<int(int)> find = [&](int x) {
    while (uf[x] != x) {
      uf[x] = uf[uf[x]];
      x = uf[x];
    }
    return x;
  };
  for (const auto b : m->bonds()) {
    if (!isLink[b->getIdx()]) {
      int ra = find(b->getBeginAtomIdx()), rb = find(b->getEndAtomIdx());
      if (ra != rb) uf[ra] = rb;
    }
  }
  std::map<int, int> fid;
  out.atomFragment.assign(n, 0);
  for (unsigned int a = 0; a < n; ++a) {
    int r = find(a);
    auto it = fid.find(r);
    int id = (it == fid.end()) ? (fid[r] = static_cast<int>(fid.size()))
                               : it->second;
    out.atomFragment[a] = id;
  }
  out.fragments.resize(fid.size());
  for (unsigned int a = 0; a < n; ++a) {
    out.fragments[out.atomFragment[a]].atoms.push_back(a);
  }

  auto sortConfs = [](std::vector<ZipFragment> &fr) {
    for (auto &f : fr) {
      std::sort(
          f.confs.begin(), f.confs.end(),
          [](const JoinFragmentConf &a, const JoinFragmentConf &b) {
            if (std::isnan(a.energy)) return false;
            if (std::isnan(b.energy)) return true;
            return a.energy < b.energy;
          });
    }
  };

  if (links.empty()) {
    if (lib) {
      // we just have one frag so the embedding is all we need
      for (auto a : m->atoms()) {
        a->setProp<int>("_asmMolIdx", static_cast<int>(a->getIdx()));
      }
      std::vector<std::vector<RDGeom::Point3D>> bufs;
      std::vector<double> energies;
      if (!lib->getConformerCoords(*m, n, "_asmMolIdx", bufs, &energies))
        return out;
      // `energy` is always an energy in kcal/mol; position in `confs` carries
      // the ordering.
      const double eBase = poolEnergyBase(energies);
      for (size_t ci = 0; ci < bufs.size(); ++ci) {
        JoinFragmentConf conf;
        conf.pos = std::move(bufs[ci]);
        conf.energy = relativePoolEnergy(energies, ci, eBase);
        out.fragments[0].confs.push_back(std::move(conf));
      }
    } else {
      RWMol pc(*m);
      auto es = embedMin(pc, nConfs, seed, ffVariant, fragUseDG,
                         fragMinimizeFF, shrugDispl, shrugK);
      for (unsigned int cid = 0; cid < pc.getNumConformers(); ++cid) {
        JoinFragmentConf conf;
        conf.pos.resize(n);
        conf.energy = es[cid];
        const Conformer &cf = pc.getConformer(cid);
        for (unsigned int a = 0; a < n; ++a) conf.pos[a] = cf.getAtomPos(a);
        out.fragments[0].confs.push_back(std::move(conf));
      }
    }
    sortConfs(out.fragments);
    return out;
  }

  // Use the MMFF94 reference bond length (r0) for the junction bonds, falling
  // back to the covalent-radius sum below when MMFF cannot type the molecule.
  std::unique_ptr<MMFF::MMFFMolProperties> mmffProps;
  try {
    mmffProps.reset(new MMFF::MMFFMolProperties(*m, "MMFF94"));
    if (!mmffProps->isValid()) mmffProps.reset();
  } catch (...) {
    mmffProps.reset();
  }
  // junctions
  for (auto bi : links) {
    const Bond *b = m->getBondWithIdx(bi);
    unsigned int a = b->getBeginAtomIdx(), c = b->getEndAtomIdx();
    ZipJunction J;
    J.fragA = out.atomFragment[a];
    J.fragB = out.atomFragment[c];
    J.atomA = a;
    J.atomB = c;

    unsigned int bondType = 0;
    ForceFields::MMFF::MMFFBond bp;
    if (mmffProps &&
        mmffProps->getMMFFBondStretchParams(*m, a, c, bondType, bp) &&
        bp.r0 > 0.1) {
      J.bondLen = bp.r0;  // MMFF ideal length for this atom-type/bond-type pair
    } else {
      // otherwise covalent bond length
      J.bondLen = covBondLen(m->getAtomWithIdx(a)->getAtomicNum(),
                             m->getAtomWithIdx(c)->getAtomicNum());
    }
    out.junctions.push_back(J);
  }

  // cut with isotope-labelled dummies (one label per cut) and a tracked mapping
  const unsigned int LB = 1000;
  std::vector<std::pair<unsigned int, unsigned int>> dummyLabels;
  for (size_t i = 0; i < links.size(); ++i)
    dummyLabels.emplace_back(LB + i, LB + i);
  std::unique_ptr<ROMol> fragged(
      MolFragmenter::fragmentOnBonds(*m, links, true, &dummyLabels));
  if (!fragged) return out;
  std::vector<int> flabels;
  std::vector<std::vector<int>> mapping;
  RWMOL_SPTR_VECT pieces =
      getRWMolFrags(*fragged, /*sanitizeFrags=*/true, &flabels, &mapping,
                    /*copyConformers=*/false);

  for (size_t p = 0; p < pieces.size(); ++p) {
    RWMol &pc = *pieces[p];  // mutate the fresh fragment in place (no copy)
    struct Dinfo {
      unsigned int pieceAtom, partner;
    };
    // This is the complicated bookkeeping bit.
    // Detect exit dummies + this piece's fragment id, and tag every atom with
    // its destination index in the full molecule: real atoms -> mapping[p][j];
    // exit dummies -> the partner atom across the cut (its exit-vector
    // stand-in).
    std::vector<Dinfo> dummies;
    int fp = -1;
    for (unsigned int j = 0; j < pc.getNumAtoms(); ++j) {
      Atom *at = pc.getAtomWithIdx(j);
      if (at->getAtomicNum() == 0) {
        unsigned int cut = at->getIsotope() - LB;
        unsigned int nb = (*pc.atomNeighbors(at).begin())->getIdx();
        int nbHeavy = mapping[p][nb];
        const Bond *lb = m->getBondWithIdx(links[cut]);
        unsigned int a = lb->getBeginAtomIdx(), c = lb->getEndAtomIdx();
        unsigned int partner = (nbHeavy == static_cast<int>(a)) ? c : a;
        dummies.push_back({j, partner});
        at->setProp<int>("_asmMolIdx", static_cast<int>(partner));
      } else {
        int mi = mapping[p][j];
        if (mi < static_cast<int>(n)) {
          at->setProp<int>("_asmMolIdx", mi);
          if (fp < 0) fp = out.atomFragment[mi];
        }
      }
    }

    if (exactGeom) {
      // XXX FIX ME slated for removal
      //  use the input's geometries so we can isolate the sampler and torsion
      //  drivers
      const Conformer &mc = m->getConformer();
      JoinFragmentConf conf;
      conf.pos.resize(n);
      conf.energy = 0.0;
      for (unsigned int j = 0; j < mapping[p].size(); ++j) {
        int fidx = mapping[p][j];
        if (fidx < static_cast<int>(n)) conf.pos[fidx] = mc.getAtomPos(fidx);
      }
      for (const auto &d : dummies)
        conf.pos[d.partner] = mc.getAtomPos(d.partner);
      if (fp >= 0) out.fragments[fp].confs.push_back(std::move(conf));
    } else if (lib) {
      std::vector<std::vector<RDGeom::Point3D>> bufs;
      std::vector<double> energies;
      if (!lib->getConformerCoords(pc, n, "_asmMolIdx", bufs, &energies))
        continue;
      // Use the fragment's REAL MMFF energy.  This slot feeds the assembled
      // node energy
      // (`leftE + rightE + crossTerms`
      const double eBase = poolEnergyBase(energies);
      for (size_t ci = 0; ci < bufs.size(); ++ci) {
        JoinFragmentConf conf;
        conf.pos = std::move(bufs[ci]);
        conf.energy = relativePoolEnergy(energies, ci, eBase);
        if (fp >= 0) out.fragments[fp].confs.push_back(std::move(conf));
      }
    } else {
      // XXX FIX ME - we should never really get here
      for (const auto &d : dummies) {
        Atom *at = pc.getAtomWithIdx(d.pieceAtom);
        at->setAtomicNum(6);  // dummies always carbon-capped
        at->setIsotope(0);
        at->setFormalCharge(0);
        at->setNoImplicit(false);
      }
      MolOps::sanitizeMol(pc);
      auto es = embedMin(pc, nConfs, seed, ffVariant, fragUseDG,
                         fragMinimizeFF, shrugDispl, shrugK);
      for (unsigned int cid = 0; cid < pc.getNumConformers(); ++cid) {
        JoinFragmentConf conf;
        conf.pos.resize(n);
        conf.energy = es[cid];
        const Conformer &cf = pc.getConformer(cid);
        for (unsigned int j = 0; j < mapping[p].size(); ++j) {
          int fidx = mapping[p][j];
          if (fidx < static_cast<int>(n)) conf.pos[fidx] = cf.getAtomPos(j);
        }
        for (const auto &d : dummies) {
          conf.pos[d.partner] = cf.getAtomPos(d.pieceAtom);
        }
        if (fp >= 0) out.fragments[fp].confs.push_back(std::move(conf));
      }
    }
  }
  sortConfs(out.fragments);
  return out;
}

std::string FragmentZipperParams::validate() const {
  if (!isValidFF(ffVariant)) {
    // isValidFF() is the same test; kept local so the zipper does not depend on
    // FragmentConfGen.  Widen BOTH when a non-MMFF family arrives.
    return "zipper.ffVariant is not a known force field (got \"" + ffVariant +
           "\")";
  }
  if (interFragVdwCutoff < 0.0) {
    return "zipper.interFragVdwCutoff must be >= 0";
  }
  if (fragShrugDisplacement < 0.0) {
    return "zipper.fragShrugDisplacement must be >= 0";
  }
  if (fragShrugForceConst < 0.0) {
    return "zipper.fragShrugForceConst must be >= 0";
  }
  return {};
}

const char *fragmentJoinerStatusMessage(FragmentJoinerStatus s) {
  switch (s) {
    case FragmentJoinerStatus::Ok:
      return "ok";
    case FragmentJoinerStatus::BadParams:
      return "the zipper parameters are not viable";
    case FragmentJoinerStatus::NoFragments:
      return "no fragments to join";
    case FragmentJoinerStatus::NoFragmentConformers:
      return "a fragment has no conformers";
    case FragmentJoinerStatus::NoScorer:
      return "no force field could be built for the molecule";
  }
  return "unknown";
}

std::vector<unsigned int> FragmentZipperContext::getIntraRotorBonds(
    const RigidRotorSearchParams &sp) const {
  std::vector<unsigned int> out;
  for (const auto &r : intraRotorBonds) {
    switch (r.type) {
      case IntraRotorType::Free:
      case IntraRotorType::PlanarAmide:
        if (sp.driveIntraFragmentTorsions) {
          out.push_back(r.bond);
        }
        break;
      case IntraRotorType::Atropisomer:
        if (sp.atropisomerSampling == AtropisomerSampling::Basin) {
          out.push_back(r.bond);
        }
        break;
    }
  }
  return out;
}

std::optional<IntraRotorType> FragmentZipperContext::getIntraType(
    unsigned int bond) const {
  for (const auto &r : intraRotorBonds) {
    if (r.bond == bond) {
      return r.type;
    }
  }
  return std::nullopt;
}

FragmentZipperContext zipFragments(const ROMol &mol, std::vector<int> atomFragment,
                                    std::vector<ZipFragment> fragments,
                                    std::vector<ZipJunction> junctions,
                                    FragmentZipperParams params,
                                    std::vector<IntraRotor> intraRotorBonds) {
  FragmentZipperContext ctx;
  ctx.mol = ROMol(mol);  // ROMol copy-ASSIGN is deleted; copy then move
  ctx.frags = std::move(fragments);
  ctx.intraRotorBonds = std::move(intraRotorBonds);
  ctx.ffVariant = params.ffVariant;
  if (!params.isValid()) {
    ctx.status = FragmentJoinerStatus::BadParams;
    return ctx;
  }
  if (ctx.mol.getNumConformers() == 0) {
    ctx.mol.addConformer(new Conformer(ctx.mol.getNumAtoms()), true);
  }
  if (ctx.frags.empty()) {
    ctx.status = FragmentJoinerStatus::NoFragments;
    return ctx;
  }

  // always root the fragment tree at the largest fragment
  //  This can save a bit of time during rotor driving
  for (unsigned int f = 0; f < ctx.frags.size(); ++f) {
    if (ctx.frags[f].atoms.size() > ctx.frags[ctx.root].atoms.size()) {
      ctx.root = f;
    }
  }

  // orient junctions into parent->child edges by BFS from the root
  std::vector<std::vector<size_t>> adj(ctx.frags.size());
  for (size_t j = 0; j < junctions.size(); ++j) {
    adj[junctions[j].fragA].push_back(j);
    adj[junctions[j].fragB].push_back(j);
  }
  std::vector<char> visited(ctx.frags.size(), 0);
  std::queue<unsigned int> q;
  q.push(ctx.root);
  visited[ctx.root] = 1;
  std::vector<std::pair<unsigned int, unsigned int>> junctionBonds;
  while (!q.empty()) {
    unsigned int f = q.front();
    q.pop();
    for (size_t ji : adj[f]) {
      const auto &J = junctions[ji];
      unsigned int other = (J.fragA == f) ? J.fragB : J.fragA;
      if (visited[other]) {
        continue;
      }
      visited[other] = 1;
      q.push(other);
      ZipEdge e;
      e.parentFrag = f;
      e.childFrag = other;
      if (J.fragA == f) {
        e.parentAtom = J.atomA;
        e.childAtom = J.atomB;
      } else {
        e.parentAtom = J.atomB;
        e.childAtom = J.atomA;
      }
      e.bondLen = J.bondLen;
      const Bond *b = ctx.mol.getBondBetweenAtoms(e.parentAtom, e.childAtom);
      e.bondIdx = b ? b->getIdx() : 0;
      ctx.edges.push_back(e);
      ctx.rotorBonds.push_back(e.bondIdx);
      junctionBonds.emplace_back(e.parentAtom, e.childAtom);
    }
  }

  // Check for embedding fails
  for (const auto &f : ctx.frags) {
    if (f.confs.empty()) {
      BOOST_LOG(rdWarningLog) << "FragmentZipper:: Fragment embedding failure" << std::endl;
      ctx.status = FragmentJoinerStatus::NoFragmentConformers;
      return ctx;
    }
  }
  
  {
    std::vector<unsigned int> lowest(ctx.frags.size(), 0);
    std::vector<double> flat =
        placeFragments(ctx.mol, ctx.frags, ctx.edges, ctx.root, lowest);
    auto *conf = new Conformer(ctx.mol.getNumAtoms());
    for (unsigned int a = 0; a < ctx.mol.getNumAtoms(); ++a) {
      conf->setAtomPos(
          a, RDGeom::Point3D(flat[3 * a], flat[3 * a + 1], flat[3 * a + 2]));
    }
    ctx.mol.clearConformers();
    ctx.mol.addConformer(conf, true);
  }

  // Scorer: default is the decomposed inter-fragment vdW (fast; assumes rigid
  // fragments).
  //  We can also score via a full FF, but this is slower.
  // NOTE: a fallback UFF embedder and scorer would widen the element coverage
  //  (MMFF cannot type boron, for one).  Slow, but possibly fine for Thompson
  //  sampling.  Until then an untypeable molecule reports NoScorer.
  if (params.useFullFFScorer) {
    ctx.scorer =
        makeFullFFScoreFn(ctx.mol, /*electrostatics=*/false, params.ffVariant,
                        /*nonBondedThresh=*/100.0);
  }
  if (!ctx.scorer) {
    ctx.scorer = makeInterFragmentScoreFn(
        ctx.mol, atomFragment, junctionBonds, false, params.ffVariant,
        params.interFragVdwCutoff, &ctx.scoreHandles);
  }
  if (!ctx.scorer) {
    ctx.status = FragmentJoinerStatus::NoScorer;
  }
  return ctx;
}

FragmentZipperContext zipFragments(FragmentZipperInput in,
                                    FragmentZipperParams params) {
  return zipFragments(*in.mol, std::move(in.atomFragment), std::move(in.fragments),
                       std::move(in.junctions), std::move(params),
                       std::move(in.intraRotorBonds));
}

namespace {
std::vector<double> placeFragments(
    const ROMol &mol, const std::vector<ZipFragment> &frags,
    const std::vector<ZipEdge> &edges, unsigned int root,
    const std::vector<unsigned int> &confChoice) {
  const unsigned int n = mol.getNumAtoms();
  std::vector<RDGeom::Point3D> buffer(n);
  std::vector<RDGeom::Point3D> parentExit(edges.size());

  // PRECONDITION: all frags have confs.
  const auto &rootPos = frags[root].confs[confChoice[root]].pos;
  for (unsigned int a : frags[root].atoms) {
    buffer[a] = rootPos[a];
  }
  for (size_t e = 0; e < edges.size(); ++e) {
    if (edges[e].parentFrag == root) {
      parentExit[e] = rootPos[edges[e].childAtom];
    }
  }

  // walk the BFS edges: each parent is already placed before its child edge is
  // hit
  for (size_t e = 0; e < edges.size(); ++e) {
    const ZipEdge &ed = edges[e];
    std::vector<RDGeom::Point3D> cc =
        frags[ed.childFrag].confs[confChoice[ed.childFrag]].pos;
    placeChildCoords(cc, ed.childAtom, ed.parentAtom, buffer[ed.parentAtom],
                     parentExit[e], ed.bondLen);
    for (unsigned int a : frags[ed.childFrag].atoms) {
      buffer[a] = cc[a];
    }
    // this child is a parent for its own downstream junctions
    for (size_t e2 = 0; e2 < edges.size(); ++e2) {
      if (edges[e2].parentFrag == ed.childFrag) {
        parentExit[e2] = cc[edges[e2].childAtom];
      }
    }
  }

  std::vector<double> flat(3 * static_cast<size_t>(n));
  for (unsigned int a = 0; a < n; ++a) {
    flat[3 * a] = buffer[a].x;
    flat[3 * a + 1] = buffer[a].y;
    flat[3 * a + 2] = buffer[a].z;
  }
  return flat;
}
}  // namespace

std::vector<double> FragmentZipperContext::placeAll(
    const std::vector<unsigned int> &confChoice) const {
  return placeFragments(mol, frags, edges, root, confChoice);
}

void FragmentZipperContext::symmetryDedup(
    std::vector<SearchResult> &v, double thr) const {
  if (thr <= 0.0 || v.size() <= 1) {
    return;
  }
  RMSDPruner acc(mol, thr);  // automorphisms computed once here
  std::vector<SearchResult> kept;
  kept.reserve(v.size());
  for (auto &r : v) {
    if (acc.add(r.coords)) {
      kept.push_back(std::move(r));
    }
  }
  v.swap(kept);
}


}  // namespace RDKit
