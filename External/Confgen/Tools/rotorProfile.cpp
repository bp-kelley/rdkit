//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  DIAGNOSTIC.  For the widest-moving rotors of each molecule:
//    * scan the junction dihedral finely and report the energy landscape
//    * census clashes by TOPOLOGICAL distance from the rotor bond -- local
//      junction strain vs a distal collision
//    * ask whether a COARSE grid predicts the fine minimum, i.e. whether a
//      cheap coarse pass could seed a prior for a sampler
//
//  usage: rotorProfile <sdf> [molName|ALL] [topRotors] [nScramble]
//
//  nScramble > 0 re-runs each rotor's profile with every OTHER rotor driven to
//  a random angle first.  The native (crystal) pose flatters a single-rotor
//  landscape; if the coarse-grid prior only works from there it is an artifact.
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <Confgen/FragmentConfGen.h>
#include <Confgen/Search/RotorDriver.h>
#include <Confgen/Search/InterFragScore.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <numeric>
#include <memory>
#include <random>
#include <string>
#include <vector>
using namespace RDKit;

namespace {

constexpr int kFineStep = 10;                 // degrees
constexpr int kNFine = 360 / kFineStep;
//! An angle is "forbidden" if it costs this much over the rotor's own best.
constexpr double kForbidden = 50.0;

struct Agg {
  const char *label = "";
  std::map<int, long> clashByDist;   // topological distance -> clashes
  long rotors = 0, mols = 0;
  long forbidden = 0, sampled = 0;   // forbidden angle fraction
  // coarse-grid predictive power
  std::map<int, long> coarseHit;     // step -> times coarse min is in the fine basin
  std::map<int, long> coarseTries;
  std::map<int, double> coarseRegret;  // step -> summed (E(coarseArgmin)-Emin)
};

double angDiff(double a, double b) {
  double d = std::fabs(a - b);
  return std::min(d, 360.0 - d);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    printf("usage: rotorProfile <sdf> [molName|ALL] [topRotors]\n");
    return 1;
  }
  const std::string want = argc > 2 ? argv[2] : "ALL";
  const int nTop = argc > 3 ? std::atoi(argv[3]) : 3;
  const int nScram = argc > 4 ? std::atoi(argv[4]) : 0;
  const bool all = (want == "ALL");
  // CHAINS mode: per-molecule rotor-topology descriptors, TSV to stdout.
  // Rotor COUNT alone conflates very different molecules: ten rotors in one
  // long tail is far floppier than ten split into five short arms between
  // rings, because a long chain's rotors move each other's downstream atoms.
  // A "chain" here = a maximal run of rotatable bonds joined through atoms
  // that are NOT in a ring; a ring system therefore terminates a chain.
  const bool chainsMode = (want == "CHAINS");
  if (chainsMode) {
    printf("name\tnRotors\tnChains\tmaxChain\tmeanChain\tmedianChain\n");
  }
  Agg agg, scram;
  agg.label = "NATIVE pose";
  scram.label = "SCRAMBLED pose";
  std::mt19937 rng(0xC0FFEE);

  SDMolSupplier supp(argv[1], true, false, false);
  while (!supp.atEnd()) {
    std::unique_ptr<ROMol> raw(supp.next());
    if (!raw) continue;
    std::string nm;
    raw->getPropIfPresent("_Name", nm);
    if (!all && !chainsMode && nm != want) continue;

    std::unique_ptr<RWMol> m(new RWMol(*raw));
    MolOps::addHs(*m, false, true);
    if (!m->getRingInfo()->isInitialized()) MolOps::fastFindRings(*m);
    const auto rb = FragmentConfGen::findRotatableBonds(*m, false, false);
    std::vector<unsigned int> bonds(rb.inter.begin(), rb.inter.end());
    if (bonds.empty()) continue;
    auto score = makeFullFFScoreFn(*m, false, "MMFF94", 100.0);
    if (!score) continue;
    RotorDriver drv(*m, bonds, -1, score);
    const unsigned int n = m->getNumAtoms();
    const double *dm = MolOps::getDistanceMat(*m);
    auto *pt = PeriodicTable::getTable();
    const std::vector<double> home = drv.positions();
    ++agg.mols;

    std::vector<unsigned int> order;
    for (unsigned int r = 0; r < drv.numRotors(); ++r) order.push_back(r);
    std::sort(order.begin(), order.end(), [&](unsigned int a, unsigned int b) {
      return drv.movingAtoms(a).size() > drv.movingAtoms(b).size();
    });

    if (chainsMode) {
      // union-find over rotatable bonds; two are in the same chain when they
      // share an atom that is not in a ring
      const size_t nb = bonds.size();
      std::vector<int> par(nb);
      for (size_t i = 0; i < nb; ++i) par[i] = static_cast<int>(i);
      std::function<int(int)> find = [&](int x) {
        while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; }
        return x;
      };
      for (size_t i = 0; i < nb; ++i) {
        const Bond *bi = m->getBondWithIdx(bonds[i]);
        for (size_t j = i + 1; j < nb; ++j) {
          const Bond *bj = m->getBondWithIdx(bonds[j]);
          for (unsigned int a : {bi->getBeginAtomIdx(), bi->getEndAtomIdx()}) {
            if ((a == bj->getBeginAtomIdx() || a == bj->getEndAtomIdx()) &&
                !m->getRingInfo()->numAtomRings(a)) {
              par[find(static_cast<int>(i))] = find(static_cast<int>(j));
            }
          }
        }
      }
      std::map<int, int> sz;
      for (size_t i = 0; i < nb; ++i) sz[find(static_cast<int>(i))]++;
      std::vector<int> lens;
      for (const auto &kv : sz) lens.push_back(kv.second);
      std::sort(lens.begin(), lens.end());
      const double mean =
          lens.empty() ? 0.0
                       : std::accumulate(lens.begin(), lens.end(), 0.0) / lens.size();
      const double med = lens.empty() ? 0.0 : lens[lens.size() / 2];
      printf("%s\t%zu\t%zu\t%d\t%.2f\t%.1f\n", nm.c_str(), nb, lens.size(),
             lens.empty() ? 0 : lens.back(), mean, med);
      continue;
    }
    if (!all) printf("%s: %u atoms, %zu rotors\n\n", nm.c_str(), n, bonds.size());

    const int top = std::min<size_t>(nTop, order.size());
    for (int t = 0; t < top; ++t) {
      const unsigned int r = order[t];
      const auto tq = drv.torsion(r);
      const unsigned int j = tq[1], k = tq[2];
      ++agg.rotors;

      // measure this rotor's landscape from the CURRENT pose in `base`
      auto measure = [&](const std::vector<double> &base, Agg &A, bool census) {
      std::vector<double> E(kNFine);
      double lo = 1e30, hi = -1e30;
      int loI = 0, hiI = 0;
      for (int i = 0; i < kNFine; ++i) {
        drv.positions() = base;
        E[i] = drv.drive(r, i * kFineStep);
        if (E[i] < lo) { lo = E[i]; loI = i; }
        if (E[i] > hi) { hi = E[i]; hiI = i; }
      }
      for (int i = 0; i < kNFine; ++i) {
        ++A.sampled;
        if (E[i] > lo + kForbidden) ++A.forbidden;
      }

      // ---- does a COARSE grid land in the fine basin?
      for (int step : {30, 60, 120}) {
        double cBest = 1e30;
        int cArg = 0;
        for (int a = 0; a < 360; a += step) {
          const int i = (a / kFineStep) % kNFine;
          if (E[i] < cBest) { cBest = E[i]; cArg = a; }
        }
        ++A.coarseTries[step];
        // "same basin" = within one coarse step of the fine minimum
        if (angDiff(cArg, loI * kFineStep) <= step) ++A.coarseHit[step];
        A.coarseRegret[step] += cBest - lo;
      }
      if (!census) return;

      // ---- clash census at both ends
      for (int which = 0; which < 2; ++which) {
        const int useI = which ? hiI : loI;
        drv.positions() = base;
        drv.drive(r, useI * kFineStep);
        const auto &p = drv.positions();
        std::map<int, int> byDist;
        int total = 0;
        for (unsigned int x = 0; x < n; ++x) {
          for (unsigned int y = x + 1; y < n; ++y) {
            if (dm[x * n + y] < 4) continue;
            const double dx = p[3*x]-p[3*y], dy = p[3*x+1]-p[3*y+1],
                         dz = p[3*x+2]-p[3*y+2];
            const double d = std::sqrt(dx*dx + dy*dy + dz*dz);
            const double rs = pt->getRvdw(m->getAtomWithIdx(x)->getAtomicNum()) +
                              pt->getRvdw(m->getAtomWithIdx(y)->getAtomicNum());
            if (d >= 0.75 * rs) continue;
            const int dj = static_cast<int>(
                std::min(std::min(dm[x*n+j], dm[x*n+k]),
                         std::min(dm[y*n+j], dm[y*n+k])));
            byDist[dj]++;
            ++total;
            if (which) A.clashByDist[dj]++;   // aggregate the WORST angle
          }
        }
        if (!all) {
          printf("  rotor %u (bond %u-%u, moves %zu/%u): span %.1f kcal\n", r, j,
                 k, drv.movingAtoms(r).size(), n, hi - lo);
          printf("    %s %3d deg E=%.1f  %d clash(es)\n", which ? "WORST" : "best ",
                 useI * kFineStep, which ? hi : lo, total);
          for (const auto &kv : byDist)
            printf("      %2d bonds away : %3d\n", kv.first, kv.second);
        }
      }
      };  // measure

      measure(home, agg, /*census=*/true);
      // ---- the same rotor, but with every OTHER rotor moved off-native
      for (int sIt = 0; sIt < nScram; ++sIt) {
        std::vector<double> base = home;
        std::uniform_int_distribution<int> pick(0, kNFine - 1);
        drv.positions() = base;
        for (unsigned int q : drv.numRotorAtoms()) {   // coarse -> fine order
          if (q == r) continue;
          drv.setDihedral(q, pick(rng) * kFineStep);
        }
        base = drv.positions();
        measure(base, scram, /*census=*/true);
      }
    }
    if (!all && nm == want) break;
  }

  auto report = [&](const Agg &A, bool clash) {
    printf("\n================ %s ================\n", A.label);
    if (clash) {
    printf("CLASH LOCALITY at each rotor's worst angle\n");
    long tot = 0;
    for (const auto &kv : A.clashByDist) tot += kv.second;
    long cum = 0;
    for (const auto &kv : A.clashByDist) {
      cum += kv.second;
      printf("  %2d bonds from the rotor bond : %6ld  (%5.1f%%, cum %5.1f%%)\n",
             kv.first, kv.second, 100.0 * kv.second / tot, 100.0 * cum / tot);
    }
    }
    if (!A.sampled) return;
    printf("ANGLE LANDSCAPE\n");
    printf("  forbidden angles (> best + %.0f kcal): %.1f%% of %ld sampled\n",
           kForbidden, 100.0 * A.forbidden / A.sampled, A.sampled);
    printf("COARSE GRID as a PRIOR\n");
    for (int step : {30, 60, 120}) {
      const long tries = A.coarseTries.count(step) ? A.coarseTries.at(step) : 0;
      if (!tries) continue;
      printf("  %3d deg grid: basin hit %5.1f%%   mean regret %7.2f kcal\n", step,
             100.0 * A.coarseHit.at(step) / tries, A.coarseRegret.at(step) / tries);
    }
  };
  if (all && !chainsMode) {
    printf("=== %ld molecules, %ld rotors (top-%d widest each) ===\n", agg.mols,
           agg.rotors, nTop);
    report(agg, true);
    if (nScram) report(scram, true);
  }
  return 0;
}
