//
//  Symmetry-aware conformer dedup (see SymmetricRmsd.h).
//  Ported by Claude from: https://github.com/pandegroup/IRMSD
//   note: not SSE enabled for portability reasons, but hopefully
//         the compiler can take care of that.
//
#include "Utils/SymmetricRmsd.h"
#include "Utils/TheobaldRmsd.h"

#include <GraphMol/ROMol.h>
#include <GraphMol/new_canon.h>

#include <functional>
#include <limits>
#include <queue>

namespace RDKit {

std::vector<std::vector<unsigned int>> heavyAtomAutomorphisms(
    const ROMol &mol, std::vector<unsigned int> &heavyOut, size_t maxAutos) {
  heavyOut.clear();
  std::vector<int> molToHeavy(mol.getNumAtoms(), -1);
  for (const auto a : mol.atoms()) {
    if (a->getAtomicNum() > 1) {
      molToHeavy[a->getIdx()] = static_cast<int>(heavyOut.size());
      heavyOut.push_back(a->getIdx());
    }
  }
  const unsigned int H = static_cast<unsigned int>(heavyOut.size());

  std::vector<std::vector<unsigned int>> autos;
  autos.emplace_back(H);  // identity, filled below
  for (unsigned int p = 0; p < H; ++p) autos[0][p] = p;
  if (H == 0) return autos;

  // Canonical ranks WITHOUT tie-breaking: equal ranks == same symmetry class.
  std::vector<unsigned int> ranks;
  Canon::rankMolAtoms(mol, ranks, /*breakTies=*/false);

  std::vector<int> rk(H);
  std::vector<std::vector<int>> adj(H);
  std::vector<std::vector<char>> isAdj(H, std::vector<char>(H, 0));
  for (unsigned int p = 0; p < H; ++p)
    rk[p] = static_cast<int>(ranks[heavyOut[p]]);
  for (unsigned int p = 0; p < H; ++p) {
    const Atom *a = mol.getAtomWithIdx(heavyOut[p]);
    for (const auto nb : mol.atomNeighbors(a)) {
      const int hq = molToHeavy[nb->getIdx()];
      if (hq >= 0) {
        adj[p].push_back(hq);
        isAdj[p][hq] = 1;
      }
    }
  }

  // BFS visitation order so every atom after a component's root has an
  // already-placed neighbour -> the bond check prunes candidates hard.
  std::vector<int> order;
  order.reserve(H);
  std::vector<char> seen(H, 0);
  for (unsigned int s = 0; s < H; ++s) {
    if (seen[s]) continue;
    std::queue<int> q;
    q.push(static_cast<int>(s));
    seen[s] = 1;
    while (!q.empty()) {
      const int u = q.front();
      q.pop();
      order.push_back(u);
      for (int v : adj[u]) {
        if (!seen[v]) {
          seen[v] = 1;
          q.push(v);
        }
      }
    }
  }

  std::vector<int> perm(H, -1);
  std::vector<char> used(H, 0);
  bool overflow = false;
  std::function<void(unsigned int)> rec = [&](unsigned int pos) {
    if (overflow) return;
    if (pos == H) {
      std::vector<unsigned int> pm(H);
      for (unsigned int p = 0; p < H; ++p)
        pm[p] = static_cast<unsigned int>(perm[p]);
      // skip the identity (already stored first)
      bool ident = true;
      for (unsigned int p = 0; p < H && ident; ++p) ident = (pm[p] == p);
      if (!ident) {
        autos.push_back(std::move(pm));
        if (autos.size() >= maxAutos) overflow = true;
      }
      return;
    }
    const int a = order[pos];
    for (unsigned int c = 0; c < H; ++c) {
      if (used[c] || rk[c] != rk[a] || adj[c].size() != adj[a].size()) continue;
      bool ok = true;
      for (int nb : adj[a]) {
        if (perm[nb] >= 0 && !isAdj[c][perm[nb]]) {  // placed neighbour of a
          ok = false;
          break;
        }
      }
      if (!ok) continue;
      perm[a] = static_cast<int>(c);
      used[c] = 1;
      rec(pos + 1);
      used[c] = 0;
      perm[a] = -1;
      if (overflow) return;
    }
  };
  rec(0);

  if (overflow) {  // pathological symmetry: fall back to identity only
    autos.resize(1);
  }
  return autos;
}

RMSDPruner::RMSDPruner(const ROMol &mol, double rmsThresh)
    : d_thresh(rmsThresh) {
  auto perms = heavyAtomAutomorphisms(mol, d_heavy);
  const unsigned int H = static_cast<unsigned int>(d_heavy.size());
  d_permB.reserve(perms.size());
  for (const auto &perm : perms) {
    std::vector<unsigned int> idxB(H);
    for (unsigned int p = 0; p < H; ++p) idxB[p] = d_heavy[perm[p]];
    d_permB.push_back(std::move(idxB));
  }
  // size the reusable QCP scratch once (see header): kept-centred [3H],
  // candidate- centred per automorphism [P*3H], and the per-automorphism inner
  // products [P].
  const size_t P = d_permB.size();
  d_scratchA.assign(static_cast<size_t>(3) * H, 0.0);
  d_candB.assign(P * static_cast<size_t>(3) * H, 0.0);
  d_candGb.assign(P, 0.0);
}

double RMSDPruner::RMSD(
    const std::vector<double> &coords) const {
  const unsigned int H = static_cast<unsigned int>(d_heavy.size());
  if (d_kept.empty() || H == 0) {
    return std::numeric_limits<double>::infinity();
  }
  const unsigned int P = static_cast<unsigned int>(d_permB.size());
  const double invH = 1.0 / static_cast<double>(H);

  // Centre the candidate's heavy atoms under EACH automorphism once
  // (independent of which kept conformer we compare against) into the reused
  // scratch d_candB/d_candGb.
  for (unsigned int p = 0; p < P; ++p) {
    const auto &idxB = d_permB[p];
    double cx = 0, cy = 0, cz = 0;
    for (unsigned int k = 0; k < H; ++k) {
      const unsigned int ib = 3 * idxB[k];
      cx += coords[ib];
      cy += coords[ib + 1];
      cz += coords[ib + 2];
    }
    cx *= invH;
    cy *= invH;
    cz *= invH;
    double *B = d_candB.data() + static_cast<size_t>(p) * 3 * H;
    double Gb = 0.0;
    for (unsigned int k = 0; k < H; ++k) {
      const unsigned int ib = 3 * idxB[k];
      const double bx = coords[ib] - cx, by = coords[ib + 1] - cy,
                   bz = coords[ib + 2] - cz;
      B[3 * k] = bx;
      B[3 * k + 1] = by;
      B[3 * k + 2] = bz;
      Gb += bx * bx + by * by + bz * bz;
    }
    d_candGb[p] = Gb;
  }

  // Compare in MSD space (theobaldMSD returns mean-square deviation): RMSD <
  // thresh
  // <=> MSD < thresh^2, so the early-outs are identical to the old sqrt-space
  // breaks.
  const double threshMSD = d_thresh * d_thresh;
  double bestMSD = std::numeric_limits<double>::infinity();
  double *A = d_scratchA.data();
  for (const auto &kept : d_kept) {
    // centre this kept conformer's heavy atoms once (shared across all
    // automorphisms)
    double cx = 0, cy = 0, cz = 0;
    for (unsigned int k = 0; k < H; ++k) {
      const unsigned int ia = 3 * d_heavy[k];
      cx += kept[ia];
      cy += kept[ia + 1];
      cz += kept[ia + 2];
    }
    cx *= invH;
    cy *= invH;
    cz *= invH;
    double Ga = 0.0;
    for (unsigned int k = 0; k < H; ++k) {
      const unsigned int ia = 3 * d_heavy[k];
      const double ax = kept[ia] - cx, ay = kept[ia + 1] - cy,
                   az = kept[ia + 2] - cz;
      A[3 * k] = ax;
      A[3 * k + 1] = ay;
      A[3 * k + 2] = az;
      Ga += ax * ax + ay * ay + az * az;
    }
    double m = std::numeric_limits<double>::infinity();
    for (unsigned int p = 0; p < P; ++p) {
      const double msd =
          theobaldMSD(H, A, d_candB.data() + static_cast<size_t>(p) * 3 * H, Ga,
                      d_candGb[p]);
      if (msd < m) m = msd;
      if (m < threshMSD) break;  // already a duplicate of this kept conformer
    }
    if (m < bestMSD) bestMSD = m;
    if (bestMSD < threshMSD) break;  // duplicate of some kept conformer
  }
  return std::sqrt(bestMSD);
}

bool RMSDPruner::add(const std::vector<double> &coords) {
  if (d_thresh > 0.0 && RMSD(coords) < d_thresh) {
    return false;
  }
  d_kept.push_back(coords);
  return true;
}

}  // namespace RDKit
