//  Copyright (C) 2026 Pande Group, Glysade Inc
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
//  Portable scalar port of IRMSD's msd_atom_major + msdFromMandG (Theobald
//  QCP).
//
#include "Utils/TheobaldRmsd.h"

#include <cmath>

namespace RDKit {
namespace {

// Newton-Raphson on the QCP characteristic quartic P(l) = l^4 + C2 l^2 + C1 l +
// C0 (C4 = 1, C3 = 0), started from the upper bound l0 = (Ga + Gb)/2 -- from
// there it converges monotonically to the largest root, which is what QCP
// needs.  This is the scheme in IRMSD's NewtonSolve, written in double
// precision.
double newtonSolve(double lambda, double C0, double C1, double C2) {
  double lambda_old = lambda;
  for (int i = 0; i < 50; ++i) {
    lambda_old = lambda;
    const double l2 = lambda_old * lambda_old;
    const double bb = (l2 + C2) * lambda_old;  // l^3 + C2 l
    const double aa = bb + C1;                 // l^3 + C2 l + C1
    const double denom =
        2.0 * l2 * lambda_old + bb + aa;  // P'(l) = 4l^3 + 2 C2 l + C1
    if (denom == 0.0) {
      break;
    }
    lambda = lambda_old - (aa * lambda_old + C0) / denom;  // l - P(l)/P'(l)
    if (std::fabs(lambda - lambda_old) < std::fabs(1e-11 * lambda)) {
      break;
    }
  }
  return lambda;
}

// Correlation matrix, column-major: M[i + 3*j] = sum_k a_k[i] * b_k[j]  (the
// layout msdFromMandG indexes as M[row + 3*col]).  This tight reduction is what
// the IRMSD SSE path (msd_atom_major, deinterleaved-load + haddps) accelerates;
// here it is a plain loop that clang -O3 auto-vectorises to SSE/AVX on x86 and
// NEON on arm64.  A hand- tuned x86 intrinsic version is the single drop-in
// point: guard it with
// `#if defined(__SSE3__)` and produce this same M (fragments are small +
// single- threaded, so the scalar form is the tested default).
void correlationMatrix(unsigned int nAtoms, const double *a, const double *b,
                       double M[9]) {
  for (int i = 0; i < 9; ++i) {
    M[i] = 0.0;
  }
  for (unsigned int k = 0; k < nAtoms; ++k) {
    const double ax = a[3 * k], ay = a[3 * k + 1], az = a[3 * k + 2];
    const double bx = b[3 * k], by = b[3 * k + 1], bz = b[3 * k + 2];
    M[0] += ax * bx;
    M[1] += ay * bx;
    M[2] += az * bx;  // column bx
    M[3] += ax * by;
    M[4] += ay * by;
    M[5] += az * by;  // column by
    M[6] += ax * bz;
    M[7] += ay * bz;
    M[8] += az * bz;  // column bz
  }
}

}  // namespace

double theobaldMSD(unsigned int nAtoms, const double *a, const double *b,
                   double Ga, double Gb) {
  if (nAtoms == 0) {
    return 0.0;
  }
  double M[9];
  correlationMatrix(nAtoms, a, b, M);

  const int m = 3;
  const double k00 = M[0 + 0 * m] + M[1 + 1 * m] + M[2 + 2 * m];
  const double k01 = M[1 + 2 * m] - M[2 + 1 * m];
  const double k02 = M[2 + 0 * m] - M[0 + 2 * m];
  const double k03 = M[0 + 1 * m] - M[1 + 0 * m];
  const double k11 = M[0 + 0 * m] - M[1 + 1 * m] - M[2 + 2 * m];
  const double k12 = M[0 + 1 * m] + M[1 + 0 * m];
  const double k13 = M[2 + 0 * m] + M[0 + 2 * m];
  const double k22 = -M[0 + 0 * m] + M[1 + 1 * m] - M[2 + 2 * m];
  const double k23 = M[1 + 2 * m] + M[2 + 1 * m];
  const double k33 = -M[0 + 0 * m] - M[1 + 1 * m] + M[2 + 2 * m];

  double C2 = 0.0;
  for (int i = 0; i < 9; ++i) {
    C2 += M[i] * M[i];
  }
  C2 *= -2.0;

  const double detM = M[0] * (M[4] * M[8] - M[5] * M[7]) +
                      M[3] * (M[7] * M[2] - M[8] * M[1]) +
                      M[6] * (M[1] * M[5] - M[2] * M[4]);
  const double detK = k01 * k01 * k23 * k23 - k22 * k33 * k01 * k01 +
                      2 * k33 * k01 * k02 * k12 - 2 * k01 * k02 * k13 * k23 -
                      2 * k01 * k03 * k12 * k23 + 2 * k22 * k01 * k03 * k13 +
                      k02 * k02 * k13 * k13 - k11 * k33 * k02 * k02 -
                      2 * k02 * k03 * k12 * k13 + 2 * k11 * k02 * k03 * k23 +
                      k03 * k03 * k12 * k12 - k11 * k22 * k03 * k03 -
                      k00 * k33 * k12 * k12 + 2 * k00 * k12 * k13 * k23 -
                      k00 * k22 * k13 * k13 - k00 * k11 * k23 * k23 +
                      k00 * k11 * k22 * k33;

  const double C1 = -8.0 * detM;
  const double C0 = detK;

  const double lambda = newtonSolve((Ga + Gb) / 2.0, C0, C1, C2);
  const double msd = (Ga + Gb - 2.0 * lambda) / static_cast<double>(nAtoms);
  return msd > 0.0 ? msd : 0.0;
}

double qcpRmsd(const std::vector<RDGeom::Point3D> &a,
               const std::vector<RDGeom::Point3D> &b,
               const std::vector<unsigned int> &idx) {
  const unsigned int n = static_cast<unsigned int>(idx.size());
  if (n == 0) {
    return 0.0;
  }
  RDGeom::Point3D ca, cb;  // default-constructed to the origin
  for (auto i : idx) {
    ca += a[i];
    cb += b[i];
  }
  ca /= n;
  cb /= n;

  std::vector<double> A(3 * n), B(3 * n);
  double Ga = 0.0, Gb = 0.0;
  for (unsigned int k = 0; k < n; ++k) {
    const RDGeom::Point3D pa = a[idx[k]] - ca;
    const RDGeom::Point3D pb = b[idx[k]] - cb;
    A[3 * k] = pa.x;
    A[3 * k + 1] = pa.y;
    A[3 * k + 2] = pa.z;
    B[3 * k] = pb.x;
    B[3 * k + 1] = pb.y;
    B[3 * k + 2] = pb.z;
    Ga += pa.x * pa.x + pa.y * pa.y + pa.z * pa.z;
    Gb += pb.x * pb.x + pb.y * pb.y + pb.z * pb.z;
  }
  return std::sqrt(theobaldMSD(n, A.data(), B.data(), Ga, Gb));
}

double qcpRmsd(const std::vector<double> &a, const std::vector<double> &b,
               const std::vector<unsigned int> &idxA,
               const std::vector<unsigned int> &idxB) {
  const unsigned int n = static_cast<unsigned int>(idxA.size());
  if (n == 0 || idxB.size() != n) {
    return 0.0;
  }
  double cax = 0, cay = 0, caz = 0, cbx = 0, cby = 0, cbz = 0;
  for (unsigned int k = 0; k < n; ++k) {
    const unsigned int ia = 3 * idxA[k], ib = 3 * idxB[k];
    cax += a[ia];
    cay += a[ia + 1];
    caz += a[ia + 2];
    cbx += b[ib];
    cby += b[ib + 1];
    cbz += b[ib + 2];
  }
  const double inv = 1.0 / n;
  cax *= inv;
  cay *= inv;
  caz *= inv;
  cbx *= inv;
  cby *= inv;
  cbz *= inv;

  std::vector<double> A(3 * n), B(3 * n);
  double Ga = 0.0, Gb = 0.0;
  for (unsigned int k = 0; k < n; ++k) {
    const unsigned int ia = 3 * idxA[k], ib = 3 * idxB[k];
    const double ax = a[ia] - cax, ay = a[ia + 1] - cay, az = a[ia + 2] - caz;
    const double bx = b[ib] - cbx, by = b[ib + 1] - cby, bz = b[ib + 2] - cbz;
    A[3 * k] = ax;
    A[3 * k + 1] = ay;
    A[3 * k + 2] = az;
    B[3 * k] = bx;
    B[3 * k + 1] = by;
    B[3 * k + 2] = bz;
    Ga += ax * ax + ay * ay + az * az;
    Gb += bx * bx + by * by + bz * bz;
  }
  return std::sqrt(theobaldMSD(n, A.data(), B.data(), Ga, Gb));
}

}  // namespace RDKit
