#!/usr/bin/env python
#  Copyright (C) 2026 Greg Landrum and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Build the data tables used by the NormalizedDescriptors extension.

descriptastorus has two separate ways of normalizing a descriptor value v,
and this tool keeps them separate:

1. Fitted distribution (RDKit2DNormalized):

    normalized(v) = clip(dist.cdf(clip(v, min, max), *shape, loc, scale), 0, 1)

   where `dist` is the scipy.stats distribution that best fits the
   descriptor's values over a reference set of molecules. The C++ code does
   not link scipy, so every fitted CDF is evaluated on an adaptive grid and
   written to a table the C++ code linearly interpolates (within 1e-5 of
   scipy). Nothing about the data itself goes into this table, only the
   fitted function.

2. Histogram (RDKit2DHistogramNormalized): no fitting at all. The raw
   descriptor values are binned with numpy.histogram (min(1000, number of
   distinct values) equal-width bins) exactly as descriptastorus'
   data/d_descriptors/make_histdists.py does, giving (left edge, cumulative
   fraction) pairs, and v is looked up the way descriptastorus does it:

    p = bisect.bisect(bins, (v,))       # number of edges < v
    normalized(v) = bins[p][1] if p < len(bins) else 1.0

   No interpolation.

Files (paths relative to External/NormalizedDescriptors):

  data/normalized_descriptor_fits.json        method 1: fitted
                                              distributions (source of truth)
  data/normalized_descriptor_cdfs.txt         method 1: the interpolation
                                              table the C++ code reads
                                              (generated, don't edit)
  data/normalized_descriptor_histograms.txt   method 2: histogram CDF bins

Subcommands:

  import-dists  seed or update the fits file from a descriptastorus dists.py
  compute       compute descriptor values for a SMILES file (by default
                descriptastorus' data/chembl_100k.smi) and write one
                d_<name>.gz sample file per descriptor, in the same format as
                descriptastorus' data/d_descriptors
  fit           fit a distribution to each sample file and store the best
                fit in the fits file
  missing       list RDKit descriptors that have no fit (or histogram) yet
  table         evaluate every fit on a dense grid and write the table
  check         compare the table's interpolation against scipy
  import-hists  seed or update the histogram table from a descriptastorus
                hists.py
  histogram     bin each sample file and add it to the histogram table
  check-hists   compare histogram table lookups against a hists.py

Regenerating the table for the current fits:

  python fit_normalized_descriptors.py table

Adding new descriptors (e.g. everything RDKit has that isn't fit yet):

  python fit_normalized_descriptors.py compute \\
      --smiles /path/to/descriptastorus/data/chembl_100k.smi \\
      --out samples $(python fit_normalized_descriptors.py missing)
  python fit_normalized_descriptors.py fit samples/d_*.gz
  python fit_normalized_descriptors.py table
  python fit_normalized_descriptors.py check --samples samples
  python fit_normalized_descriptors.py histogram samples/d_*.gz

Requires numpy and scipy; `compute` also needs the RDKit python wrappers.
"""

import argparse
import ast
import bisect
import gzip
import json
import math
import multiprocessing
import os
import sys
import time

import numpy as np
import scipy
import scipy.stats as st

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(HERE, "..", "data"))
DEFAULT_FITS = os.path.join(DATA_DIR, "normalized_descriptor_fits.json")
DEFAULT_TABLE = os.path.join(DATA_DIR, "normalized_descriptor_cdfs.txt")
DEFAULT_HISTS = os.path.join(DATA_DIR, "normalized_descriptor_histograms.txt")

TABLE_FORMAT_VERSION = 1
FITS_FORMAT_VERSION = 1

# scipy families descriptastorus' dists.py uses; the default candidates for
# new fits so new descriptors are drawn from the same pool as the old ones
DEFAULT_FAMILIES = [
  "alpha", "beta", "betaprime", "burr", "cauchy", "chi", "dgamma", "dweibull",
  "exponnorm", "exponpow", "exponweib", "fisk", "foldcauchy", "foldnorm",
  "genexpon", "genhalflogistic", "genlogistic", "gennorm", "genpareto",
  "gengamma", "gibrat", "gompertz", "halfgennorm", "halflogistic", "halfnorm",
  "hypsecant", "invgamma", "invweibull", "johnsonsb", "johnsonsu", "laplace",
  "logistic", "loggamma", "lognorm", "lomax", "mielke", "nct", "ncx2",
  "pareto", "pearson3", "powerlaw", "recipinvgauss", "t", "tukeylambda",
  "wald",
  # slow to fit, use --families to opt in: gausshyper, ncf
]

# descriptastorus make_histdists.py uses at most this many bins
MAX_HIST_BINS = 1000


def log(*args):
  print(*args, file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# distributions


def scipy_dist(name):
  """Look up a scipy.stats distribution, allowing for renames."""
  if name in ("gilbrat", "gibrat"):
    name = "gibrat" if hasattr(st, "gibrat") else "gilbrat"
  return getattr(st, name)


def make_cdf(fit):
  """The descriptastorus normalization for one fit, vectorized."""
  lo, hi = fit["min"], fit["max"]
  dist = scipy_dist(fit["dist"])
  params = fit["params"]
  args, loc, scale = params[:-2], params[-2], params[-1]

  def cdf(x):
    x = np.clip(np.asarray(x, dtype=float), lo, hi)
    with np.errstate(all="ignore"):
      v = dist.cdf(x, *args, loc=loc, scale=scale)
    return np.clip(v, 0.0, 1.0)

  return cdf


# ---------------------------------------------------------------------------
# fits file


def load_fits(path):
  if not os.path.exists(path):
    return {"version": FITS_FORMAT_VERSION, "descriptors": {}}
  with open(path) as f:
    fits = json.load(f)
  if fits.get("version") != FITS_FORMAT_VERSION:
    raise ValueError(f"{path}: unsupported fits version {fits.get('version')}")
  return fits


def save_fits(fits, path):
  fits["descriptors"] = dict(sorted(fits["descriptors"].items()))
  with open(path, "w") as f:
    json.dump(fits, f, indent=1)
    f.write("\n")
  log(f"Wrote {len(fits['descriptors'])} fits to {path}")


def read_dists_py(path):
  """Read a descriptastorus dists.py without importing descriptastorus."""
  ns = {"inf": math.inf}
  with open(path) as f:
    exec(f.read(), ns)
  return ns["dists"]


def cmd_import_dists(args):
  fits = load_fits(args.fits)
  dists = read_dists_py(args.dists_py)
  source = args.source or f"descriptastorus {os.path.basename(args.dists_py)}"
  n = 0
  for name, (dist, params, lo, hi, mean, std) in dists.items():
    if name in fits["descriptors"] and not args.replace:
      continue
    fits["descriptors"][name] = {
      "dist": dist,
      "params": [float(p) for p in params],
      "min": float(lo),
      "max": float(hi),
      "mean": float(mean),
      "std": float(std),
      "source": source,
    }
    n += 1
  log(f"Imported {n} of {len(dists)} fits from {args.dists_py}")
  save_fits(fits, args.fits)


# ---------------------------------------------------------------------------
# samples (descriptastorus data/d_descriptors format: gzipped python list)


def read_sample(path):
  with gzip.open(path, "rt") as f:
    txt = f.read()
  try:
    vals = ast.literal_eval(txt)
  except ValueError:
    # older files can contain bare inf / nan
    vals = eval(txt, {"__builtins__": {}}, {"inf": math.inf, "nan": math.nan})
  vals = np.asarray(vals, dtype=float)
  return vals[np.isfinite(vals)]


def sample_name(path):
  base = os.path.basename(path)
  if base.endswith(".gz"):
    base = base[:-3]
  if base.startswith("d_"):
    base = base[2:]
  return base


def write_sample(path, values):
  with gzip.open(path, "wt") as f:
    f.write(repr([float(v) for v in values]))


# ---------------------------------------------------------------------------
# compute


def rdkit_descriptor_functions():
  from rdkit.Chem import Descriptors
  return dict(Descriptors._descList)


_worker_names = None


def _compute_init(names):
  global _worker_names
  _worker_names = names


def _compute_chunk(smiles):
  from rdkit import Chem, RDLogger
  RDLogger.DisableLog("rdApp.*")
  funcs = rdkit_descriptor_functions()
  res = {n: [] for n in _worker_names}
  nfailed = 0
  for smi in smiles:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
      nfailed += 1
      continue
    for n in _worker_names:
      try:
        res[n].append(float(funcs[n](mol)))
      except Exception:
        pass
  return res, nfailed


def read_smiles(path):
  smiles = []
  with open(path) as f:
    for line in f:
      line = line.strip()
      if not line or line.startswith("#"):
        continue
      smi = line.split()[0]
      if smi.lower() == "smiles":
        continue
      smiles.append(smi)
  return smiles


def cmd_compute(args):
  funcs = rdkit_descriptor_functions()
  unknown = [n for n in args.descriptors if n not in funcs]
  if unknown:
    sys.exit(f"Unknown RDKit descriptors: {' '.join(unknown)}")
  smiles = read_smiles(args.smiles)
  log(f"Computing {len(args.descriptors)} descriptors for {len(smiles)} "
      f"molecules from {args.smiles}")
  chunks = [smiles[i:i + 500] for i in range(0, len(smiles), 500)]
  values = {n: [] for n in args.descriptors}
  nfailed = 0
  with multiprocessing.Pool(args.nprocs or None, _compute_init,
                            (args.descriptors, )) as pool:
    for i, (res, nf) in enumerate(pool.imap(_compute_chunk, chunks)):
      if i % 20 == 19:
        log(f"  {min((i + 1) * 500, len(smiles))} molecules done")
      nfailed += nf
      for n, vs in res.items():
        values[n].extend(vs)
  if nfailed:
    log(f"{nfailed} SMILES failed to parse and were skipped")
  os.makedirs(args.out, exist_ok=True)
  for n, vs in values.items():
    path = os.path.join(args.out, f"d_{n}.gz")
    write_sample(path, vs)
    log(f"Wrote {len(vs)} values to {path}")


def cmd_missing(args):
  if args.histograms:
    have = read_hist_table(args.table)
  else:
    have = load_fits(args.fits)["descriptors"]
  for n in rdkit_descriptor_functions():
    if n not in have:
      print(n)


# ---------------------------------------------------------------------------
# fit


def ks_distance(sorted_vals, cdf_vals):
  """max |cdf(v) - fraction of values <= v| over the observed values.

  This is the KS statistic without the jump term, so count descriptors are
  judged on the values they actually take (a continuous CDF can't match the
  jumps of a discrete one between integers, and nothing is normalized there).
  """
  n = len(sorted_vals)
  ecdf = np.searchsorted(sorted_vals, sorted_vals, side="right") / n
  return float(np.max(np.abs(ecdf - cdf_vals)))


def _fit_one(task):
  family, fit_vals, all_vals = task
  t0 = time.time()
  try:
    dist = scipy_dist(family)
    with np.errstate(all="ignore"):
      import warnings
      with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        params = dist.fit(fit_vals)
        cdf = dist.cdf(all_vals, *params)
    if not np.all(np.isfinite(params)) or not np.all(np.isfinite(cdf)):
      return family, None, math.inf, time.time() - t0
    return family, [float(p) for p in params], ks_distance(
      all_vals, np.clip(cdf, 0, 1)), time.time() - t0
  except Exception as e:  # many families fail on some data
    return family, None, math.inf, time.time() - t0


def fit_sample(vals, families, nfit, timeout, pool, rng):
  vals = np.sort(vals)
  fit_vals = vals if len(vals) <= nfit else rng.choice(vals, nfit,
                                                       replace=False)
  pending = {
    f: pool.apply_async(_fit_one, ((f, fit_vals, vals), ))
    for f in families
  }
  results = []
  deadline = time.time() + timeout
  for f, r in pending.items():
    try:
      results.append(r.get(max(0.1, deadline - time.time())))
    except multiprocessing.TimeoutError:
      results.append((f, None, math.inf, timeout))
  results.sort(key=lambda r: r[2])
  return results


def cmd_fit(args):
  fits = load_fits(args.fits)
  families = args.families.split(",") if args.families else DEFAULT_FAMILIES
  rng = np.random.default_rng(args.seed)
  source = args.source
  for path in args.samples:
    name = sample_name(path)
    if name in fits["descriptors"] and not args.replace:
      log(f"Skipping {name}: already fit (use --replace)")
      continue
    vals = read_sample(path)
    if len(vals) == 0:
      log(f"Skipping {name}: no finite values in {path}")
      continue
    # a fresh pool per descriptor so fits that time out don't hold workers
    pool = multiprocessing.Pool(args.nprocs or None)
    try:
      results = fit_sample(vals, families, args.nfit, args.timeout, pool, rng)
    finally:
      pool.terminate()
      pool.join()
    best = results[0]
    entry = {
      "min": float(vals.min()),
      "max": float(vals.max()),
      "mean": float(vals.mean()),
      "std": float(vals.std()),
    }
    if best[1] is None:
      log(f"No family could be fit to {name}")
      continue
    entry.update(dist=best[0], params=best[1])
    entry["source"] = source or (
      f"fit to {os.path.basename(path)} (n={len(vals)}, KS={best[2]:.4g}, "
      f"scipy {scipy.__version__})")
    log(f"{name}: " + ", ".join(f"{r[0]} {r[2]:.4g}" for r in results[:3]))
    fits["descriptors"][name] = entry
    # save as we go, fitting everything can take a while
    save_fits(fits, args.fits)


# ---------------------------------------------------------------------------
# table


def build_grid(cdf, lo, hi, dist, params, tol, max_points):
  """Points (x, cdf(x)) on [lo, hi] such that linear interpolation between
  neighbours stays within tol of the CDF (checked at the quarter points of
  every interval, failing intervals are bisected)."""
  if not hi > lo:
    return np.array([lo]), cdf([lo])
  span = hi - lo
  seeds = [np.linspace(lo, hi, 33)]
  # cluster points where the probability mass is, whatever the scale
  q = np.linspace(0, 1, 65)[1:-1]
  with np.errstate(all="ignore"):
    try:
      ppf = dist.ppf(q, *params[:-2], loc=params[-2], scale=params[-1])
      seeds.append(ppf[np.isfinite(ppf) & (ppf > lo) & (ppf < hi)])
    except Exception:
      pass
  # CDFs are often steep right at the clip limits
  offsets = span * np.logspace(-9, -1, 20)
  seeds += [lo + offsets, hi - offsets]
  xs = np.unique(np.concatenate(seeds))
  xs = xs[(xs >= lo) & (xs <= hi)]
  ys = cdf(xs)

  fracs = np.array([0.25, 0.5, 0.75])
  while True:
    x0, x1 = xs[:-1], xs[1:]
    y0, y1 = ys[:-1], ys[1:]
    tx = x0[:, None] + (x1 - x0)[:, None] * fracs
    ty = cdf(tx)
    interp = y0[:, None] + (y1 - y0)[:, None] * fracs
    err = np.abs(ty - interp).max(axis=1)
    # stop splitting intervals that are already at float resolution
    bad = (err > tol) & ((x1 - x0) > 4 * np.spacing(np.maximum(abs(x0),
                                                               abs(x1))))
    if not bad.any():
      break
    if len(xs) + bad.sum() > max_points:
      raise RuntimeError(f"grid needs more than {max_points} points")
    xs = np.concatenate([xs, tx[bad, 1]])
    ys = np.concatenate([ys, ty[bad, 1]])
    order = np.argsort(xs, kind="stable")
    xs, ys = xs[order], ys[order]
  return xs, ys


def table_entries(fits, tol, max_points):
  for name, fit in fits["descriptors"].items():
    cdf = make_cdf(fit)
    dist = scipy_dist(fit["dist"])
    xs, ys = build_grid(cdf, fit["min"], fit["max"], dist, fit["params"], tol,
                        max_points)
    if not np.all(np.isfinite(ys)):
      raise RuntimeError(f"{name}: CDF is not finite on [{fit['min']}, "
                         f"{fit['max']}]")
    if np.any(np.diff(ys) < -1e-12):
      log(f"warning: {name}: CDF decreases by up to "
          f"{-np.diff(ys).min():.3g} (scipy numerical noise)")
    yield name, fit, xs, ys


def cmd_table(args):
  fits = load_fits(args.fits)
  lines = [
    f"# NormalizedDescriptors CDF table, format version "
    f"{TABLE_FORMAT_VERSION}",
    "# Generated by tools/fit_normalized_descriptors.py from "
    f"{os.path.basename(args.fits)}; do not edit.",
    f"# scipy {scipy.__version__}, numpy {np.__version__}, "
    f"interpolation tolerance {args.tol:g}",
    "#",
    "# normalized(v) = linear interpolation of the (x, cdf) points below at",
    "# clip(v, min, max). Points are sorted by x, the first is at min and the",
    "# last at max. A descriptor with min == max has a single point.",
    "#",
    "# descriptor <name> <scipy distribution> <min> <max> <npoints>",
    "# <x> <cdf>   (npoints lines)",
  ]
  total = 0
  for name, fit, xs, ys in table_entries(fits, args.tol, args.max_points):
    lines.append(f"descriptor {name} {fit['dist']} {fit['min']!r} "
                 f"{fit['max']!r} {len(xs)}")
    # x keeps full precision (points can be very close to min/max), the
    # cdf only needs to be well below the interpolation tolerance
    lines.extend(f"{float(x)!r} {float(y):.9g}" for x, y in zip(xs, ys))
    total += len(xs)
    if args.verbose:
      log(f"{name}: {len(xs)} points")
  with open(args.out, "w") as f:
    f.write("\n".join(lines))
    f.write("\n")
  log(f"Wrote {len(fits['descriptors'])} descriptors, {total} points to "
      f"{args.out}")


# ---------------------------------------------------------------------------
# check


def read_table(path):
  table = {}
  with open(path) as f:
    it = (l for l in f if l.strip() and not l.startswith("#"))
    for line in it:
      tag, name, dist, lo, hi, n = line.split()
      assert tag == "descriptor", line
      pts = np.array([[float(v) for v in next(it).split()]
                      for _ in range(int(n))])
      table[name] = (float(lo), float(hi), pts[:, 0], pts[:, 1])
  return table


def table_value(entry, v):
  lo, hi, xs, ys = entry
  return np.interp(np.clip(v, lo, hi), xs, ys)


def cmd_check(args):
  fits = load_fits(args.fits)
  table = read_table(args.table)
  rng = np.random.default_rng(0)
  worst = 0.0
  missing = sorted(set(fits["descriptors"]) - set(table))
  if missing:
    log(f"In the fits file but not the table: {' '.join(missing)}")
  for name, fit in fits["descriptors"].items():
    if name not in table:
      continue
    cdf = make_cdf(fit)
    lo, hi = fit["min"], fit["max"]
    pad = max(hi - lo, 1.0) * 0.1
    probes = [rng.uniform(lo - pad, hi + pad, args.n)]
    if hi - lo < 1e5:
      # count descriptors only ever see integers
      probes.append(np.arange(math.floor(lo), math.ceil(hi) + 1))
    sample = os.path.join(args.samples, f"d_{name}.gz") if args.samples else ""
    if os.path.exists(sample):
      probes.append(read_sample(sample))
    v = np.concatenate(probes)
    err = float(np.abs(table_value(table[name], v) - cdf(v)).max())
    worst = max(worst, err)
    if args.verbose or err > args.tol:
      log(f"{name}: max error {err:.3g}" +
          (" (sample checked)" if os.path.exists(sample) else ""))
  log(f"Checked {len(table)} descriptors, worst error {worst:.3g}")
  if worst > args.tol:
    sys.exit(1)


# ---------------------------------------------------------------------------
# histograms (method 2, no fitting)


def make_histogram(vals):
  """descriptastorus data/d_descriptors/make_histdists.py, unchanged:
  (left bin edge, cumulative fraction of values up to and including that
  bin) for min(1000, distinct values) equal-width bins."""
  vals = np.asarray(vals, dtype=float)
  vals = vals[np.isfinite(vals)]
  n = min(MAX_HIST_BINS, len(set(vals)))
  hist, xaxis = np.histogram(vals, bins=n)
  total = hist.sum()
  bins = []
  last = 0.0
  for value, x in zip(hist, xaxis):
    last += value
    bins.append((float(x), float(last / total)))
  return bins


def hist_lookup(bins, v):
  """descriptastorus rdNormalizedDescriptors.histcdf, unchanged."""
  p = bisect.bisect(bins, (v, ))
  if p < len(bins):
    return bins[p][1]
  return 1.0


def read_hist_table(path):
  hists = {}
  if not os.path.exists(path):
    return hists
  with open(path) as f:
    it = (l for l in f if l.strip() and not l.startswith("#"))
    for line in it:
      tag, name, n = line.split()
      assert tag == "histogram", line
      bins = []
      for _ in range(int(n)):
        x, c = next(it).split()
        bins.append((float(x), float(c)))
      hists[name] = bins
  return hists


def write_hist_table(hists, path):
  lines = [
    "# NormalizedDescriptors histogram table, format version "
    f"{TABLE_FORMAT_VERSION}",
    "# Generated by tools/fit_normalized_descriptors.py; do not edit.",
    "#",
    "# Raw-data histograms (no fitted distribution), as descriptastorus'",
    "# RDKit2DHistogramNormalized. For each descriptor, nbins lines of",
    "# <left edge> <cumulative fraction>, edges ascending. normalized(v) is",
    "# the fraction of the first bin whose edge is >= v, or 1.0 when v is",
    "# greater than every edge (python: bins[bisect(bins, (v,))][1]).",
    "# No interpolation and no clipping.",
    "#",
    "# histogram <name> <nbins>",
    "# <edge> <cumulative fraction>   (nbins lines)",
  ]
  for name in sorted(hists):
    bins = hists[name]
    lines.append(f"histogram {name} {len(bins)}")
    lines.extend(f"{x!r} {c!r}" for x, c in bins)
  with open(path, "w") as f:
    f.write("\n".join(lines))
    f.write("\n")
  log(f"Wrote {len(hists)} histograms to {path}")


def read_hists_py(path):
  """Read a descriptastorus hists.py without importing descriptastorus."""
  ns = {"inf": math.inf, "nan": math.nan}
  with open(path) as f:
    exec(f.read(), ns)
  return ns["hists"]


def cmd_import_hists(args):
  hists = read_hist_table(args.table)
  src = read_hists_py(args.hists_py)
  n = 0
  for name, bins in src.items():
    if name in hists and not args.replace:
      continue
    hists[name] = [(float(x), float(c)) for x, c in bins]
    n += 1
  log(f"Imported {n} of {len(src)} histograms from {args.hists_py}")
  write_hist_table(hists, args.table)


def cmd_histogram(args):
  hists = read_hist_table(args.table)
  for path in args.samples:
    name = sample_name(path)
    if name in hists and not args.replace:
      log(f"Skipping {name}: already has a histogram (use --replace)")
      continue
    vals = read_sample(path)
    if len(vals) == 0:
      log(f"Skipping {name}: no finite values in {path}")
      continue
    hists[name] = make_histogram(vals)
    log(f"{name}: {len(hists[name])} bins from {len(vals)} values")
  write_hist_table(hists, args.table)


def cmd_check_hists(args):
  table = read_hist_table(args.table)
  ref = read_hists_py(args.hists_py)
  rng = np.random.default_rng(0)
  bad = 0
  for name, bins in ref.items():
    if name not in table:
      log(f"{name}: in {args.hists_py} but not the table")
      bad += 1
      continue
    edges = np.array([b[0] for b in bins])
    probes = np.concatenate([
      edges, edges + 1e-9,
      rng.uniform(edges[0] - 1, edges[-1] + 1, 2000)
    ])
    diff = max(
      abs(hist_lookup(table[name], v) - hist_lookup(bins, v)) for v in probes)
    if diff:
      log(f"{name}: lookups differ by up to {diff}")
      bad += 1
  log(f"Checked {len(ref)} histograms, {bad} differ")
  if bad:
    sys.exit(1)


# ---------------------------------------------------------------------------


def main(argv=None):
  p = argparse.ArgumentParser(
    description=__doc__.split("\n\n")[0],
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__.split("\n\n", 1)[1])
  p.add_argument("--fits", default=DEFAULT_FITS,
                 help="fits file (default: %(default)s)")
  sub = p.add_subparsers(dest="cmd", required=True)

  s = sub.add_parser("import-dists", help="import a descriptastorus dists.py")
  s.add_argument("dists_py")
  s.add_argument("--replace", action="store_true",
                 help="overwrite descriptors already in the fits file")
  s.add_argument("--source", help="provenance note stored with each fit")
  s.set_defaults(func=cmd_import_dists)

  s = sub.add_parser("compute", help="compute descriptor samples")
  s.add_argument("descriptors", nargs="+", help="RDKit descriptor names")
  s.add_argument("--smiles", required=True,
                 help="SMILES file, e.g. descriptastorus/data/chembl_100k.smi")
  s.add_argument("--out", default="samples",
                 help="directory for d_<name>.gz files (default: %(default)s)")
  s.add_argument("--nprocs", type=int, default=0)
  s.set_defaults(func=cmd_compute)

  s = sub.add_parser("missing", help="list RDKit descriptors without a fit")
  s.add_argument("--histograms", action="store_true",
                 help="list those without a histogram instead")
  s.add_argument("--table", default=DEFAULT_HISTS,
                 help="histogram table for --histograms")
  s.set_defaults(func=cmd_missing)

  s = sub.add_parser("fit", help="fit distributions to sample files")
  s.add_argument("samples", nargs="+", help="d_<name>.gz sample files")
  s.add_argument("--families",
                 help="comma separated scipy.stats distributions to try "
                 "(default: the families descriptastorus uses)")
  s.add_argument("--nfit", type=int, default=20000,
                 help="fit on a random subset of this many values; the KS "
                 "distance used to rank fits always uses every value "
                 "(default: %(default)s)")
  s.add_argument("--timeout", type=float, default=300,
                 help="seconds allowed for all families on one descriptor "
                 "(default: %(default)s)")
  s.add_argument("--seed", type=int, default=42)
  s.add_argument("--replace", action="store_true",
                 help="refit descriptors already in the fits file")
  s.add_argument("--source", help="provenance note stored with each fit")
  s.add_argument("--nprocs", type=int, default=0)
  s.set_defaults(func=cmd_fit)

  s = sub.add_parser("table", help="write the CDF interpolation table")
  s.add_argument("--out", default=DEFAULT_TABLE,
                 help="table file (default: %(default)s)")
  s.add_argument("--tol", type=float, default=1e-5,
                 help="max interpolation error; 1e-6 roughly triples the "
                 "table size (default: %(default)s)")
  s.add_argument("--max-points", type=int, default=200000)
  s.add_argument("-v", "--verbose", action="store_true")
  s.set_defaults(func=cmd_table)

  s = sub.add_parser("check", help="check the table against scipy")
  s.add_argument("--table", default=DEFAULT_TABLE)
  s.add_argument("--samples",
                 help="directory of d_<name>.gz files to use as probe values")
  s.add_argument("--n", type=int, default=100000,
                 help="random probe values per descriptor")
  s.add_argument("--tol", type=float, default=1.1e-5)
  s.add_argument("-v", "--verbose", action="store_true")
  s.set_defaults(func=cmd_check)

  s = sub.add_parser("import-hists",
                     help="import a descriptastorus hists.py (method 2)")
  s.add_argument("hists_py")
  s.add_argument("--table", default=DEFAULT_HISTS)
  s.add_argument("--replace", action="store_true")
  s.set_defaults(func=cmd_import_hists)

  s = sub.add_parser("histogram",
                     help="add raw-data histograms of sample files (method 2)")
  s.add_argument("samples", nargs="+", help="d_<name>.gz sample files")
  s.add_argument("--table", default=DEFAULT_HISTS)
  s.add_argument("--replace", action="store_true")
  s.set_defaults(func=cmd_histogram)

  s = sub.add_parser("check-hists",
                     help="check histogram table lookups against a hists.py")
  s.add_argument("hists_py")
  s.add_argument("--table", default=DEFAULT_HISTS)
  s.set_defaults(func=cmd_check_hists)

  args = p.parse_args(argv)
  args.func(args)


if __name__ == "__main__":
  main()
