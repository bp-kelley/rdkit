#!/usr/bin/env python
#
#  Copyright (C) 2026 The RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Convert descriptastorus-style scipy distribution fits into CDF tables.

The input is a python file defining a ``dists`` dictionary in the format used
by descriptastorus (descriptastorus/descriptors/dists.py)::

  dists = {
    name: (scipy_dist_name, (shape..., loc, scale), minV, maxV, mean, std),
    ...
  }

For every descriptor the fitted scipy CDF is evaluated on a grid covering
[minV, maxV] and written as one line of the table format read by
RDKit::NormalizedDescriptors::CDFTableSet::loadFromStream::

  name minV maxV npts x_0 cdf_0 x_1 cdf_1 ... x_{npts-1} cdf_{npts-1}

The grid is the union of evenly spaced points, log-spaced points close to minV,
the quantiles of the fitted distribution, and every integer in the range when the range is small, so that
linear interpolation reproduces the scipy CDF closely both in the dense part
of the distribution and in long tails.

Usage::

  python dists_to_cdf_tables.py path/to/dists.py output.txt [--npts 256]
"""
import argparse
import math
import runpy
import sys

import numpy as np
import scipy.stats as st


def get_dist(name):
  if name in ('gilbrat', 'gibrat'):
    # scipy renamed gilbrat -> gibrat
    name = 'gibrat' if hasattr(st, 'gibrat') else 'gilbrat'
  return getattr(st, name)


def make_grid(dist, args, loc, scale, minV, maxV, npts):
  if maxV <= minV:
    return np.array([minV], dtype=float)
  span = maxV - minV
  # evenly spaced points plus points that are dense close to minV, where many
  # of the fitted distributions (e.g. wald for the fr_ counts) are very steep
  pts = [np.linspace(minV, maxV, npts), minV + span * np.logspace(-9, 0, npts // 2)]
  with np.errstate(all='ignore'):
    qs = dist.ppf(np.linspace(0, 1, npts + 1)[1:-1], *args, loc=loc, scale=scale)
  qs = np.asarray(qs, dtype=float)
  pts.append(qs[np.isfinite(qs)])
  if maxV - minV <= 4 * npts:
    pts.append(np.arange(math.ceil(minV), math.floor(maxV) + 1, dtype=float))
  grid = np.unique(np.clip(np.concatenate(pts), minV, maxV))
  return grid


def make_table(dist_name, params, minV, maxV, npts):
  dist = get_dist(dist_name)
  args = params[:-2]
  loc = params[-2]
  scale = params[-1]
  grid = make_grid(dist, args, loc, scale, minV, maxV, npts)
  with np.errstate(all='ignore'):
    cdf = dist.cdf(grid, *args, loc=loc, scale=scale)
  cdf = np.clip(np.nan_to_num(cdf, nan=0.0), 0.0, 1.0)
  # guard against tiny numerical non-monotonicity in some scipy cdfs
  cdf = np.maximum.accumulate(cdf)
  return grid, cdf


def write_tables(dists, out, npts):
  out.write('# RDKit normalized descriptor CDF tables\n')
  out.write('# name minV maxV npts x_0 cdf_0 ... x_{npts-1} cdf_{npts-1}\n')
  for name in sorted(dists):
    dist_name, params, minV, maxV = dists[name][:4]
    grid, cdf = make_table(dist_name, params, float(minV), float(maxV), npts)
    fields = [name, repr(float(minV)), repr(float(maxV)), str(len(grid))]
    for x, c in zip(grid, cdf):
      fields.append(f'{x:.10g}')
      fields.append(f'{c:.8g}')
    out.write(' '.join(fields))
    out.write('\n')


def main(argv=None):
  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('dists', help='python file defining a "dists" dictionary')
  parser.add_argument('output', help='output table file ("-" for stdout)')
  parser.add_argument('--npts', type=int, default=256,
                      help='number of evenly spaced and quantile grid points (default 256)')
  opts = parser.parse_args(argv)
  dists = runpy.run_path(opts.dists)['dists']
  if opts.output == '-':
    write_tables(dists, sys.stdout, opts.npts)
  else:
    with open(opts.output, 'w') as out:
      write_tables(dists, out, opts.npts)


if __name__ == '__main__':
  main()
