# Building the normalized descriptor tables

`fit_normalized_descriptors.py` maintains the data behind the
NormalizedDescriptors extension, a port of descriptastorus' normalized
descriptors. descriptastorus has two separate methods, and so do these tables.
They are never mixed: a descriptor's fitted value never comes from its data
histogram, and its histogram value never comes from a fit.

1. **Fitted distribution** (`RDKit2DNormalized`). Each value `v` becomes

       clip(dist.cdf(clip(v, min, max), *shape, loc, scale), 0, 1)

   where `dist` is the `scipy.stats` distribution that best fits the
   descriptor over a reference set of molecules.

2. **Histogram** (`RDKit2DHistogramNormalized`). No fitting. The raw values
   are binned with `numpy.histogram` exactly as descriptastorus'
   `data/d_descriptors/make_histdists.py` does, and `v` is looked up in the
   cumulative bins exactly as descriptastorus does.

## Files

| file | method | what it is |
| --- | --- | --- |
| `../data/normalized_descriptor_fits.json` | fitted | One fit per descriptor: scipy family, parameters (`shape..., loc, scale`), clip `min`/`max`, sample `mean`/`std`, and where the fit came from. The source of truth for method 1. |
| `../data/normalized_descriptor_cdfs.txt` | fitted | Generated from the fits. Each fitted CDF evaluated on an adaptive grid so that linear interpolation stays within `1e-5` of scipy. |
| `../data/normalized_descriptor_histograms.txt` | histogram | The cumulative histogram bins. |

The first 201 fits are imported unchanged from descriptastorus'
`descriptastorus/descriptors/dists.py`, and the first 201 histograms unchanged
from `descriptastorus/descriptors/hists.py`, so the C++ results match
descriptastorus.

## Fitted table format

    # comment lines start with #
    descriptor <name> <scipy distribution> <min> <max> <npoints>
    <x> <cdf>
    ... (npoints lines, x ascending, first x == min, last x == max)

To normalize `v`: clip it to `[min, max]`, find the bracketing points and
interpolate linearly. When `min == max` there is a single point and the result
is that constant.

## Histogram table format

    # comment lines start with #
    histogram <name> <nbins>
    <left edge> <cumulative fraction>
    ... (nbins lines, edges ascending)

To normalize `v`: the cumulative fraction of the first bin whose edge is
`>= v`, or `1.0` when `v` is greater than every edge. In python this is
`bins[bisect.bisect(bins, (v,))][1]`. There is no interpolation and no
clipping.

## Regenerating the tables

    python fit_normalized_descriptors.py table
    python fit_normalized_descriptors.py check --samples /path/to/descriptastorus/data/d_descriptors
    python fit_normalized_descriptors.py check-hists /path/to/descriptastorus/descriptastorus/descriptors/hists.py

`check` compares the interpolated table with scipy on random values, on every
integer in range, and on the reference samples when given. `check-hists`
compares lookups in the histogram table with lookups in a descriptastorus
`hists.py`.

## Adding new or missing descriptors

The reference set is descriptastorus' `data/chembl_100k.smi`.

    # RDKit descriptors that don't have a fit yet
    # (add --histograms for those without a histogram)
    python fit_normalized_descriptors.py missing

    # compute them (one d_<name>.gz per descriptor, the same format as
    # descriptastorus' data/d_descriptors)
    python fit_normalized_descriptors.py compute \
        --smiles /path/to/descriptastorus/data/chembl_100k.smi \
        --out samples $(python fit_normalized_descriptors.py missing)

    # method 1: fit each sample, then regenerate the interpolation table
    python fit_normalized_descriptors.py fit samples/d_*.gz
    python fit_normalized_descriptors.py table
    python fit_normalized_descriptors.py check --samples samples

    # method 2: bin each sample into the histogram table
    python fit_normalized_descriptors.py histogram samples/d_*.gz

`fit` tries the 45 scipy families descriptastorus already uses (pass
`--families` to change that) on a random subset of 20000 values, then keeps the
family with the smallest distance between its CDF and the fraction of values
below each observed value (the KS statistic without the jump term, so count
descriptors are judged at the integers they take). Use `--replace` to refit a
descriptor that already has a fit, or to rebuild a histogram.

Requires numpy and scipy; `compute` also needs the RDKit python wrappers.
