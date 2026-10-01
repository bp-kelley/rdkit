# Fitting the normalized descriptors

`fit_normalized_descriptors.py` maintains the data behind the
NormalizedDescriptors extension, a port of descriptastorus' `RDKit2DNormalized`
descriptors. Each descriptor value `v` is normalized as

    clip(dist.cdf(clip(v, min, max), *shape, loc, scale), 0, 1)

where `dist` is a `scipy.stats` distribution fitted to the descriptor over a
reference set of molecules.

## Files

| file | what it is |
| --- | --- |
| `../data/normalized_descriptor_fits.json` | One fit per descriptor: scipy family, parameters (`shape..., loc, scale`), clip `min`/`max`, sample `mean`/`std`, and where the fit came from. This is the source of truth. |
| `../data/normalized_descriptor_cdfs.txt` | Generated. Each fitted CDF evaluated on an adaptive grid so that linear interpolation stays within `1e-5` of scipy. This is what the C++ code reads. |

The first 201 fits are imported unchanged from descriptastorus'
`descriptastorus/descriptors/dists.py`, so the C++ results match
`RDKit2DNormalized`.

## Table format

    # comment lines start with #
    descriptor <name> <scipy distribution> <min> <max> <npoints>
    <x> <cdf>
    ... (npoints lines, x ascending, first x == min, last x == max)

To normalize `v`: clip it to `[min, max]`, find the bracketing points and
interpolate linearly. When `min == max` there is a single point and the result
is that constant.

## Regenerating the table

    python fit_normalized_descriptors.py table
    python fit_normalized_descriptors.py check --samples /path/to/descriptastorus/data/d_descriptors

`check` compares the interpolated table with scipy on random values, on every
integer in range, and on the reference samples when given.

## Adding new or missing descriptors

The reference set is descriptastorus' `data/chembl_100k.smi`.

    # RDKit descriptors that don't have a fit yet
    python fit_normalized_descriptors.py missing

    # compute them (one d_<name>.gz per descriptor, the same format as
    # descriptastorus' data/d_descriptors)
    python fit_normalized_descriptors.py compute \
        --smiles /path/to/descriptastorus/data/chembl_100k.smi \
        --out samples $(python fit_normalized_descriptors.py missing)

    # fit each sample and store the best fit in the fits file
    python fit_normalized_descriptors.py fit samples/d_*.gz

    python fit_normalized_descriptors.py table
    python fit_normalized_descriptors.py check --samples samples

`fit` tries the 45 scipy families descriptastorus already uses (pass
`--families` to change that) on a random subset of 20000 values, then picks the
family with the smallest distance between its CDF and the fraction of values
below each observed value (the KS statistic without the jump term, so count
descriptors are judged at the integers they take). If no family gets within
`--max-ks` (default 0.05), the empirical CDF of the sample is stored instead
(`"dist": "empirical"` with `knots`), which suits multimodal descriptors
such as BCUT2D. Use
`--replace` to refit a descriptor that already has a fit.

Requires numpy and scipy; `compute` also needs the RDKit python wrappers.
