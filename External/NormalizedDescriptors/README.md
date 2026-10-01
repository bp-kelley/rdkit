# NormalizedDescriptors

A C++ port of the `RDKit2DNormalized` descriptors from
[descriptastorus](https://github.com/bp-kelley/descriptastorus).

Each descriptor value is mapped to the fraction of a reference population
(the descriptastorus fits were made on ~100k compounds) with a smaller value,
using a tabulated cumulative distribution function (CDF):

1. clip the value to the descriptor's `[minV, maxV]`,
2. linearly interpolate the CDF table,
3. clip the result to `[0, 1]`.

Descriptors without a table, and non-finite values, normalize to `0.0`, as in
descriptastorus.

The extension is built when `RDK_BUILD_NORMALIZED_DESCRIPTORS` is `ON` (the
default).

## Tables

`data/rdkit2d_cdf_v1.txt` holds one table per descriptor in the 201-entry
`dists.py` of descriptastorus. It was produced by sampling each fitted
`scipy.stats` CDF with `Scripts/dists_to_cdf_tables.py`:

```
python Scripts/dists_to_cdf_tables.py path/to/descriptastorus/descriptors/dists.py \
    data/rdkit2d_cdf_v1.txt
```

Linear interpolation of these tables matches the scipy CDFs to within ~0.004
everywhere (median max error ~1e-4), and exactly at integer values for
descriptors with a small range (the counts).

The table format is one line per descriptor, whitespace separated, with `#`
comment lines:

```
name minV maxV npts x_0 cdf_0 x_1 cdf_1 ... x_{npts-1} cdf_{npts-1}
```

## Usage

C++:

```c++
#include <GraphMol/NormalizedDescriptors/NormalizedDescriptors.h>

using namespace RDKit::NormalizedDescriptors;
const auto &tables = getDefaultTables();  // reads $RDBASE/External/...
double v = tables.normalize("MolLogP", 2.5);
```

Python:

```python
from rdkit.Chem import rdNormalizedDescriptors as rdnd
rdnd.NormalizeDescriptor('MolLogP', 2.5)
tables = rdnd.CDFTableSet()
tables.LoadFromFile('my_tables.txt')
tables.Normalize('MolLogP', 2.5)
```
