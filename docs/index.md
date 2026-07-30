# IsoGen

IsoGen generates isotope distributions for protein, RNA, and DNA analytes. It
accepts either a neutral mass or a sequence and returns a two-column NumPy
array containing neutral mass and relative intensity.

```python
import isogen

distribution = isogen.isodist(
    "ACDEFGHIK",
    type="PEPTIDE",
    isolen=64,
    method="FFT",
)
```

The first mass is the monoisotopic mass for sequence input. Subsequent masses
use a configurable isotope spacing: 1.0033 Da for peptides and 1.0027 Da for
nucleic acids by default.

## Features

- FFT and neural-network isotope intensity calculations
- Neutral-mass and sequence input
- Protein, RNA, and DNA mass axes
- Intact proteins and a, b, c, x, y, and z fragment ions
- Configurable RNA and DNA terminal chemistry
- A command-line interface for CSV output and example stick plots

!!! note "DNA intensity model"
    DNA sequence intensities currently use the RNA sequence model after
    thymine is mapped to uracil. The returned mass axis still uses the
    DNA-specific residue masses.

Continue with [Getting started](getting-started.md), or go directly to the
[API reference](api.md).
