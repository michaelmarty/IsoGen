# Getting started

## Installation

Install a published wheel:

```shell
python -m pip install isogen
```

IsoGen requires Python 3.13 or newer.

IsoGen publishes precompiled native libraries for 64-bit Windows and Linux.
Linux installations require the FFTW 3 runtime; released manylinux wheels
bundle that dependency. Other platforms require a native build from source.

## Generate distributions

Use a numeric neutral mass with any analyte model:

```python
import isogen

protein_by_mass = isogen.isodist(10_000, type="PEPTIDE")
rna_by_mass = isogen.isodist(10_000, type="RNA")
```

For sequence input, IsoGen calculates the monoisotopic first mass:

```python
protein = isogen.isodist("ACDEFGHIK", type="PEPTIDE", isolen=64)
rna = isogen.isodist("AUGCAGUACGUA", type="RNA", isolen=64)
dna = isogen.isodist("ATGCAGTACGTA", type="DNA", isolen=64)
```

For an elemental formula, select the FFT-only `ATOM` input type:

```python
glucose = isogen.isodist("C6H12O6", type="ATOM", isolen=32)
```

The result has shape `(isolen, 2)`. Column zero is neutral mass and column one
is relative intensity:

```python
masses = protein[:, 0]
intensities = protein[:, 1]
```

Supported `type` values are `PEPTIDE`, `RNA`, `DNA`, and `ATOM`; names are
case-insensitive in the Python API. `FFT`, `BRAIN`, and `NN` are available for
biopolymer inputs. `ATOM` supports `FFT` only.

## Peptide ion types

Peptide sequence input represents an intact neutral peptide by default, using
`ion_type="H2O"`. For a fragment, pass only the residues present in that
fragment and choose one of the supported ion types:

| Ion type | Sequence to pass | Neutral terminal composition |
| --- | --- | --- |
| `H2O` | Intact peptide | +H2O |
| `a` | N-terminal fragment | -CO |
| `b` | N-terminal fragment | No terminal shift |
| `c` | N-terminal fragment | +NH3 |
| `x` | C-terminal fragment | +CO2 |
| `y` | C-terminal fragment | +H2O |
| `z` | C-terminal fragment | +H2O-NH3 |

For example, split `PEPTIDE` into the appropriate N- or C-terminal sequence:

```python
b6 = isogen.isodist(
    "PEPTID",
    type="PEPTIDE",
    ion_type="b",
)
y6 = isogen.isodist(
    "EPTIDE",
    type="PEPTIDE",
    ion_type="y",
)
```

IsoGen returns neutral masses. It does not add protons, assign charge, or
convert these values to m/z.

## RNA ions and terminal chemistry

`ion_type` is a peptide-only option. Named RNA fragmentation series such as
a/b/c/d and w/x/y/z are not currently calculated automatically. To describe
an intact or manually truncated RNA sequence, use the available terminal
chemistry options:

| Keyword | Supported values |
| --- | --- |
| `threeend` | `OH` (default), or no terminal adjustment |
| `fiveend` | `OH`, `MP` (default), or `TP` |

```python
rna_5_hydroxyl = isogen.isodist(
    "AUGC",
    type="RNA",
    threeend="OH",
    fiveend="OH",
)
rna_5_triphosphate = isogen.isodist(
    "AUGC",
    type="RNA",
    threeend="OH",
    fiveend="TP",
)
```

These options are terminal mass adjustments; they are not aliases for named
RNA fragment-ion series.

!!! important
    `ion_type`, `threeend`, and `fiveend` are forwarded to the mass-axis
    calculation. They change the monoisotopic origin, but the FFT, BRAIN, or NN
    sequence-model intensity vector retains its standard terminal
    composition.

## Choose the calculation engine

FFT is the default:

```python
fft_result = isogen.isodist("PEPTIDE", method="FFT")
brain_result = isogen.isodist("PEPTIDE", method="BRAIN")
nn_result = isogen.isodist("PEPTIDE", method="NN")
```

FFT performs the direct isotope calculation, while BRAIN uses a polynomial
recurrence. The neural-network engine uses the packaged pretrained model and
can be useful for rapid approximation.
Elemental formulas always use the direct FFT calculation.
