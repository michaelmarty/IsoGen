# Getting started

## Installation

Install a published wheel:

```shell
python -m pip install isogen
```

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

The result has shape `(isolen, 2)`. Column zero is neutral mass and column one
is relative intensity:

```python
masses = protein[:, 0]
intensities = protein[:, 1]
```

`type` and `method` currently use uppercase names. Supported values are
`PEPTIDE`, `RNA`, and `DNA`, and `FFT` and `NN`, respectively.

## Choose the calculation engine

FFT is the default:

```python
fft_result = isogen.isodist("PEPTIDE", method="FFT")
nn_result = isogen.isodist("PEPTIDE", method="NN")
```

FFT performs the direct isotope calculation. The neural-network engine uses
the packaged pretrained model and can be useful for rapid approximation.
