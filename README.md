# IsoGen
IsoGen is a toolbox for predicting isotope distributions in mass spectrometry data.

It includes both an absolute FFT-based calculation and a neural network prediction.

Pretrained models are included for both peptides and RNA based on either average mass or sequence.

The FFT methods are absolute and are limited only by the accuracy of the data you put in. They are a little faster, especially on larger species.

The NN methods are very accurate and can be faster on smaller species. The primary advantage of these is that they can be retrained on non-standard isotope distributions.

## Installation

Install a published wheel from PyPI:

```shell
python -m pip install isogen
```

Precompiled native libraries are provided for 64-bit Windows and Linux. Linux
requires the FFTW 3 runtime; published Linux wheels bundle it during the
manylinux repair step. For other platforms, build the native library from
source using CMake.

## Usage

A simple native test executable is provided. It can be run with:

```shell
isogen_test.exe -mass 10000
```

This will create an isotope distribution printout for a mass of 10 kDa.

From Python:

```python
import isogen

protein = isogen.isodist("ACDEFGHIK", type="PEPTIDE", isolen=64)
rna = isogen.isodist("AUGCAGUACGUA", type="RNA", isolen=64)
dna = isogen.isodist("ATGCAGTACGTA", type="DNA", isolen=64)
```

### Peptide ions and RNA termini

For peptide fragments, pass the fragment sequence and select its neutral
terminal composition with `ion_type`. IsoGen supports intact `H2O` (the
default) and the peptide `a`, `b`, `c`, `x`, `y`, and `z` ion types:

```python
b6 = isogen.isodist("PEPTID", type="PEPTIDE", ion_type="b")
y6 = isogen.isodist("EPTIDE", type="PEPTIDE", ion_type="y")
```

Supply the N-terminal subsequence for a/b/c ions and the C-terminal subsequence
for x/y/z ions. Returned values are neutral masses, not charge-adjusted m/z.

RNA does not currently accept named RNA fragment-ion series through
`ion_type`. For an intact or manually truncated RNA sequence, configure the
supported terminal chemistry with `threeend` and `fiveend`:

```python
rna_5_triphosphate = isogen.isodist(
    "AUGC",
    type="RNA",
    threeend="OH",
    fiveend="TP",
)
```

The available 5' settings are hydroxyl (`OH`), monophosphate (`MP`, default),
and triphosphate (`TP`); the supported explicit 3' setting is hydroxyl (`OH`,
default). These peptide-ion and RNA-terminal options adjust the mass-axis
origin. The sequence-model intensity vector retains its standard terminal
composition.

From the command line:

```shell
isogen dist ACDEFGHIK --type PEPTIDE --isolen 64
isogen plot
```

See `python -m isogen --help` for all options.

## Documentation

The full user guide and API reference are in the
[`docs`](docs/index.md) directory. To preview them locally:

```shell
python -m pip install -e ".[docs]"
python -m mkdocs serve
```

## Tests

The test suite uses Pyteomics as an independent mass reference. Pyteomics is
only part of the optional test dependencies and is not installed with IsoGen:

```shell
python -m pip install -e ".[test]"
python -m pytest
```

## License

IsoGen is released under the BSD 3-Clause License. See LICENSE for details.

PLEASE CITE THIS SOFTWARE IN ANY PUBLICATIONS THAT USE IT (publication to follow).

## Contact

If you have any questions, please email mtmarty@utexas.edu or open a ticket on GitHub.

TODO:

