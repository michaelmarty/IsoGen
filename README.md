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

TODO: Fix isojim and other import issues. Remap models to be in subfolder to clean up the isogen folder.



