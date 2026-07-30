# IsoGen
IsoGen is a toolbox for predicting isotope distributions from protein, RNA,
DNA, neutral-mass, and elemental-formula inputs.

It includes both an absolute FFT-based calculation and a neural network prediction.

Pretrained models are included for both peptides and RNA based on either average mass or sequence. DNA prediction uses the RNA model due to the similarity of their elemental compositions.

The FFT methods are absolute and are limited only by the accuracy of the data you put in. They are a little faster, especially on larger species.

The NN methods are very accurate and can be faster on smaller species. The primary advantage of these is that they can be retrained on non-standard isotope distributions.

## Installation

Install a published wheel from PyPI:

```shell
python -m pip install pyisogen
```

IsoGen requires Python 3.13 or newer.

Precompiled native libraries are provided for 64-bit Windows and Linux. Linux
requires the FFTW 3 runtime; published Linux wheels bundle it during the
manylinux repair step. For other platforms, build the native library from
source using CMake.

## Usage

From Python:

```python
import isogen

protein = isogen.isodist("ACDEFGHIK", type="PEPTIDE", isolen=64)
rna = isogen.isodist("AUGCAGUACGUA", type="RNA", isolen=64)
dna = isogen.isodist("ATGCAGTACGTA", type="DNA", isolen=64)
glucose_mass_dist = isogen.isodist("C6H12O6", type="ATOM", isolen=32)
```

The output is a numpy array of shape `(isolen, 2)` with the first column containing the monoisotopic mass and the second column containing the relative intensity. The `isolen` parameter controls the number of isotopic peaks returned.

In addition to FFT methods, IsoGen provides neural-network models for peptides and RNA. The default is to use the exact "FFT" model with averagines to calculate the isotope distribution from an average protein, RNA, or DNA mass. Turning on "NN" mode uses the neural-network model to predict the isotope distribution from a peptide or RNA sequence. 

The `PEPTIDE` model is trained on peptide sequences, while the `RNA` model is trained on RNA sequences. The `DNA` type uses the RNA model, and

The public `ATOM` type uses the FFT method; no neural-network formula model is
available.

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
isogen dist C6H12O6 --type ATOM --isolen 32
isogen plot
```

See `python -m isogen --help` for all options.

The source repository also builds a native development executable named
`isogen_test.exe` on Windows (`isogen_test` on Linux). It can be run from the
repository's `bin` directory with `isogen_test.exe -mass 10000`, but it is not
installed by the Python wheel. Use the `isogen` console command for installed
packages.

## Documentation

Read the [full IsoGen documentation](https://michaelmarty.github.io/IsoGen/).
The documentation sources are also available in the repository's `docs`
directory. To preview them locally:

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

Development and model-training modules have additional dependencies:

```shell
python -m pip install -e ".[training]"
```

## License

IsoGen is released under the BSD 3-Clause License. See LICENSE for details.

PLEASE CITE THIS SOFTWARE IN ANY PUBLICATIONS THAT USE IT (publication to follow).

## Contact

If you have any questions, please email mtmarty@utexas.edu or open a ticket on GitHub.


## CHANGELOG

### 1.0.0

Initial release. Rewrote significantly from UniDec build using AI tool to improve the release and add in atomic formula support.

