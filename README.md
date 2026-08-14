# IsoGen
IsoGen is a toolbox for predicting isotope distributions from protein, RNA,
DNA, neutral-mass, and elemental-formula inputs.

It includes absolute FFT and BRAIN calculations plus a neural network prediction.

Pretrained models are included for both peptides and RNA based on either average mass or sequence. DNA prediction uses the RNA model due to the similarity of their elemental compositions.

The FFT methods are absolute and are limited only by the accuracy of the data you put in. They are a little faster, especially on larger species.

The BRAIN method is an absolute calculation based on a polynomial recurrence.
It provides an alternative to FFT for peptide, RNA, and DNA sequence or
neutral-mass inputs.

The NN methods are very accurate and can be faster on smaller species. The primary advantage of these is that they can be retrained on non-standard isotope distributions.

## License

IsoGen is released under the BSD 3-Clause License. See LICENSE for details.

Additional license information is available in THIRD_PARTY_NOTICES.md for FFTW and Intel compiler runtime libraries.

PLEASE CITE THIS SOFTWARE IN ANY PUBLICATIONS THAT USE IT AND RELEVANT BACKGROUND LITERATURE BELOW. Machine-readable citation metadata is available in [CITATION.cff](CITATION.cff).

IsoGen is currently in preprint, please cite:

Pavek, J.G.; Grimes, J.; Frey, B.L.; Welham, N.V.; Smith, L.M.; Marty, M.T. "[Neural Network Prediction of Isotopic Distributions](https://chemrxiv.org/doi/abs/10.26434/chemrxiv.15006709/v1)" ChemRxiv 2026, doi:10.26434/chemrxiv.15006709/v1

The FFT method is derived from the following citation:

Rockwood, A.L.; Palmblad, M. [Isotopic Distributions. In Mass Spectrometry Data Analysis in Proteomics](https://doi.org/10.1007/978-1-62703-392-3_3), Matthiesen, R. Ed.; Humana Press, 2013; pp 65–99.

The BRAIN method is derived from these citations:

Dittwald, P.; Claesen, J.; Burzykowski, T.; Valkenborg, D.; Gambin, A. "[BRAIN: A Universal Tool for High-Throughput Calculations of the Isotopic Distribution for Mass Spectrometry](https://doi.org/10.1021/ac303439m)" Analytical Chemistry 2013, 85, 1991–1994

Dittwald, P.; Valkenborg, D. "[BRAIN 2.0: Time and Memory Complexity Improvements in the Algorithm for Calculating the Isotope Distribution](https://doi.org/10.1007/s13361-013-0796-5)" Journal of the American Society for Mass Spectrometry 2014, 25, 588–594

## Contact

If you have any questions, please email mtmarty@utexas.edu or open a ticket on GitHub.

## Installation

Install a published wheel from PyPI:

```shell
python -m pip install pyisogen
```

IsoGen requires Python 3.9 or newer. The native library is loaded through
`ctypes` and does not depend on a particular CPython minor-version ABI.

Published platform wheels include a native library built from the bundled C
sources for 64-bit Windows, Linux, and macOS (Intel and Apple Silicon). Linux
and macOS wheels bundle the required FFTW 3 runtime during wheel repair, so a
compiler, CMake, and a separate FFTW installation are not needed when a
compatible wheel is available.

If pip cannot find a compatible wheel, it falls back to the source
distribution and automatically builds the native library with CMake. A source
build requires a C/C++ compiler, CMake 3.22.1 or newer, and the FFTW 3
development libraries. Installation stops with a native-build error when
those prerequisites are unavailable.

## Usage

From Python:

```python
import isogen

protein = isogen.isodist("ACDEFGHIK", type="PEPTIDE", isolen=64)
protein_brain = isogen.isodist(
    "ACDEFGHIK", type="PEPTIDE", isolen=64, method="BRAIN"
)
rna = isogen.isodist("AUGCAGUACGUA", type="RNA", isolen=64)
dna = isogen.isodist("ATGCAGTACGTA", type="DNA", isolen=64)
glucose_mass_dist = isogen.isodist("C6H12O6", type="ATOM", isolen=32)
```

The output is a numpy array of shape `(isolen, 2)` with the first column containing the monoisotopic mass and the second column containing the relative intensity. The `isolen` parameter controls the number of isotopic peaks returned.

IsoGen provides FFT, BRAIN, and neural-network methods for peptides and RNA.
The default is the exact `FFT` calculation. `BRAIN` selects the polynomial
recurrence calculation, while `NN` uses the neural-network model to predict
the distribution from a peptide or RNA sequence or neutral mass.

The `PEPTIDE` model is trained on peptide sequences, while the `RNA` model is trained on RNA sequences. The `DNA` type uses the RNA model, and

The public `ATOM` type uses the FFT method; no neural-network formula model is
available.

### Custom neural-network models

Use `isodist_custom` to generate a distribution from a binary model file rather
than one of IsoGen's bundled neural-network models:

```python
from pathlib import Path

import isogen

model_file = Path("models/my_peptide_model_64.bin")
custom = isogen.isodist_custom(
    "ACDEFGHIK",
    model_file=model_file,
    isolen=64,
    type="PEPTIDE",
)
```

The function accepts peptide, RNA, and DNA sequences or numeric neutral masses.
It always uses the neural-network method. The model must have the correct input
size for the selected input and type, and its output size must equal `isolen`.
Peptide sequence models have 20 inputs, RNA/DNA sequence models have 4 inputs,
and neutral-mass models have 5 inputs. Invalid, unreadable, or incompatible
model files raise `ValueError`. As with `isodist`, the result has shape
`(isolen, 2)`, containing neutral masses and relative intensities.

#### Training custom models

Install the training dependencies before importing the training modules:

```shell
python -m pip install -e ".[training]"
```

Training data is stored in NumPy `.npz` archives. Sequence models expect a
`seqs` array and mass models expect a `masses` array. Every archive also needs
a `dists` array with shape `(number_of_examples, isolen)`. Each row of `dists`
is the target relative-intensity distribution for its corresponding sequence
or neutral mass. For example:

```python
import numpy as np

np.savez_compressed(
    "peptide_training.npz",
    seqs=np.asarray(["ACDE", "PEPTIDE", "MARTY"]),
    dists=np.asarray(peptide_target_distributions, dtype=np.float32),
)

np.savez_compressed(
    "mass_training.npz",
    masses=np.asarray([1_000.0, 5_000.0, 10_000.0]),
    dists=np.asarray(mass_target_distributions, dtype=np.float32),
)
```

Use the engine matching the kind of input the model will receive. The helper
below directs generated models to a separate directory instead of overwriting
the models installed with IsoGen:

```python
from pathlib import Path

from isogen.isogenmass import IsoGenMassEngine
from isogen.isogenpep import IsoGenPepEngine
from isogen.isogenrna import IsoGenRNAEngine
from isogen.isogenrna_averagine import IsoGenRNAveragineEngine


model_dir = Path("trained_models")
model_dir.mkdir(exist_ok=True)


def set_model_directory(engine):
    """Set the output directory before a model is initialized or loaded."""
    engine.model.working_dir = str(model_dir)
    for model in engine.models:
        model.working_dir = str(model_dir)


# Peptide sequences: 20-element amino-acid composition input.
pep = IsoGenPepEngine(isolen=64)
set_model_directory(pep)
pep.train("peptide_training.npz", epochs=20, forcenew=True)

# RNA sequences: 4-element A/C/G/U composition input. This model is also
# used for DNA inference after IsoGen converts thymine to uracil.
rna = IsoGenRNAEngine(isolen=64)
set_model_directory(rna)
rna.train("rna_training.npz", epochs=20, forcenew=True)

# Peptide-like neutral masses: 5-element mass encoding.
mass = IsoGenMassEngine(isolen=64)
set_model_directory(mass)
mass.train_multiple(
    ["mass_training.npz"],
    inputname="masses",
    epochs=20,
    forcenew=True,
)

# RNA-like neutral masses: 5-element mass encoding.
rna_mass = IsoGenRNAveragineEngine(isolen=64)
set_model_directory(rna_mass)
rna_mass.train_multiple(
    ["rna_mass_training.npz"],
    inputname="masses",
    epochs=20,
    forcenew=True,
)
```

`IsoGenPepEngine` supports output lengths 16, 64, and 128;
`IsoGenRNAEngine` supports 64 and 128; `IsoGenMassEngine` models intended for
`isodist_custom` support 8, 32, 64, and 128; and
`IsoGenRNAveragineEngine` supports 32, 64, and 128. The output length used to
construct the engine must match the width of `dists` and the `isolen` passed to
`isodist_custom`.

After training, each engine saves a PyTorch `.pth` checkpoint and a raw `.bin`
model in `trained_models`. The `.pth` file is used to resume Python training;
pass the `.bin` file to `isodist_custom`. The generated filenames are
`isogenpep_model_<isolen>.bin`, `isogenrna_model_<isolen>.bin`,
`isogenmass_model_<isolen>.bin`, and
`isogen_rnaveragine_model<isolen>.bin`, respectively:

```python
custom = isogen.isodist_custom(
    "ACDEFGHIK",
    model_file=model_dir / "isogenpep_model_64.bin",
    isolen=64,
    type="PEPTIDE",
)
```

Passing `forcenew=True` starts from newly initialized weights. Use
`forcenew=False` to resume from a matching `.pth` checkpoint in the configured
model directory. `IsoGenMassEngine.train(...)` and
`IsoGenRNAveragineEngine.train(...)` can also generate standard FFT targets
from random masses when a custom target archive is not needed.

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


## CHANGELOG

### 1.0.5

Updated citation info.

Reformatted build to pull version from _version.py file.

### 1.0.4

Moved native-library packaging to scikit-build-core so Windows, Linux, and
macOS wheels compile the bundled C sources automatically with CMake.

Added native macOS wheels for Intel and Apple Silicon, including bundled FFTW
runtime dependencies.

Expanded Python compatibility from Python 3.13-only to Python 3.9 and newer.

Added automatic native compilation when pip falls back to the source
distribution, with clearer errors for missing build prerequisites or runtime
libraries.

### 1.0.3

Added the BRAIN polynomial-recurrence isotope calculation for peptide, RNA,
  and DNA sequence and neutral-mass inputs. 

Added `method="BRAIN"` to the Python API and command-line interface. 

Dramatically improved BRAIN performance by about double using some computational tricks the AI found.

Added side-by-side FFT, NN, and BRAIN protein-sequence example plots. 

Added runtime-dispatched AVX2/FMA neural-network acceleration on supported
  x86 processors, with a portable scalar fallback. 

Improved native normalization and large-input regression coverage.

Added a timing test script for internal use.

### 1.0.2

Added support for custom models with isogen_custom function and new C bindings for custom models.

### 1.0.1 

Small updates to README.md

### 1.0.0

Initial release. Rewrote significantly from UniDec build using AI tool to improve the release and add in atomic formula support.
