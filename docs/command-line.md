# Command line

The console command and module entry point are equivalent:

```shell
isogen --help
python -m isogen --help
```

## CSV distributions

Print a distribution to standard output:

```shell
isogen dist ACDEFGHIK --type PEPTIDE --isolen 64
```

Write it to a file:

```shell
isogen dist AUGCAGUACGUA --type RNA --output rna.csv
```

Generate a distribution from an elemental formula:

```shell
isogen dist C6H12O6 --type ATOM --isolen 32 --output glucose.csv
```

The output contains a header followed by `mass,intensity` rows. Use
`--method FFT`, `--method BRAIN`, or `--method NN` to choose the calculation
engine. `ATOM` inputs support `FFT` only.

## Example stick plots

Open the example figure:

```shell
isogen plot
```

Save without opening a window:

```shell
isogen plot --save isodist_examples.png --no-show
```

The example figure compares the same protein sequence with FFT, NN, and BRAIN.
The `--method` option selects the engine used by the protein-mass and RNA
panels; the three protein-sequence comparison panels remain fixed.
