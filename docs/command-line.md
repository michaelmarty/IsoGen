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

The output contains a header followed by `mass,intensity` rows. Use
`--method FFT` or `--method NN` to choose the calculation engine.

## Example stick plots

Open the example figure:

```shell
isogen plot
```

Save without opening a window:

```shell
isogen plot --save isodist_examples.png --no-show
```
