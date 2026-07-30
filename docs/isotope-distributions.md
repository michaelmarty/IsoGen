# Isotope distributions

The top-level `isodist` function combines an intensity vector from the native
library with a matching neutral-mass axis.

```python
isogen.isodist(
    input,
    type="PEPTIDE",
    isolen=128,
    method="FFT",
    **mass_kwargs,
)
```

## Numeric mass input

When `input` is numeric, it becomes the first point on the mass axis:

```python
distribution = isogen.isodist(25_000, isolen=32)
assert distribution[0, 0] == 25_000
```

Numeric input selects an averagine-style intensity model. It does not infer an
elemental composition from a sequence.

## Sequence input

String input selects the corresponding sequence model. The mass-axis origin is
calculated from the sequence:

```python
distribution = isogen.isodist("PEPTIDE", type="PEPTIDE")
assert distribution[0, 0] == isogen.calc_pep_monoisotopic_mass("PEPTIDE")
```

RNA accepts `T` as uracil, and DNA accepts `U` as thymine for mass
calculations.

## Elemental-formula input

Select `type="ATOM"` to calculate a direct isotope distribution from a
formula:

```python
distribution = isogen.isodist(
    "C6H12O6",
    type="ATOM",
    isolen=32,
)
```

ATOM uses the FFT method only. Its first mass is the formula's light-isotope
mass, and its intensity vector is base-peak normalized. Formulas use standard
element symbols followed by optional nonnegative integer counts, such as
`H2O`, `C6H12O6`, or `NaCl`.

## Axis length and spacing

`isolen` controls both columns, so the mass axis always matches the intensity
vector length. Override the default spacing through a mass keyword:

```python
distribution = isogen.isodist(
    "AUGC",
    type="RNA",
    isolen=16,
    isotope_spacing=1.0029,
)
```

Formula spacing is editable in the same way:

```python
formula_distribution = isogen.isodist(
    "C6H12O6",
    type="ATOM",
    isotope_spacing=1.0033,
)
```

## Fragment ions and termini

Mass options are forwarded by `isodist`. For example:

```python
b_ion = isogen.isodist("PEPTID", ion_type="b")
triphosphate_rna = isogen.isodist(
    "AUGC",
    type="RNA",
    fiveend="TP",
    threeend="OH",
)
```

For a/b/c ions, pass the N-terminal fragment sequence. For x/y/z ions, pass
the C-terminal fragment sequence. `ion_type` applies only to peptides; named
RNA fragment-ion series are not currently implemented. RNA `threeend` and
`fiveend` values describe supported terminal chemistry for an intact or
manually truncated sequence.

These keyword arguments modify the monoisotopic mass-axis origin. They do not
change the terminal composition used by the FFT or NN sequence intensity
model.
