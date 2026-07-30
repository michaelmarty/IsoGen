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
the C-terminal fragment sequence.
