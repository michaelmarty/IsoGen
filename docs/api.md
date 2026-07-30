# API reference

All functions below are available directly from `import isogen`.

## Distribution

### `isodist`

```python
isodist(input, type="PEPTIDE", isolen=128, method="FFT", **mass_kwargs)
```

Generate a two-column `(mass, relative intensity)` array. `input` may be a
numeric neutral mass or sequence. `mass_kwargs` are forwarded to
`gen_mass_axis`, including `isotope_spacing`, peptide `ion_type`, and nucleic
acid `threeend` and `fiveend`.

## Protein masses

### `calc_pep_mass`

```python
calc_pep_mass(
    sequence,
    allow_float=True,
    remove_nan=True,
    all_cyst_ox=False,
    pyroglu=False,
    round_to=2,
    ion_type="H2O",
)
```

Return the average neutral mass. Numeric input is returned as an existing mass
when `allow_float` is true.

### `calc_pep_monoisotopic_mass`

```python
calc_pep_monoisotopic_mass(
    sequence,
    allow_float=True,
    remove_nan=True,
    all_cyst_ox=False,
    pyroglu=False,
    ion_type="H2O",
)
```

Return the monoisotopic neutral mass.

### `calc_pep_mass_axis`

```python
calc_pep_mass_axis(
    sequence,
    isolen=128,
    isotope_spacing=1.0033,
    ion_type="H2O",
    **mass_kwargs,
)
```

Return a peptide mass axis beginning at the monoisotopic mass.

## Nucleic-acid masses

### RNA

```python
calc_rna_mass(sequence, threeend="OH", fiveend="MP")
calc_rna_monoisotopic_mass(sequence, threeend="OH", fiveend="MP")
calc_rna_mass_axis(
    sequence,
    isolen=128,
    isotope_spacing=1.0027,
    threeend="OH",
    fiveend="MP",
)
```

### DNA

```python
calc_dna_mass(sequence, threeend="OH", fiveend="MP")
calc_dna_monoisotopic_mass(sequence, threeend="OH", fiveend="MP")
calc_dna_mass_axis(
    sequence,
    isolen=128,
    isotope_spacing=1.0027,
    threeend="OH",
    fiveend="MP",
)
```

## Generic mass-axis functions

### `calc_mass_axis`

```python
calc_mass_axis(monoisotopic_mass, isolen=128, isotope_spacing=1.0033)
```

Return an evenly spaced one-dimensional NumPy array.

### `gen_mass_axis`

```python
gen_mass_axis(
    input,
    type="PEPTIDE",
    isolen=128,
    isotope_spacing=None,
    **mass_kwargs,
)
```

Dispatch numeric mass or sequence input to the appropriate axis calculator.
