# API reference

All functions below are available directly from `import isogen`.

## Distribution

### `isodist`

```python
isodist(input, type="PEPTIDE", isolen=128, method="FFT",
        use_modifications=True, **mass_kwargs)
```

Generate a two-column `(mass, relative intensity)` array. `input` may be a
numeric neutral mass, biopolymer sequence, or elemental formula. Select
`method="FFT"`, `method="BRAIN"`, or `method="NN"` for biopolymer inputs.
BRAIN is an absolute polynomial-recurrence calculation. Select `type="ATOM"`
for a formula; ATOM supports the FFT method only. `mass_kwargs` are forwarded
to `gen_mass_axis`, including `isotope_spacing`, peptide `ion_type`, and
nucleic acid `threeend` and `fiveend`.

Peptide ProForma composition is included by default. Set
`use_modifications=False` to use the legacy unmodified envelope while retaining
known modification mass shifts. Unsupported or mass-only annotations print a
warning to standard error and are ignored for intensities; other resolvable
modifications on the same peptide still apply. Unknown annotations also have no
mass contribution. Malformed syntax and cross-linked/multiple peptide expressions
still raise an error.

FFT supports natural-abundance elements and fixed isotope labels. BRAIN accepts
CHNOS compositions and phosphorus, and prints a warning before using FFT for
other elements or explicit isotopes. Bundled residue-count NN models likewise
warn and use FFT for composition-changing modifications. Custom NN models retain
their requested model and warn that modifications are ignored for intensities.
`isodist_custom` accepts the same `use_modifications` flag.

```python
modified = isodist("M[Oxidation]")
legacy = isodist("M[Oxidation]", use_modifications=False)
labeled = isodist("K[UNIMOD:259]")  # six 13C and two 15N replacements
global_label = isodist("<13C>PEPTIDE")
```

The flag is also available on `gen_isodist`, the FFT/NN/BRAIN sequence and generic
wrappers, `brain_pep_seq_to_dist`, and `isogen_tools.peptide_to_dist`. For the CLI,
use `python -m isogen dist "M[Oxidation]" --ignore-modifications`.

## Elemental formulas

```python
calc_atom_monoisotopic_mass(formula)
calc_atom_mass_axis(formula, isolen=128, isotope_spacing=1.0033)
```

These functions parse formulas with the native ATOM parser. The mass axis
begins at the light-isotope mass used for isotope index zero by the FFT
calculation.

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

Both peptide mass functions automatically recognize supported ProForma
annotations. The lower-level helpers are also public:

```python
calc_proforma_mass(sequence, monoisotopic=True, ion_type="H2O")
strip_proforma(sequence)
```

`calc_proforma_mass` resolves UniMod, PSI-MOD, and RESID annotations and
returns a neutral mass. `strip_proforma` returns the unannotated sequence.

### `calc_pep_fragments`

```python
calc_pep_fragments(
    sequence,
    ion_types=None,
    monoisotopic=True,
    ambiguous_rule="reject",
    fragmentation_type=None,
)
```

Return all backbone-cleavage fragment masses in a dictionary keyed by ion
name. `ion_types` accepts any combination of `a`, `a+1`, `b`, `c`, `x`,
`x+1`, `y`, `y-1`, `z`, and `z'`. The aliases `z+1`, `z•`, `z·`, and
`z.` normalize to `z'`. `z` matches Pyteomics; `z'` is one neutral hydrogen
heavier.

When `ion_types` is omitted, `fragmentation_type` selects the conventional
series for `CID`, `HCD`, `SID`, `IRMPD`, `ETD`, `ECD`, `EThcD`, `BYCZ*`,
`UVPD`, `UVPD4`, `UVPD6`, or `UVPD9`. Explicit `ion_types` take precedence.
If both arguments are omitted, the default remains `b`/`y`.
ProForma modifications are retained when their site is part of a fragment.
By default, ions whose mass depends on an ambiguous localization or residue
are omitted. Set `ambiguous_rule="both"` to return their possible masses as a
series of numbered keys ordered by mass, such as `b2#1` and `b2#2`.

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

Dispatch numeric mass, sequence, or ATOM formula input to the appropriate axis
calculator.
