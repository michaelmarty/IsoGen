# Mass calculations

IsoGen exposes average and monoisotopic neutral-mass functions for proteins,
RNA, DNA, and elemental formulas.

## Elemental formulas

Calculate a formula's light-isotope mass or matching mass axis:

```python
formula_mass = isogen.calc_atom_monoisotopic_mass("C6H12O6")
formula_axis = isogen.calc_atom_mass_axis("C6H12O6", isolen=8)
```

The formula parser supports the same 109 elements as the native ATOM FFT
calculation. The count vector is available from
`isogen.isogenwrapper.atom_formula_to_vector`.

## Proteins

An intact protein uses `H2O` terminal composition by default:

```python
average = isogen.calc_pep_mass("PEPTIDE")
monoisotopic = isogen.calc_pep_monoisotopic_mass("PEPTIDE")
```

The supported `ion_type` values describe neutral terminal compositions:

| Ion type | Sequence to supply | Terminal shift from residue sum |
| --- | --- | --- |
| `H2O` | Full protein | +H2O |
| `a` | N-terminal fragment | -CO |
| `b` | N-terminal fragment | none |
| `c` | N-terminal fragment | +NH3 |
| `x` | C-terminal fragment | +CO2 |
| `y` | C-terminal fragment | +H2O |
| `z` | C-terminal fragment | +H2O-NH3 |

```python
b6_mass = isogen.calc_pep_monoisotopic_mass("PEPTID", ion_type="b")
y6_mass = isogen.calc_pep_monoisotopic_mass("EPTIDE", ion_type="y")
```

These are neutral masses. Charge and proton/adduct masses are not applied.
When used through `isodist`, `ion_type` changes the mass-axis origin. The
sequence intensity calculation retains its standard intact-sequence terminal
composition.

## RNA and DNA termini

Nucleic-acid calculations default to a 3'-OH and a 5'-monophosphate:

```python
rna_mass = isogen.calc_rna_monoisotopic_mass("AUGC")
dna_mass = isogen.calc_dna_monoisotopic_mass("ATGC")
```

`threeend` accepts `OH` or an unrecognized/no-adjustment value. `fiveend`
accepts:

| Value | Meaning |
| --- | --- |
| `OH` | 5'-hydroxyl |
| `MP` | 5'-monophosphate (default) |
| `TP` | 5'-triphosphate |

```python
triphosphate = isogen.calc_rna_monoisotopic_mass(
    "AUGC",
    threeend="OH",
    fiveend="TP",
)
```

RNA `a`, `b`, `c`, `d`, `w`, `x`, `y`, and `z` fragment-series names are not
accepted by `ion_type`; that keyword applies only to peptides. Pass a manually
truncated RNA sequence and select from the terminal chemistries above, or pass
a known neutral fragment mass as numeric RNA input. Numeric input uses the RNA
averagine intensity model.

When `threeend` or `fiveend` is passed through `isodist`, it adjusts the mass
axis but does not alter the sequence-model intensity vector's terminal
composition.

## Standalone mass axes

Create an axis from a known first mass:

```python
axis = isogen.calc_mass_axis(1000.0, isolen=8, isotope_spacing=1.0033)
```

Or calculate the origin from a sequence:

```python
peptide_axis = isogen.calc_pep_mass_axis("PEPTIDE", isolen=8)
rna_axis = isogen.calc_rna_mass_axis("AUGC", isolen=8)
dna_axis = isogen.calc_dna_mass_axis("ATGC", isolen=8)
formula_axis = isogen.calc_atom_mass_axis("C6H12O6", isolen=8)
```

`gen_mass_axis` dispatches between numeric, sequence, and formula input in the
same way as `isodist`.
