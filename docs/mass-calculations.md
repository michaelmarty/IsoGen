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

### Modified proteins (ProForma)

Protein mass functions accept the mass-bearing parts of the current
[HUPO-PSI ProForma notation](https://github.com/HUPO-PSI/ProForma). UniMod,
PSI-MOD, and RESID names and accessions are resolved from databases bundled
with IsoGen. Numeric shifts, elemental formulas, terminal modifications,
labile or unlocalized modifications, and global fixed modifications are also
supported:

```python
oxidized = isogen.calc_pep_monoisotopic_mass("EM[Oxidation]E")
psi_mod_oxidized = isogen.calc_pep_monoisotopic_mass(
    "EM[MOD:00719]E"
)
resid_oxidized = isogen.calc_pep_monoisotopic_mass(
    "EM[RESID:AA0581]E"
)
n_terminal = isogen.calc_pep_monoisotopic_mass("[Acetyl]-PEPTIDE")
fixed_cysteine = isogen.calc_pep_monoisotopic_mass(
    "<[Carbamidomethyl]@C>ACDC"
)
```

`J`, `O`, and `U` have defined residue masses. `B` and `Z` use the arithmetic
midpoint of their D/N and E/Q possibilities. Because `X` has no unique mass,
it must carry a complete signed mass gap, for example
`RTAAX[+367.0537]WT`; bare `X` raises `ValueError`.

XL-MOD, GNO/GNOme, cross-links, branched peptides, isotope replacement, and
multi-peptidoform expressions are not yet supported.

For `isodist`, ProForma annotations affect the mass-axis origin but are
removed before FFT, NN, or BRAIN calculates intensities. Thus the current
isotope envelope is the unmodified-sequence approximation. An exact modified
envelope will require applying every modification formula to an elemental
composition; a numeric mass shift alone is insufficient to determine one.

Use the dedicated helpers when parsing is useful independently:

```python
mass = isogen.calc_proforma_mass("EM[UNIMOD:35]E")
sequence = isogen.strip_proforma("EM[UNIMOD:35]E")  # "EME"
```

The supported `ion_type` values describe neutral terminal compositions:

| Ion type | Sequence to supply | Terminal shift from residue sum |
| --- | --- | --- |
| `H2O` | Full protein | +H2O |
| `a` | N-terminal fragment | -CO |
| `a+1` | N-terminal fragment | -CO+H |
| `b` | N-terminal fragment | none |
| `c` | N-terminal fragment | +NH3 |
| `x` | C-terminal fragment | +CO2 |
| `x+1` | C-terminal fragment | +CO2+H |
| `y` | C-terminal fragment | +H2O |
| `y-1` | C-terminal fragment | +H2O-H |
| `z` | C-terminal fragment | +H2O-NH3 (Pyteomics `z`) |
| `z'` | C-terminal fragment | +H2O-NH2 (one H above `z`) |

`z+1`, `z•`, `z·`, and `z.` are accepted as aliases for `z'`.

`calc_pep_fragments` can select ion series by fragmentation method:

| `fragmentation_type` | Ion series |
| --- | --- |
| `CID`, `HCD`, `SID`, `IRMPD` | `b`, `y` |
| `ETD`, `ECD` | `c`, `z'` |
| `EThcD`, `BYCZ*` | `b`, `y`, `c`, `z'` |
| `UVPD` | `a`, `b`, `c`, `x`, `y`, `z'` |
| `UVPD4` | `a`, `a+1`, `x+1`, `y-1` |
| `UVPD6` | `a`, `a+1`, `x+1`, `x`, `y-1`, `z'` |
| `UVPD9` | `a`, `a+1`, `b`, `c`, `x`, `x+1`, `y`, `y-1`, `z'` |

```python
b6_mass = isogen.calc_pep_monoisotopic_mass("PEPTID", ion_type="b")
y6_mass = isogen.calc_pep_monoisotopic_mass("EPTIDE", ion_type="y")
fragments = isogen.calc_pep_fragments("PEPTIDE")
# {"b1": ..., "b2": ..., ..., "y1": ..., "y2": ..., ...}
ecd_fragments = isogen.calc_pep_fragments("PEPTIDE", ion_types=("c", "z'"))
uvpd_fragments = isogen.calc_pep_fragments("PEPTIDE", fragmentation_type="UVPD9")

modified = isogen.calc_pep_fragments("EM[Oxidation]E")
ambiguous = isogen.calc_pep_fragments(
    "AS[#g1]T[Phospho#g1]K", ambiguous_rule="both"
)
# A cleavage separating the candidate sites produces "b2#1" and "b2#2".
```

Supplying `ion_types` explicitly overrides `fragmentation_type`.

Localized, terminal, global fixed, region-localized, labile, and unlocalized
ProForma modifications are supported. The default `ambiguous_rule="reject"`
omits any fragment with more than one possible mass. With
`ambiguous_rule="both"`, each distinct possibility is returned as a float
under a numbered key such as `b2#1`, ordered by increasing mass. Ambiguous
`B` and `Z` residues follow the same rule. `J` remains unambiguous by mass
because isoleucine and leucine have equal residue masses.

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
