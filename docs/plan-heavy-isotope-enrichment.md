# Plan: Heavy-isotope enrichment option for `type="ATOM"`

## Goal

Let a user describe a partially labeled sample, e.g.

```python
isogen.isodist("C6H12O6", type="ATOM", isolen=32, heavy_isotopes={"2H": 0.05})
```

`{"2H": 0.05}` means each hydrogen is 5% deuterium. The per-element isotope
distribution is replaced by a linear combination of the natural distribution
and a pure heavy-isotope delta:

```
p_new = (1 - f) * p_natural + f * delta_heavy
```

Supported labels at minimum: `2H`, `13C`, `15N`, `18O` (also `17O`, `33S`,
`34S`, and any isotope present in the native isotope table). The mixed
single-atom distribution is FFT'd and used **instead of** the precomputed FFT
for that element when forming the power/product in the convolution.

## Background (existing code)

- `src/isogenatom.c::formula_vector_to_probability_dist` multiplies, per
  element, `setup_ft(Z)^count` (precomputed FFT via `fft_from_precomputed`
  when `length == precomputed_isotope_length`, else built from
  `construct_isotope_array` + FFTW r2c).
- `fft_atom_counts_to_dist` already exists with sparse *fixed* isotopes
  (`setup_fixed_isotope_ft`); that is 100% labeling of specific atoms, not
  fractional enrichment. It is the closest precedent for API shape.
- Python: `isogen/isogen.py::isodist` -> `isogenwrapper.gen_isodist` ->
  `fft_gen_atom_isodist`. Mass axis comes from `mass.gen_mass_axis`.
- Reference algorithm: `UniDecDev\Scripts\JGP\IsoGen_Analysis\strange_isodist_generator.py`
  builds per-element vectors (`c = [1, 0.0108, 0, 0]`), zero-pads, FFTs, and
  multiplies `ft**count`; heavy vectors are deltas (`c_heavy = [0,1,0,0]`).
  Here the same idea is applied via a mixture, in the native code.

## Design

### Semantics
- Input: `dict[str, float]`, key = `"<mass_number><Symbol>"` (`"2H"`, `"13C"`,
  `"15N"`, `"18O"`); also accept `"H2"`-style? **No** (keep one format);
  also accept `"D"` as alias of `"2H"`.
- Value: fraction `f` in `[0, 1]` of that element's atoms that are the heavy
  isotope. `f == 0` is a no-op; `f == 1` is fully labeled.
- Interpretation: **fraction of the element**, not extra to natural
  abundance. The natural distribution is scaled by `(1-f)` and the heavy
  isotope delta is added with weight `f`. The heavy isotope's natural
  abundance is thus counted only in the `(1-f)` part.
- Multiple isotopes of one element (`{"18O": 0.1, "17O": 0.05}`) are allowed
  if the sum is <= 1: `p = (1 - sum f) * p_nat + sum f_i * delta_i`.
- Validation: unknown isotope, element absent from formula (warn/ignore vs
  error -> **ValueError**, to catch typos), fraction outside [0,1] or
  non-finite, sum > 1, non-dict -> `TypeError`/`ValueError`.
- Mass axis: the shift of the labeled mass is *not* a fixed offset; the
  intensity vector stays indexed from the light-isotope monoisotopic mass
  (nominal-mass steps), consistent with the present ATOM behaviour. The
  axis spacing is unchanged. Optional follow-up: report the intensity-weighted
  mean mass shift in docs.

### Native (C) changes — `src/isogenatom.c/.h`
1. New helper `build_mixed_element_ft(Z, heavy_mass_numbers[], fractions[], n, length, ftlen, outft)`:
   - Build time-domain array via `construct_isotope_array(Z, array, length)`
     (natural distribution, base index = lightest isotope, offset by
     `round(mass - reference)` as already done in `fft_atom_counts_to_dist`).
   - Scale by `(1 - sum f)`, then add `f_i` at index
     `round(mass_i - lightest_mass)` for each heavy isotope. Index must be
     `< length` (else error).
   - r2c FFT with FFTW -> `outft`. Bypasses `fft_from_precomputed` entirely.
   - Note: heavy delta must be exact for the isotope's nominal-mass
     index; natural table must be normalized (sum to 1) before mixing.
2. In `formula_vector_to_probability_dist`, take an extra optional
   "enrichment" list (`enrich_elements`, `enrich_mass_numbers`,
   `enrich_fractions`, `enrich_count`). For each element with enrichment use
   the helper instead of `setup_ft`; otherwise unchanged (precomputed path
   retained -> no perf regression for unlabeled calls).
3. New exported function:
   ```c
   float fft_atom_counts_to_dist_enriched(
       const int atom_counts[109],
       const int* fixed_elements, const int* fixed_mass_numbers,
       const int* fixed_counts, int fixed_count,
       const int* enrich_elements, const int* enrich_mass_numbers,
       const double* enrich_fractions, int enrich_count,
       float* isodist, int isolen, int offset);
   ```
   Existing `fft_atom_counts_to_dist` stays ABI-compatible (wrapper calling the
   new function with `enrich_count = 0`). Add the declaration to `isogenatom.h`.
4. Length selection: extend the `mean`/`variance` estimate used to size the
   zero-padded FFT so it uses the *mixed* distribution moments for enriched
   elements (heavy labels widen/shift the envelope, e.g. 13C at 50% on large
   molecules). Keep the 12-sigma + 32 margin.
5. Numerical care: negative/denormal FFT noise is already handled by
   `normalize_isodist`; verify no `abs` hack hides sign problems for large
   counts (check `complex_power` for large exponents with mixed FFT bins where
   `|z|` is small).

### Python changes
- `isogen/isogenwrapper.py`: ctypes signature for the new function;
  `fft_gen_atom_isodist(formula, isolen, offset, heavy_isotopes=None)`;
  parser `_parse_heavy_isotopes` (symbol -> atomic number via
  `protein_mods._ELEMENTS`, mass number validation against the native isotope
  table via existing `isotope_mass`), fallback error if the library lacks the
  new symbol (like `_require_atom_formula_function`).
- Route through `gen_isodist` and `isogen.isodist(..., heavy_isotopes=None)`;
  only valid for `type="ATOM"` (and optionally formula inputs through
  existing composition resolution). Other types: `ValueError` if supplied.
- Pass-through from `fft_atom_formula_to_dist` and `isogen_atom`
  (compat alias) as optional keyword.
- CLI (`docs\command-line.md`): `--heavy-isotopes "2H=0.05,13C=0.02"` (or JSON)
  if an ATOM CLI path exists; verify in `isogen` CLI module.
- Mass axis: add mean-shift note only; no behaviour change.

## Tests (`tests/test_atom_heavy_isotopes.py`, pytest)

Pure-Python reference implemented in the test (NumPy, float64, mirrors the
strange_isodist_generator approach): per-element mixed vector -> rfft -> `**count`
-> irfft -> normalize by max.

1. **Identity**: `heavy_isotopes=None`, `{}`, and `{"2H": 0.0}` equal the
   current output (bit-for-bit or `atol=1e-7`).
2. **Reference match** for `{"2H": 0.05}` on `C6H12O6`, `H2O`, large `C100H200`
   plus `13C`, `15N`, `18O` cases, multi-label combos, `atol≈1e-6`.
3. **Full labeling (f=1)**: equals pure delta convolution, i.e. matches the
   existing fixed-isotope path (`fft_atom_counts_to_dist` with fixed isotopes)
   shifted appropriately — cross-validates the two native paths.
4. **Closed-form checks**: single atom `H` with f=0.05 -> `[~1, p/(1-p)...]`
   ratio; `H2` binomial: M+1/M+0 and M+2/M+0 analytic; a 13C-only monocarbon
   check against `(1-f)*0.9893 : f + (1-f)*0.0107`.
5. **Mixture linearity**: for one atom, output before normalization is exactly
   `(1-f) p_nat + f δ` (verify via native return `max_probability`).
6. **Monotonicity / mean shift**: mean of the distribution increases linearly
   with `f` (`Δmean = count * f * (1 - nat_mean)`) within tolerance.
7. **Validation errors**: bad key (`"3H"` unsupported by table? `"XX"`),
   negative / >1 / NaN fractions, sum > 1, element not in formula, wrong
   types, used with non-ATOM type, `isolen` too short, heavy index beyond
   FFT length.
8. **Truncation/padding**: short `isolen` (e.g. 4) with heavy labels equals the
   head of a long run (guards against FFT wrap-around).
9. **Public API/doc test**: `tests/test_public_api.py` signature check,
   and `isodist(..., heavy_isotopes=...)` axis shape `(isolen, 2)`, `dist_only`.
10. **Regression**: existing `tests/test_atom_wrapper.py` and fixed-isotope
    tests unchanged and passing.

### External validation
- **pyteomics**: `pyteomics.mass.isotopologues`/`isotopic_composition_abundance`
  with a custom `mass.Composition` and a modified `nist_mass` table where each
  element's abundances are replaced by the mixture (`{1: (1.00782503, 1-f_total...), 2: (2.01410178, f)}`).
  The repo already has `tests/test_pyteomics_proforma_reference.py`, so
  pyteomics is an accepted optional (skip if missing) dev dependency.
- **BRAIN / other calculators** (optional): IsoSpec (`pip install IsoSpecPy`)
  supports custom abundances per element; compare envelopes for
  `C6H12O6` with `{2H: 0.05}` and `C100H160N30O30S` with `{13C: 0.02, 15N: 0.98 mix}`
  — exact tolerance `1e-5` relative on top peaks. Tests skip cleanly if the
  package is absent. Also a one-off manual comparison with enviPat/Mercury
  (documented, not automated).
- Brute-force multinomial enumeration for tiny molecules (H2O, CH4) as an
  independent third check with no external package.

## Documentation
- `docs\isotope-distributions.md`: new "Partially labeled (enriched) samples"
  section: semantics (fraction of the element, replaces natural abundance by a
  linear mixture), accepted key format, examples for `2H`, `13C`, `15N`,
  `18O`, multi-label, errors, performance note (enriched elements skip the
  precomputed FFT table).
- `docs\api.md`: new parameter on `isodist`, `fft_gen_atom_isodist`,
  `fft_atom_formula_to_dist`; C API section for the new native function.
- `README.md` feature bullet + short example; `docs\command-line.md` flag
  if CLI applies; changelog entry (if present).
- Docstrings updated (Args/Raises).

## Implementation order
1. Pure-Python reference + failing tests (identity, reference, closed form).
2. C helper + extended probability function + new export; rebuild library
   (CMake, FFTW) and confirm export.
3. ctypes/Python API + validation; make tests pass.
4. External validation tests (pyteomics, IsoSpec, brute force; skip if absent).
5. Docs, README, CLI; run full `pytest`.
6. Perf check: unlabeled calls unaffected; enriched call timing for a
   `C1000H2000...` formula (FFT length growth).

## Risks / open decisions
- **Fraction semantics** (assumed: fraction of the element replaced,
  natural scaled by `1-f`). Alternative "added on top of natural" is
  rejected as it breaks normalization; confirm with the user.
- 2H nominal-index relies on `round(mass - lightest)`; fine for H/C/N/O/S;
  elements with non-integer-spaced isotopes (e.g. Br/Cl are fine) need the
  same rounding rule already used.
- Prebuilt binaries in `bin\` and wheels need rebuild; Python wrapper must
  degrade with a clear error if running against an old library.
- Large counts: bins of the mixed FFT with tiny `|z|` may underflow;
  verify against the reference at count ~10^4.
