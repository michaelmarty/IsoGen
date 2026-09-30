"""Reproducible Windows benchmark and terminal/NN input experiments.

Run from the repository root: python -m benchmarks.proforma_performance
Times are medians of five warmed batches; no model training is performed.
"""

import json
from pathlib import Path
import platform
from statistics import median
from time import perf_counter

import numpy as np
from pyteomics import mass

import isogen
from isogen.protein_composition import resolve_proforma_composition, _lookup_composition
from isogen.isogenwrapper import _composition_isodist
from benchmarks.proforma_reference import independently_modified, pyteomics_envelope


def timing(call, repeats=50):
    call()
    samples = []
    for _ in range(5):
        start = perf_counter()
        for _ in range(repeats):
            call()
        samples.append((perf_counter() - start) * 1e6 / repeats)
    return median(samples)


def elemental_spectra(length=2048):
    spectra = {}
    for element in "CHNOS":
        isotopes = {n: abundance for n, (_, abundance) in mass.nist_mass[element].items()
                    if n and abundance > 0}
        light = min(isotopes)
        polynomial = np.zeros(length)
        for number, abundance in isotopes.items():
            polynomial[number - light] = abundance
        spectra[element] = np.fft.rfft(polynomial)
    return spectra


def spectrum(counts, spectra):
    result = np.ones_like(spectra["C"])
    for element, count in counts.items():
        if count:
            result *= spectra[element] ** count
    return result


def terminal_experiment():
    """Compare algorithms in NumPy under identical isotope data/FFT length."""
    size = 2048
    spectra = elemental_spectra(size)
    intact = dict(mass.Composition(sequence="ACDEFGHIKLMNPQRSTVWY" * 10))
    ions = ("a", "a+1", "b", "c", "x", "x+1", "y", "z", "z+1")
    deltas = [dict(mass.std_ion_comp[ion]) for ion in ions]
    final_counts = [{e: intact.get(e, 0) + d.get(e, 0) for e in "CHNOS"} for d in deltas]
    factors = [spectrum(delta, spectra) for delta in deltas]
    def direct():
        return [np.fft.irfft(spectrum(counts, spectra), n=size) for counts in final_counts]
    def shared():
        base = spectrum(intact, spectra)
        return [np.fft.irfft(base * factor, n=size) for factor in factors]
    maximum_error = max(float(np.max(np.abs(a - b))) for a, b in zip(direct(), shared()))
    return {"ions": len(ions), "fft_length": size,
            "direct_counts_numpy_us": timing(direct),
            "shared_intact_spectrum_numpy_us": timing(shared),
            "shared_spectrum_max_probability_error": maximum_error,
            "minimum_element_spectrum_magnitude": min(float(np.min(np.abs(v))) for v in spectra.values())}


def hydrogen_experiment():
    spectra = elemental_spectra(2048)
    results = []
    for copies in (1, 10, 100):
        counts = dict(mass.Composition(sequence="ACDEFGHIKLMNPQRSTVWY" * copies))
        complete = np.fft.irfft(spectrum(counts, spectra), n=2048)
        reduced = np.fft.irfft(spectrum({e: n for e, n in counts.items() if e != "H"}, spectra), n=2048)
        complete /= complete.max()
        reduced /= reduced.max()
        results.append({"residues": 20 * copies, "hydrogens": counts["H"],
                        "max_basepeak_normalized_error_omitting_H": float(np.max(np.abs(complete - reduced)))})
    return results


def main():
    output = {"platform": platform.platform(), "python": platform.python_version(),
              "units": "microseconds per call; median of five warmed batches", "proteins": []}
    pattern = "ACDEFGHIKLMNPQRSTVWY"
    for length, isolen in ((100, 128), (1000, 256), (10000, 1024)):
        plain = (pattern * (length // len(pattern)))
        annotated = "<[Carbamidomethyl]@C><[Oxidation]@M><13C>" + plain
        composition = resolve_proforma_composition(annotated)
        envelope = _composition_isodist(composition, isolen, 0)
        if not np.isfinite(envelope).all() or envelope.max() == 0:
            raise RuntimeError("Benchmark envelope invalid")
        output["proteins"].append({"residues": length, "output_length": isolen,
            "resolve_us": timing(lambda: resolve_proforma_composition(annotated), 20),
            "native_fft_and_marshalling_us": timing(lambda: _composition_isodist(composition, isolen, 0), 20),
            "end_to_end_us": timing(lambda: isogen.isodist(annotated, isolen=isolen), 20)})
    output["lookup_cached_us"] = timing(lambda: _lookup_composition("UNIMOD:35", "M"), 200)
    def uncached_lookup():
        _lookup_composition.cache_clear()
        return _lookup_composition("UNIMOD:35", "M")
    output["lookup_uncached_us_database_already_loaded"] = timing(uncached_lookup, 100)
    output["fragment_batch_100_residues_by_us"] = timing(
        lambda: isogen.calc_pep_fragment_isodists(pattern * 5), 20)
    output["terminal_options"] = terminal_experiment()
    output["nn_hydrogen_ablation"] = hydrogen_experiment()
    # Independent mass/envelope errors for the six complex test cases.
    from tests.test_pyteomics_proforma_reference import CASES
    output["pyteomics_reference"] = []
    for proforma, plain, deltas, replacements, labels in CASES:
        natural, fixed = independently_modified(plain, deltas, replacements, labels)
        origin, expected = pyteomics_envelope(natural, fixed)
        observed = isogen.isodist(proforma)
        output["pyteomics_reference"].append({"proforma": proforma,
            "mass_origin_error_Da": float(observed[0, 0] - origin),
            "max_basepeak_normalized_intensity_error": float(np.max(np.abs(observed[:, 1] - expected)))})
    destination = Path(__file__).with_name("proforma_results.json")
    destination.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
