import numpy as np

if __package__:
    from . import isogenwrapper as wrapper
    from . import mass
else:
    import isogenwrapper as wrapper
    import mass


def _apply_charge(mass_axis, charge, polarity):
    """Convert a neutral mass axis to m/z when a positive charge is supplied."""
    if charge is None or charge < 1:
        return mass_axis

    if not isinstance(polarity, str):
        raise TypeError("polarity must be 'positive' or 'negative'")

    polarity = polarity.lower()
    if polarity not in ("positive", "negative"):
        raise ValueError("polarity must be 'positive' or 'negative'")

    proton_mass = mass.mass_proton
    if polarity == "negative":
        proton_mass *= -1
    return (mass_axis + charge * proton_mass) / charge


def isodist(
        input,
        type="PEPTIDE",
        isolen=128,
        method="FFT",
        dist_only=False,
        charge=None,
        polarity="positive",
        use_modifications=True,
        **mass_kwargs,
):
    """Generate a mass/intensity isotope distribution.

    Args:
        input: Numeric neutral mass, protein/RNA/DNA sequence, or elemental
            formula.
        type: Input type: ``PEPTIDE``, ``RNA``, ``DNA``, or ``ATOM``. DNA
            uses the RNA intensity model with a DNA-specific mass axis.
            ATOM accepts formulas and uses the FFT method.
        isolen: Number of isotope values to return.
        method: Distribution engine: ``FFT``, ``NN``, or ``BRAIN``.
        dist_only: Return only the relative-intensity vector.
        charge: Ion charge. ``None`` (the default) returns a neutral mass
            axis; values greater than or equal to one return an m/z axis.
        polarity: ``"positive"`` (default) uses ``(M + zH) / z``;
            ``"negative"`` uses ``(M - zH) / z``.
        use_modifications: Include resolvable ProForma atomic composition in
            peptide intensities (default ``True``). Unsupported annotations
            print a warning and are ignored for intensities. ``False`` uses
            the legacy unmodified envelope; known mass shifts still apply.
        **mass_kwargs: Options forwarded to :func:`mass.gen_mass_axis`, such
            as ``ion_type``, ``isotope_spacing``, ``threeend``, or ``fiveend``.

    Returns:
        A ``(isolen, 2)`` NumPy array containing neutral masses or m/z values
        in column zero and relative intensities in column one.
    """
    type = type.upper() if isinstance(type, str) else type
    method = method.upper() if isinstance(method, str) else method
    composition = _resolve_distribution_input(input, type, use_modifications, mass_kwargs)
    int_dist = wrapper.gen_isodist(input, type=type, isolen=isolen, method=method,
                                  use_modifications=use_modifications, _composition=composition)
    if dist_only:
        return int_dist

    axis_input = composition.mass if composition is not None else input
    mass_axis = mass.gen_mass_axis(axis_input, type=type, isolen=isolen, **mass_kwargs)
    mass_axis = _apply_charge(mass_axis, charge, polarity)
    output = np.transpose(np.vstack((mass_axis, int_dist)))
    return output


def isodist_custom(
        input,
        model_file,
        isolen,
        type="PEPTIDE",
        charge=None,
        polarity="positive",
        use_modifications=True,
        **mass_kwargs,
):
    """Generate a mass/intensity distribution using a custom NN model.

    Args:
        input: Numeric neutral mass or protein/RNA/DNA sequence.
        model_file: String or path-like filename of a custom binary model.
        isolen: Number of isotope values to return. The model output size must
            match this value.
        type: Input type: ``PEPTIDE``, ``RNA``, or ``DNA``. DNA uses the RNA
            intensity model with a DNA-specific mass axis.
        charge: Ion charge. ``None`` (the default) returns a neutral mass
            axis; values greater than or equal to one return an m/z axis.
        polarity: ``"positive"`` (default) uses ``(M + zH) / z``;
            ``"negative"`` uses ``(M - zH) / z``.
        use_modifications: Default ``True``. Custom residue-count models print
            a warning and ignore composition changes in intensities. ``False``
            suppresses this warning and explicitly selects the legacy envelope.
        **mass_kwargs: Options forwarded to :func:`mass.gen_mass_axis`.

    Returns:
        A ``(isolen, 2)`` NumPy array containing neutral masses or m/z values
        in column zero and relative intensities in column one.

    Raises:
        ValueError: If ``type`` is unsupported or the model cannot be loaded.
    """
    type = type.upper() if isinstance(type, str) else type
    if type not in ("PEPTIDE", "RNA", "DNA"):
        raise ValueError("Custom models support PEPTIDE, RNA, and DNA inputs")
    composition = _resolve_distribution_input(input, type, use_modifications, mass_kwargs)

    int_dist = wrapper.gen_isodist(
        input,
        type=type,
        isolen=isolen,
        method="NN",
        model_path=model_file,
        use_modifications=use_modifications,
        _composition=composition,
    )
    mass_axis = mass.gen_mass_axis(
        composition.mass if composition is not None else input,
        type=type,
        isolen=isolen,
        **mass_kwargs,
    )
    mass_axis = _apply_charge(mass_axis, charge, polarity)
    return np.transpose(np.vstack((mass_axis, int_dist)))


def _resolve_distribution_input(input, type, use_modifications, mass_kwargs):
    if not isinstance(use_modifications, (bool, np.bool_)):
        raise TypeError("use_modifications must be a boolean")
    if type == "PEPTIDE" and isinstance(input, str) and (
        wrapper.needs_proforma_parser(input) or str(mass_kwargs.get("ion_type", "H2O")).lower() != "h2o"
    ):
        if __package__:
            from .protein_composition import resolve_proforma_composition
        else:
            from protein_composition import resolve_proforma_composition
        resolved = resolve_proforma_composition(
            input, use_modifications=use_modifications,
            ion_type=mass_kwargs.get("ion_type", "H2O"))
        if mass_kwargs.get("all_cyst_ox", False):
            resolved.mass -= resolved.sequence.count("C") * mass.mass_H_monoisotopic
        if mass_kwargs.get("pyroglu", False):
            if resolved.sequence.startswith("E"):
                resolved.mass -= mass.mass_water_monoisotopic
            elif resolved.sequence.startswith("Q"):
                resolved.mass -= mass.mass_OH_monoisotopic
        return resolved
    return None


if __name__ == "__main__":
    print("Launching IsoGen")
    test_mass = 10000
    test_pep_seq = "ACDEFGHIKLMNPQRSTVWY"
    test_rna_seq = "AUGCAGUACGUA"

    print(isodist(10000, type="PEPTIDE", isolen=128, method="FFT"))
    print(isodist("PdW2CO3", type="ATOM", isolen=16, method="FFT"))
    print(isodist(test_pep_seq, type="PEPTIDE", isolen=128, method="FFT"))
    print(isodist(test_pep_seq, type="PEPTIDE", isolen=128, method="BRAIN"))
    print(isodist(test_rna_seq, type="RNA", isolen=128, method="FFT"))
