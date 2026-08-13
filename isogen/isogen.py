import numpy as np

if __package__:
    from . import isogenwrapper as wrapper
    from . import mass
else:
    import isogenwrapper as wrapper
    import mass


def isodist(input, type="PEPTIDE", isolen=128, method="FFT", dist_only=False, **mass_kwargs):
    """Generate a mass/intensity isotope distribution.

    Args:
        input: Numeric neutral mass, protein/RNA/DNA sequence, or elemental
            formula.
        type: Input type: ``PEPTIDE``, ``RNA``, ``DNA``, or ``ATOM``. DNA
            uses the RNA intensity model with a DNA-specific mass axis.
            ATOM accepts formulas and uses the FFT method.
        isolen: Number of isotope values to return.
        method: Distribution engine: ``FFT``, ``NN``, or ``BRAIN``.
        **mass_kwargs: Options forwarded to :func:`mass.gen_mass_axis`, such
            as ``ion_type``, ``isotope_spacing``, ``threeend``, or ``fiveend``.

    Returns:
        A ``(isolen, 2)`` NumPy array containing neutral masses in column zero
        and relative intensities in column one.
    """
    type = type.upper() if isinstance(type, str) else type
    method = method.upper() if isinstance(method, str) else method
    int_dist = wrapper.gen_isodist(input, type=type, isolen=isolen, method=method)
    if dist_only:
        return int_dist

    mass_axis = mass.gen_mass_axis(input, type=type, isolen=isolen, **mass_kwargs)
    output = np.transpose(np.vstack((mass_axis, int_dist)))
    return output

def isodist_custom(input, model_file, isolen, type="PEPTIDE", **mass_kwargs):
    """Generate a mass/intensity distribution using a custom NN model.

    Args:
        input: Numeric neutral mass or protein/RNA/DNA sequence.
        model_file: String or path-like filename of a custom binary model.
        isolen: Number of isotope values to return. The model output size must
            match this value.
        type: Input type: ``PEPTIDE``, ``RNA``, or ``DNA``. DNA uses the RNA
            intensity model with a DNA-specific mass axis.
        **mass_kwargs: Options forwarded to :func:`mass.gen_mass_axis`.

    Returns:
        A ``(isolen, 2)`` NumPy array containing neutral masses in column zero
        and relative intensities in column one.

    Raises:
        ValueError: If ``type`` is unsupported or the model cannot be loaded.
    """
    type = type.upper() if isinstance(type, str) else type
    if type not in ("PEPTIDE", "RNA", "DNA"):
        raise ValueError("Custom models support PEPTIDE, RNA, and DNA inputs")

    int_dist = wrapper.gen_isodist(
        input,
        type=type,
        isolen=isolen,
        method="NN",
        model_path=model_file,
    )
    mass_axis = mass.gen_mass_axis(
        input,
        type=type,
        isolen=isolen,
        **mass_kwargs,
    )
    return np.transpose(np.vstack((mass_axis, int_dist)))

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
