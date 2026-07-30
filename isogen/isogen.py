import numpy as np

if __package__:
    from . import isogenwrapper as wrapper
    from . import mass
else:
    import isogenwrapper as wrapper
    import mass


def isodist(input, type="PEPTIDE", isolen=128, method="FFT", **mass_kwargs):
    """Generate a mass/intensity isotope distribution.

    Args:
        input: Numeric neutral mass or a protein/RNA/DNA sequence.
        type: Input type: ``PEPTIDE``, ``RNA``, or ``DNA``. DNA uses the
            RNA intensity model with a DNA-specific mass axis.
        isolen: Number of isotope values to return.
        method: Distribution engine, either ``FFT`` or ``NN``.
        **mass_kwargs: Options forwarded to :func:`mass.gen_mass_axis`, such
            as ``ion_type``, ``isotope_spacing``, ``threeend``, or ``fiveend``.

    Returns:
        A ``(isolen, 2)`` NumPy array containing neutral masses in column zero
        and relative intensities in column one.
    """
    int_dist = wrapper.gen_isodist(input, type=type, isolen=isolen, method=method)
    mass_axis = mass.gen_mass_axis(input, type=type, isolen=isolen, **mass_kwargs)
    output = np.transpose(np.vstack((mass_axis, int_dist)))
    return output

if __name__ == "__main__":
    print("Launching IsoGen")
    test_mass = 10000
    test_pep_seq = "ACDEFGHIKLMNPQRSTVWY"
    test_rna_seq = "AUGCAGUACGUA"

    print(isodist(10000, type="PEPTIDE", isolen=128, method="FFT"))
