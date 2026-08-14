import ctypes
import os
import numpy as np
import platform
from pathlib import Path


_system = platform.system()
_library_names = {
    "Windows": "isogen.dll",
    "Linux": "isogen.so",
    "Darwin": "isogen.dylib",
}

try:
    dllname = _library_names[_system]
except KeyError as error:
    raise ImportError(
        f"IsoGen does not support {_system!r} on this installation"
    ) from error

_package_dir = Path(__file__).resolve().parent
current_path = str(_package_dir)
_packaged_library_path = _package_dir / "bin" / dllname
_source_library_path = _package_dir.parent / "bin" / dllname

# The source-tree fallback keeps local development convenient. Installed
# distributions always load the library shipped inside the package.
if _packaged_library_path.is_file():
    _library_path = _packaged_library_path
elif (
    (_package_dir.parent / "pyproject.toml").is_file()
    and _source_library_path.is_file()
):
    _library_path = _source_library_path
else:
    raise ImportError(
        f"IsoGen's native library is missing: {_packaged_library_path}. "
        "Reinstall pyisogen using a compatible wheel, or install CMake, "
        "a native compiler, and FFTW before building from source."
    )

_dll_directory_handle = None
if _system == "Windows":
    _dll_directory_handle = os.add_dll_directory(str(_library_path.parent))

try:
    isogen_c_lib = ctypes.CDLL(str(_library_path))
except OSError as error:
    raise ImportError(
        f"Unable to load IsoGen's native library {_library_path}: {error}"
    ) from error

dllpath = str(_library_path)

isodist = ctypes.c_float * 64

isogen_c_lib.fft_rna_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                            ctypes.c_int]
isogen_c_lib.fft_rna_mass_to_dist.restype = ctypes.c_float

isogen_c_lib.nn_rna_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_int]
isogen_c_lib.nn_rna_mass_to_dist.restype = ctypes.c_float

isogen_c_lib.fft_pep_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                            ctypes.c_int]
isogen_c_lib.fft_pep_mass_to_dist.restype = ctypes.c_float

isogen_c_lib.nn_pep_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_int]
isogen_c_lib.nn_pep_mass_to_dist.restype = ctypes.c_float

isogen_c_lib.fft_rna_seq_to_dist.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                             ctypes.c_int]
isogen_c_lib.fft_rna_seq_to_dist.restype = ctypes.c_float

isogen_c_lib.nn_rna_seq_to_dist.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                            ctypes.c_int]
isogen_c_lib.nn_rna_seq_to_dist.restype = ctypes.c_float

isogen_c_lib.fft_pep_seq_to_dist.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                             ctypes.c_int]
isogen_c_lib.fft_pep_seq_to_dist.restype = ctypes.c_float

isogen_c_lib.nn_pep_seq_to_dist.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                            ctypes.c_int]
isogen_c_lib.nn_pep_seq_to_dist.restype = ctypes.c_float

for _brain_mass_function_name in (
    "brain_rna_mass_to_dist",
    "brain_pep_mass_to_dist",
):
    _brain_mass_function = getattr(isogen_c_lib, _brain_mass_function_name)
    _brain_mass_function.argtypes = [
        ctypes.c_float,
        ctypes.POINTER(ctypes.c_float),
        ctypes.c_int,
        ctypes.c_int,
    ]
    _brain_mass_function.restype = ctypes.c_float

for _brain_sequence_function_name in (
    "brain_rna_seq_to_dist",
    "brain_pep_seq_to_dist",
):
    _brain_sequence_function = getattr(
        isogen_c_lib, _brain_sequence_function_name
    )
    _brain_sequence_function.argtypes = [
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_float),
        ctypes.c_int,
        ctypes.c_int,
    ]
    _brain_sequence_function.restype = ctypes.c_float

isogen_c_lib.nn_rna_mass_to_dist_custom.argtypes = [
    ctypes.c_float,
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_char_p,
]
isogen_c_lib.nn_rna_mass_to_dist_custom.restype = ctypes.c_float

isogen_c_lib.nn_pep_mass_to_dist_custom.argtypes = [
    ctypes.c_float,
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_char_p,
]
isogen_c_lib.nn_pep_mass_to_dist_custom.restype = ctypes.c_float

isogen_c_lib.nn_rna_seq_to_dist_custom.argtypes = [
    ctypes.c_char_p,
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_char_p,
]
isogen_c_lib.nn_rna_seq_to_dist_custom.restype = ctypes.c_float

isogen_c_lib.nn_pep_seq_to_dist_custom.argtypes = [
    ctypes.c_char_p,
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_char_p,
]
isogen_c_lib.nn_pep_seq_to_dist_custom.restype = ctypes.c_float

ATOM_ELEMENT_COUNT = 109

_atom_formula_to_vector_c = getattr(
    isogen_c_lib, "atom_formula_to_vector", None
)
_fft_atom_formula_to_dist_c = getattr(
    isogen_c_lib, "fft_atom_formula_to_dist", None
)

if _atom_formula_to_vector_c is not None:
    _atom_formula_to_vector_c.argtypes = [
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
    ]
    _atom_formula_to_vector_c.restype = ctypes.c_int

if _fft_atom_formula_to_dist_c is not None:
    _fft_atom_formula_to_dist_c.argtypes = [
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_float),
        ctypes.c_int,
        ctypes.c_int,
    ]
    _fft_atom_formula_to_dist_c.restype = ctypes.c_float


def _require_atom_formula_function(function, name):
    """Raise an informative error for a native library without atom support."""
    if function is None:
        raise RuntimeError(
            f"{name} is unavailable in {dllpath!r}. Rebuild or reinstall "
            "the IsoGen native library with isogenatom.c included."
        )


def _encode_formula(formula):
    """Validate and encode an elemental formula for the native API."""
    if not isinstance(formula, str):
        raise TypeError("formula must be a string")
    return formula.encode("utf-8")


def _encode_model_path(model_path):
    """Encode a string or path-like model filename for the native API."""
    if model_path is None:
        return None
    try:
        return os.fsencode(os.fspath(model_path))
    except TypeError as error:
        raise TypeError("model_path must be a string or path-like object") from error


def _check_custom_model_result(result, model_path):
    """Raise when the native library could not load a custom model."""
    if result < 0:
        raise ValueError(
            f"Unable to use custom model file {os.fspath(model_path)!r}; "
            "check that it exists and matches the requested input and output sizes"
        )


def atom_formula_to_vector(formula):
    """Convert an elemental formula to an atomic-count vector.

    The returned 109-element vector is indexed by atomic number minus one.
    For example, hydrogen is at index 0, carbon at index 5, and oxygen at
    index 7. Repeated elements in a formula are combined by the native
    parser.

    Args:
        formula: Elemental formula such as ``"C6H12O6"``.

    Returns:
        An int32 NumPy array containing the count of each element.

    Raises:
        TypeError: If ``formula`` is not a string.
        ValueError: If the formula is empty, malformed, or contains an
            unsupported element.
        RuntimeError: If the loaded native library lacks atom support.
    """
    _require_atom_formula_function(
        _atom_formula_to_vector_c, "atom_formula_to_vector"
    )
    formula_bytes = _encode_formula(formula)
    atom_counts = np.zeros(ATOM_ELEMENT_COUNT, dtype=np.int32)
    ptr = atom_counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    result = _atom_formula_to_vector_c(formula_bytes, ptr)
    if result != 0:
        raise ValueError(f"Invalid elemental formula: {formula!r}")
    return atom_counts


def fft_gen_atom_isodist(formula, isolen=128, offset=0):
    """Generate FFT isotope intensities from an elemental formula.

    Args:
        formula: Elemental formula such as ``"C6H12O6"``.
        isolen: Output vector length.
        offset: Number of leading zero-intensity isotope positions.

    Returns:
        A base-peak-normalized float32 NumPy intensity vector.

    Raises:
        TypeError: If an argument has the wrong type.
        ValueError: If the formula or output geometry is invalid.
        RuntimeError: If the loaded native library lacks atom support.
    """
    _require_atom_formula_function(
        _fft_atom_formula_to_dist_c, "fft_atom_formula_to_dist"
    )
    formula_bytes = _encode_formula(formula)
    if not isinstance(isolen, (int, np.integer)):
        raise TypeError("isolen must be an integer")
    if not isinstance(offset, (int, np.integer)):
        raise TypeError("offset must be an integer")
    if isolen <= 0:
        raise ValueError("isolen must be greater than zero")
    if offset < 0 or offset >= isolen:
        raise ValueError("offset must satisfy 0 <= offset < isolen")
    isolen = int(isolen)
    offset = int(offset)

    isodist = np.zeros(isolen, dtype=np.float32)
    ptr = isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    result = _fft_atom_formula_to_dist_c(
        formula_bytes,
        ptr,
        ctypes.c_int(isolen),
        ctypes.c_int(offset),
    )
    if result < 0:
        raise ValueError(f"Invalid elemental formula: {formula!r}")
    return isodist


def fft_atom_formula_to_dist(formula, isolen=128, offset=0):
    """Alias for :func:`fft_gen_atom_isodist` using the native API name."""
    return fft_gen_atom_isodist(formula, isolen=isolen, offset=offset)


def isogen_atom(formula, isolen=128):
    """Generate formula intensities through the compatibility API name."""
    return fft_gen_atom_isodist(formula, isolen=isolen, offset=0)


def nn_gen_seq_isodist(
    sequence, type="PEPTIDE", isolen=64, offset=0, model_path=None
):
    """Generate neural-network isotope intensities from a sequence.

    DNA sequences use the RNA model after replacing thymine with uracil.

    Args:
        sequence: Protein, RNA, or DNA sequence.
        type: ``PEPTIDE``, ``RNA``, or ``DNA``.
        isolen: Output vector length.
        offset: Number of leading zero-intensity isotope positions.
        model_path: Optional string or path-like filename of a custom binary
            model. The bundled model is used when omitted.

    Returns:
        A float32 NumPy intensity vector, or ``None`` for an unknown type.
    """
    if type == "DNA":
        sequence = sequence.upper().replace("T", "U")
    sequence_bytes = sequence.encode("utf-8")
    model_path_bytes = _encode_model_path(model_path)
    isodist = np.zeros(isolen, dtype=np.float32)
    ptr = isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    if type in ("RNA", "DNA"):
        if model_path_bytes is None:
            result = isogen_c_lib.nn_rna_seq_to_dist(
                sequence_bytes, ptr, ctypes.c_int(isolen), ctypes.c_int(offset)
            )
        else:
            result = isogen_c_lib.nn_rna_seq_to_dist_custom(
                sequence_bytes,
                ptr,
                ctypes.c_int(isolen),
                ctypes.c_int(offset),
                model_path_bytes,
            )
    elif type == "PEPTIDE":
        if model_path_bytes is None:
            result = isogen_c_lib.nn_pep_seq_to_dist(
                sequence_bytes, ptr, ctypes.c_int(isolen), ctypes.c_int(offset)
            )
        else:
            result = isogen_c_lib.nn_pep_seq_to_dist_custom(
                sequence_bytes,
                ptr,
                ctypes.c_int(isolen),
                ctypes.c_int(offset),
                model_path_bytes,
            )
    else:
        print("Unknown type for NN generation:", type)
        return None

    if model_path is not None:
        _check_custom_model_result(result, model_path)

    return isodist


def fft_gen_seq_isodist(sequence, type="PEPTIDE", isolen=128, offset=0):
    """Generate FFT isotope intensities from a sequence.

    DNA sequences use the RNA model after replacing thymine with uracil.
    ``ATOM`` and ``FORMULA`` inputs are parsed as elemental formulas.

    Args:
        sequence: Protein, RNA, or DNA sequence, or an elemental formula.
        type: ``PEPTIDE``, ``RNA``, ``DNA``, ``ATOM``, or ``FORMULA``.
        isolen: Output vector length.
        offset: Number of leading zero-intensity isotope positions.

    Returns:
        A float32 NumPy intensity vector, or ``None`` for an unknown type.
    """
    type = type.upper() if isinstance(type, str) else type
    if type in ("ATOM", "FORMULA"):
        return fft_gen_atom_isodist(sequence, isolen=isolen, offset=offset)
    if type == "DNA":
        sequence = sequence.upper().replace("T", "U")
    sequence_bytes = sequence.encode("utf-8")
    isodist = np.zeros(isolen, dtype=np.float32)
    ptr = isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    if type in ("RNA", "DNA"):
        isogen_c_lib.fft_rna_seq_to_dist(sequence_bytes, ptr, ctypes.c_int(isolen), ctypes.c_int(offset))
    elif type == "PEPTIDE":
        isogen_c_lib.fft_pep_seq_to_dist(sequence_bytes, ptr, ctypes.c_int(isolen), ctypes.c_int(offset))
    else:
        print("Unknown type for FFT generation:", type)
        return None

    return isodist


def nn_gen_isodist(
    input, type="PEPTIDE", isolen=64, offset=0, model_path=None
):
    """Generate neural-network isotope intensities from a mass or sequence.

    Args:
        input: Numeric neutral mass or sequence string.
        type: ``PEPTIDE``, ``RNA``, or ``DNA``.
        isolen: Output vector length.
        offset: Number of leading zero-intensity isotope positions.
        model_path: Optional string or path-like filename of a custom binary
            model. The bundled model is used when omitted.

    Returns:
        A float32 NumPy intensity vector, or ``None`` for an unknown type.
    """
    if isinstance(input, str):
        return nn_gen_seq_isodist(
            input, type=type, isolen=isolen, offset=offset, model_path=model_path
        )

    # Create empty array
    isodist = np.zeros(isolen).astype(np.float32)
    ptr = isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    model_path_bytes = _encode_model_path(model_path)
    if type is None:
        print("Unknown type for NN generation:", type)
        return None
    if type in ("RNA", "DNA"):
        function = (
            isogen_c_lib.nn_rna_mass_to_dist
            if model_path_bytes is None
            else isogen_c_lib.nn_rna_mass_to_dist_custom
        )
    elif type == "PEPTIDE":
        function = (
            isogen_c_lib.nn_pep_mass_to_dist
            if model_path_bytes is None
            else isogen_c_lib.nn_pep_mass_to_dist_custom
        )
    else:
        print("Unknown type for NN generation:", type)
        return None

    arguments = [
        ctypes.c_float(input),
        ptr,
        ctypes.c_int(isolen),
        ctypes.c_int(offset),
    ]
    if model_path_bytes is not None:
        arguments.append(model_path_bytes)
    result = function(*arguments)
    if model_path is not None:
        _check_custom_model_result(result, model_path)

    # Convert isodist to numpy
    isodist = np.ctypeslib.as_array(isodist)
    return isodist


def fft_gen_isodist(input, type="PEPTIDE", isolen=128, offset=0):
    """Generate FFT isotope intensities from a mass or sequence.

    Args:
        input: Numeric neutral mass, sequence, or elemental formula string.
        type: ``PEPTIDE``, ``RNA``, ``DNA``, ``ATOM``, or ``FORMULA``.
        isolen: Output vector length.
        offset: Number of leading zero-intensity isotope positions.

    Returns:
        A float32 NumPy intensity vector, or ``None`` for an unknown type.
    """
    type = type.upper() if isinstance(type, str) else type
    if isinstance(input, str):
        return fft_gen_seq_isodist(input, type=type, isolen=isolen, offset=offset)

    # Create empty array
    isodist = np.zeros(isolen).astype(np.float32)
    ptr = isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    if type in ("RNA", "DNA"):
        isogen_c_lib.fft_rna_mass_to_dist(
            ctypes.c_float(input), ptr, ctypes.c_int(isolen), ctypes.c_int(offset)
        )
        isodist = np.ctypeslib.as_array(isodist)

    elif type == "PEPTIDE":
        isogen_c_lib.fft_pep_mass_to_dist(
            ctypes.c_float(input), ptr, ctypes.c_int(isolen), ctypes.c_int(offset)
        )
        isodist = np.ctypeslib.as_array(isodist)

    else:
        print("Unknown type for FFT generation:", type)
        return None

    return np.array(isodist)


def brain_gen_seq_isodist(sequence, type="PEPTIDE", isolen=128, offset=0):
    """Generate BRAIN isotope intensities from a sequence.

    DNA sequences use the RNA calculation after replacing thymine with
    uracil.

    Args:
        sequence: Protein, RNA, or DNA sequence.
        type: ``PEPTIDE``, ``RNA``, or ``DNA``.
        isolen: Output vector length.
        offset: Number of leading zero-intensity isotope positions.

    Returns:
        A float32 NumPy intensity vector, or ``None`` for an unknown type.
    """
    type = type.upper() if isinstance(type, str) else type
    if type == "DNA":
        sequence = sequence.upper().replace("T", "U")
    sequence_bytes = sequence.encode("utf-8")
    isodist = np.zeros(isolen, dtype=np.float32)
    ptr = isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    if type in ("RNA", "DNA"):
        function = isogen_c_lib.brain_rna_seq_to_dist
    elif type == "PEPTIDE":
        function = isogen_c_lib.brain_pep_seq_to_dist
    else:
        print("Unknown type for BRAIN generation:", type)
        return None

    function(
        sequence_bytes,
        ptr,
        ctypes.c_int(isolen),
        ctypes.c_int(offset),
    )
    return isodist


def brain_gen_isodist(input, type="PEPTIDE", isolen=128, offset=0):
    """Generate BRAIN isotope intensities from a mass or sequence.

    Args:
        input: Numeric neutral mass or protein, RNA, or DNA sequence.
        type: ``PEPTIDE``, ``RNA``, or ``DNA``.
        isolen: Output vector length.
        offset: Number of leading zero-intensity isotope positions.

    Returns:
        A float32 NumPy intensity vector, or ``None`` for an unknown type.
    """
    type = type.upper() if isinstance(type, str) else type
    if isinstance(input, str):
        return brain_gen_seq_isodist(
            input, type=type, isolen=isolen, offset=offset
        )

    isodist = np.zeros(isolen, dtype=np.float32)
    ptr = isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    if type in ("RNA", "DNA"):
        function = isogen_c_lib.brain_rna_mass_to_dist
    elif type == "PEPTIDE":
        function = isogen_c_lib.brain_pep_mass_to_dist
    else:
        print("Unknown type for BRAIN generation:", type)
        return None

    function(
        ctypes.c_float(input),
        ptr,
        ctypes.c_int(isolen),
        ctypes.c_int(offset),
    )
    return isodist


def brain_pep_seq_to_dist(sequence, isolen=128, offset=0):
    """Generate BRAIN intensities from a peptide sequence."""
    return brain_gen_seq_isodist(
        sequence, type="PEPTIDE", isolen=isolen, offset=offset
    )


def brain_rna_seq_to_dist(sequence, isolen=128, offset=0):
    """Generate BRAIN intensities from an RNA sequence."""
    return brain_gen_seq_isodist(
        sequence, type="RNA", isolen=isolen, offset=offset
    )


def brain_pep_mass_to_dist(mass, isolen=128, offset=0):
    """Generate BRAIN intensities from a peptide-like neutral mass."""
    return brain_gen_isodist(
        mass, type="PEPTIDE", isolen=isolen, offset=offset
    )


def brain_rna_mass_to_dist(mass, isolen=128, offset=0):
    """Generate BRAIN intensities from an RNA-like neutral mass."""
    return brain_gen_isodist(
        mass, type="RNA", isolen=isolen, offset=offset
    )

def gen_isodist(
    input,
    type="PEPTIDE",
    isolen=128,
    offset=0,
    method="FFT",
    model_path=None,
):
    """Dispatch a mass or sequence to an isotope-distribution method.

    Args:
        input: Numeric neutral mass or sequence string.
        type: ``PEPTIDE``, ``RNA``, ``DNA``, or ``ATOM``.
        isolen: Output vector length.
        offset: Number of leading zero-intensity isotope positions.
        method: ``FFT``, ``NN``, or ``BRAIN``.
        model_path: Optional custom binary model filename for the NN method.

    Returns:
        A float32 NumPy intensity vector, or ``None`` for an unknown method or
        analyte type.
    """
    type = type.upper() if isinstance(type, str) else type
    method = method.upper() if isinstance(method, str) else method
    if type in ("ATOM", "FORMULA") and method != "FFT":
        raise ValueError("ATOM inputs support only the FFT method")
    if model_path is not None and method != "NN":
        raise ValueError("model_path is supported only by the NN method")
    if method == "FFT":
        return fft_gen_isodist(input, type=type, isolen=isolen, offset=offset)
    elif method == "NN":
        return nn_gen_isodist(
            input,
            type=type,
            isolen=isolen,
            offset=offset,
            model_path=model_path,
        )
    elif method == "BRAIN":
        return brain_gen_isodist(
            input, type=type, isolen=isolen, offset=offset
        )
    else:
        print("Unknown method for generating isotope distribution:", method)
        return None


if __name__ == "__main__":
    import time
    m =100000
    n=1000
    starttime = time.perf_counter()
    for i in range(n):
        m2 = m + np.random.uniform(-1000,1000)
        d1 = nn_gen_isodist(m2, type="PEPTIDE", isolen=128)
    print("NN Time:", (time.perf_counter()-starttime)/n*1e6, "microseconds per call")

    starttime = time.perf_counter()
    for i in range(n):
        m2 = m + np.random.uniform(-1000,1000)
        d1 = fft_gen_isodist(m2, type="PEPTIDE", isolen=128)
    print("FFT Time:", (time.perf_counter()-starttime)/n*1e6, "microseconds per call")
