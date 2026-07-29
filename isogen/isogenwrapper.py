import ctypes
import os
import numpy as np
import platform


def find_dll(targetfile, dir):
    if dir is None:
        return ""

    for entry in os.scandir(dir):
        if entry.is_file() and entry.name == targetfile:
            # print("Found DLL within:", entry.path)
            return entry.path

        elif entry.is_dir():
            result = find_dll(targetfile, entry.path)
            if result:
                return result

    return ""

def start_at_iso(targetfile, guess=None):
    if guess is not None:
        if os.path.isdir(guess):
            result = find_dll(targetfile, guess)
            if result:
                return result

    return targetfile


current_path = os.path.dirname(os.path.realpath(__file__))
bin_path = os.path.realpath(os.path.join(current_path, "..", "bin"))
packaged_bin_path = os.path.join(current_path, "bin")

if platform.system() == "Windows":
    dllname = "isogen.dll"
elif platform.system() == "Linux":
    dllname = "isogen.so"

else:
    print("Not yet implemented for MacOS")
    dllname = "isogen.dylib"


dllpath = start_at_iso(dllname, guess=bin_path)
if dllpath == dllname:
    dllpath = start_at_iso(dllname, guess=packaged_bin_path)

if not dllpath:
    print("DLL not found anywhere:", dllname)

isodist = ctypes.c_float * 64

isogen_c_lib = ctypes.CDLL(dllpath)

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


def nn_gen_seq_isodist(sequence, type="PEPTIDE", isolen=64, offset=0):
    """Generate an isotope distribution from a peptide, RNA, or DNA sequence."""
    if type == "DNA":
        sequence = sequence.upper().replace("T", "U")
    sequence_bytes = sequence.encode("utf-8")
    isodist = np.zeros(isolen, dtype=np.float32)
    ptr = isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    if type in ("RNA", "DNA"):
        isogen_c_lib.nn_rna_seq_to_dist(sequence_bytes, ptr, ctypes.c_int(isolen), ctypes.c_int(offset))
    elif type == "PEPTIDE":
        isogen_c_lib.nn_pep_seq_to_dist(sequence_bytes, ptr, ctypes.c_int(isolen), ctypes.c_int(offset))
    else:
        print("Unknown type for NN generation:", type)
        return None

    return isodist


def fft_gen_seq_isodist(sequence, type="PEPTIDE", isolen=128, offset=0):
    """Generate an isotope distribution from a peptide, RNA, or DNA sequence."""
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


def nn_gen_isodist(input, type="PEPTIDE", isolen=64, offset=0):
    if isinstance(input, str):
        return nn_gen_seq_isodist(input, type=type, isolen=isolen, offset=offset)

    # Create empty array
    isodist = np.zeros(isolen).astype(np.float32)
    ptr = isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    res = b""
    if type is None:
        print("Unknown type for NN generation:", type)
        return None
    if type in ("RNA", "DNA"):
        res = b"Rna"
        # Call the C function
        isogen_c_lib.nn_rna_mass_to_dist(
            ctypes.c_float(input),
            ptr,
            ctypes.c_int(isolen),
            ctypes.c_int(offset)
        )
    elif type == "PEPTIDE":
        res = b"Peptide"
        # Call the C function
        isogen_c_lib.nn_pep_mass_to_dist(
            ctypes.c_float(input),
            ptr,
            ctypes.c_int(isolen),
            ctypes.c_int(offset)
        )

    # Convert isodist to numpy
    isodist = np.ctypeslib.as_array(isodist)
    return isodist


def fft_gen_isodist(input, type="PEPTIDE", isolen=128, offset=0):
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

def gen_isodist(input, type="PEPTIDE", isolen=128, offset=0, method="FFT"):
    if method == "FFT":
        return fft_gen_isodist(input, type=type, isolen=isolen, offset=offset)
    elif method == "NN":
        return nn_gen_isodist(input, type=type, isolen=isolen, offset=offset)
    else:
        print("Unknown method for generating isotope distribution:", method)
        return None


if __name__ == "__main__":
    import time
    m =100000
    n=10000
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
    # d1 = fft_gen_isodist(m, type="PEPTIDE")
    # print(d1)

    # import matplotlib.pyplot as plt
    # plt.plot(d1/np.amax(d1), label="Isogen")
    # plt.legend()
    # plt.show()
    exit()
    #
    #
    # eng = IsoGenWrapper(dllpath=dllpath)
    # n = 10000
    # random_masses = np.random.uniform(1000, 60000, n)
    # starttime = time.perf_counter()
    # for mass in random_masses:
    #     dist = eng.gen_isodist(mass)
    #
    # print("Isogen Time:", time.perf_counter()-starttime)
    # print("Microseconds Per:", (time.perf_counter()-starttime)/n* 1e6)
    #
    # starttime = time.perf_counter()
    # for mass in random_masses:
    #     dist = eng.gen_isomike(mass)
    #
    # print("Isomike Time:", time.perf_counter()-starttime)
    # print("Microseconds Per:", (time.perf_counter()-starttime)/n* 1e6)
    #
    # starttime = time.perf_counter()
    # for mass in random_masses:
    #     dist = eng.gen_isofft(mass)
    #
    # print("FFT Time:", time.perf_counter()-starttime)
    # print("Microseconds Per:", (time.perf_counter()-starttime)/n * 1e6)
