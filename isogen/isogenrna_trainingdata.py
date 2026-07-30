"""Generate RNA datasets for training and assessing IsoGen models.

The utilities in this module create exhaustive or randomized RNA sequence
collections, calculate composition vectors and masses, and obtain reference
isotope intensities from :mod:`isogen.isogenwrapper`'s native FFT engine.
They are development utilities rather than part of IsoGen's top-level API.
"""

import multiprocessing
import os
from itertools import combinations_with_replacement

if __package__:
    from .isogen_tools import *
    from .isogenwrapper import fft_gen_seq_isodist
else:
    from isogen_tools import *
    from isogenwrapper import fft_gen_seq_isodist

def get_big_seq(seqs):
    """Flatten RNA sequences into one list of nucleotide characters.

    Args:
        seqs: Iterable of RNA sequence strings.

    Returns:
        Nucleotides from every input sequence, preserving their order.
    """
    big_seq = [char for seq in seqs for char in seq]
    return big_seq

def gen_random_bigseq(a_frac, c_frac, g_frac, u_frac, len=100000):
    """Create a shuffled nucleotide pool with requested base fractions.

    Each fraction is multiplied by ``len`` and independently rounded, so the
    final pool length can differ slightly from the requested value.

    Args:
        a_frac: Fraction of adenine residues.
        c_frac: Fraction of cytosine residues.
        g_frac: Fraction of guanine residues.
        u_frac: Fraction of uracil residues.
        len: Approximate total number of residues.

    Returns:
        A shuffled list of one-letter RNA residue codes.
    """
    seq = ""
    seq += round(a_frac*len) * "A"
    seq += round(c_frac*len) * "C"
    seq += round(g_frac*len) * "G"
    seq += round(u_frac*len) * "U"

    seq = [char for char in seq]
    np.random.shuffle(seq)
    return seq

def seq_to_dist_vecs_seqs(seq):
    """Calculate the training values associated with one RNA sequence.

    Args:
        seq: RNA sequence using the A/C/G/U alphabet.

    Returns:
        A tuple containing the 128-point FFT intensity distribution,
        nucleotide-count vector, original sequence, and monoisotopic mass.
    """
    dist = fft_gen_seq_isodist(seq, type="RNA", isolen=128)
    vec = rnaseq_to_vector(seq)
    mass = rnaseq_to_mass(seq)
    return dist, vec, seq, mass

# Create all possible combinations for RNAs of lengths 1 to 5
def create_rnas():
    """Create exhaustive short RNAs and longer homopolymer examples.

    Returns:
        Every ordered A/C/G/U sequence of lengths one through five, followed
        by single-base homopolymers of lengths six through twenty.
    """
    print("Creating synthetic short RNAs...")
    rnaseqs = []
    for a in rnas:
        rnaseqs.append(a)

    for a in rnas:
        for b in rnas:
            rnaseqs.append(a + b)

    for a in rnas:
        for b in rnas:
            for c in rnas:
                rnaseqs.append(a + b + c)

    for a in rnas:
        for b in rnas:
            for c in rnas:
                for d in rnas:
                    rnaseqs.append(a + b + c + d)

    for a in rnas:
        for b in rnas:
            for c in rnas:
                for d in rnas:
                    for e in rnas:
                        rnaseqs.append(a + b + c + d + e)

    for a in rnas:
        rnaseqs.append(a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a + a + a + a + a + a + a)
        rnaseqs.append(a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a)
        rnaseqs.append(
            a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a
        )
        rnaseqs.append(
            a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a
        )
        rnaseqs.append(
            a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a
        )
        rnaseqs.append(
            a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
            + a
        )

    return rnaseqs

def create_all_rnas(min_length=1, max_length=10):
    """Create one sorted RNA sequence for each nucleotide composition.

    This uses combinations with replacement rather than permutations, so
    sequences with the same composition but different residue order are not
    generated separately.

    Args:
        min_length: Shortest sequence length to include.
        max_length: Longest sequence length to include, inclusively.

    Returns:
        A list of composition-unique RNA sequence strings.
    """
    nucleotides = ["A", "C", "G", "U"]
    seqs = []

    for i in range(min_length, max_length+1):
        for combo in combinations_with_replacement(nucleotides, i):
            seqs.append(''.join(combo))
    print("Produced", str(len(seqs)), "sequences.")
    return seqs

def create_rand_rnas(n=1000, start=6, maxlen=500):
    """Generate RNA sequences with random residues and lengths.

    Args:
        n: Number of sequences to generate.
        start: Inclusive minimum sequence length.
        maxlen: Exclusive maximum sequence length.

    Returns:
        A list of random A/C/G/U sequence strings.
    """
    seqs = []
    for i in range(n):
        if i % 10000 == 0:
            print("Creating random RNA number", i)
        length = np.random.randint(start, maxlen)
        seq = "".join(np.random.choice(["A", "C", "G", "U"], length))
        seqs.append(seq)
    return seqs

def seqs_to_vectors(seqs):
    """Convert RNA sequences into FFT distributions and count vectors.

    Sequences that raise an exception during conversion are reported and
    skipped.

    Args:
        seqs: Iterable of RNA sequence strings.

    Returns:
        NumPy arrays containing accepted sequences, isotope distributions, and
        nucleotide-count vectors. Returns ``None`` if no sequence succeeds.
    """
    goodseqs = []
    dists = []
    vecs = []

    for i, seq in enumerate(seqs):
        if i% 10000 == 0:
            print("Processing RNA number", i)
        try:
            dist = rnaseq_to_dist(seq)
            vec = rnaseq_to_vector(seq)
            goodseqs.append(seq)
            dists.append(dist)
            vecs.append(vec)
        except:
            print("Failed to process", seq)
            pass

    dists = np.array(dists)
    vecs = np.array(vecs)
    goodseqs = np.array(goodseqs)
    if len(goodseqs) < 1:
        return None
    return goodseqs, dists, vecs

def gen_random_seqs_even_length(n=100, min_length=1, max_length=200):
    """Generate and save a length-balanced random RNA assessment dataset.

    For each integer length in the inclusive range, the function creates
    ``n`` sequences, calculates their distributions, vectors, and masses in a
    multiprocessing pool, and saves the results in the current directory.

    Args:
        n: Number of random sequences generated at each length.
        min_length: Inclusive minimum sequence length.
        max_length: Inclusive maximum sequence length.

    Returns:
        ``None``. Data are written to an
        ``assessment_random_RNAs_<n>_min_<min>_max_<max>.npz`` file.
    """
    big_seq = gen_random_bigseq(0.295, 0.204, 0.205, 0.296)
    np.random.shuffle(big_seq)

    random_seqs = []
    good_seqs = []
    vectors = []
    masses = []
    dists = []


    for length in range(min_length, max_length+1):
        current = 0
        index = 0

        while current < n:
            seq = ''.join(big_seq[index:index+length])

            random_seqs.append(seq)
            current += 1
            index += length

            if index + length > len(big_seq):
                np.random.shuffle(big_seq)
                index = 0

    print("Processing sequences...")
    with multiprocessing.Pool(processes=8) as pool:
        results = pool.map(seq_to_dist_vecs_seqs, random_seqs)


    print("Parsing results...")
    for r in results:
        if r[0] is not None and r[1] is not None and r[2] is not None and r[3] is not None:
            dists.append(np.array(r[0]))
            vectors.append(r[1])
            good_seqs.append(str(r[2]))
            masses.append(r[3])


    np.savez_compressed("assessment_random_RNAs_"+str(n)+ "_min_" + str(min_length) + "_max_" +
                        str(max_length) + ".npz",
                        dists=dists, vecs=vectors, seqs=good_seqs,masses=masses)



if __name__ == "__main__":
    os.chdir(r"C:\Users\Admin\Documents\martylab\RNA_SeqData\Assessment")

    if True:
        n = 5000
        min_length = 10
        max_length = 500
        gen_random_seqs_even_length(n=n, min_length=min_length, max_length=max_length)

    if False:
        all_seqs = create_all_rnas(min_length=2, max_length=20)

        dists = []
        vectors = []
        good_seqs = []
        masses = []

        print("Processing sequences...")
        with multiprocessing.Pool(processes=8) as pool:
            results = pool.map(seq_to_dist_vecs_seqs, all_seqs)

        print("Parsing results...")
        for r in results:
            if r[0] is not None and r[1] is not None and r[2] is not None and r[3] is not None:
                dists.append(np.array(r[0]))
                vectors.append(r[1])
                good_seqs.append(str(r[2]))
                masses.append(r[3])

        np.savez_compressed("synthetic_RNAs_" + str(len(good_seqs)) + ".npz",
                            dists=dists, vecs=vectors, seqs=good_seqs, masses=masses)



