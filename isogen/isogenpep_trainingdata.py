"""Generate peptide datasets for training and assessing IsoGen models.

This development module parses peptide sequences from tabular or MGF data,
creates randomized sequence collections, and calculates composition vectors,
masses, and native FFT isotope distributions. It is not part of IsoGen's
top-level public API.
"""

import pandas as pd
import fileinput
import multiprocessing
import matplotlib as mpl
import random
import os

if __package__:
    from .isogen_tools import *
    from .isogenwrapper import fft_gen_seq_isodist
else:
    from isogen_tools import *
    from isogenwrapper import fft_gen_seq_isodist


atoms_to_ignore = ["Br", "I", "Cl", "F", "Hg", "Mo", "Se", "B", "Cu", "Si", "As", ]


def get_big_seq(seqs):
    """Flatten peptide sequences into one amino-acid character list.

    Args:
        seqs: Iterable of peptide or protein sequence strings.

    Returns:
        Amino-acid characters from every input sequence in order.
    """
    big_seq = [char for seq in seqs for char in seq]
    return big_seq


def seq_to_dist_vecs_seqs(seq):
    """Calculate the training values associated with one peptide sequence.

    Args:
        seq: Peptide sequence accepted by the native sequence FFT engine.

    Returns:
        A tuple containing the 128-point FFT intensity distribution,
        amino-acid-count vector, original sequence, and monoisotopic mass.
    """
    dist = fft_gen_seq_isodist(seq, type="PEPTIDE", isolen=128)
    vec = peptide_to_vector(seq)
    mass = peptide_to_mass(seq)

    return dist, vec, seq, mass


def parse_tsv_file(fname):
    """Read unique peptide sequences from a tab-separated file.

    Args:
        fname: Path to a TSV file containing a ``Sequence`` column.

    Returns:
        A NumPy array of unique sequence values.
    """
    df = pd.read_csv(fname, sep="\t")
    seqs = df["Sequence"].to_numpy()
    seqs = np.unique(seqs)
    return seqs


def parse_mgf_file(fname, maxn=1000000):
    """Extract unique peptide sequences from MGF ``SEQ=`` records.

    Numeric characters and the punctuation ``+``, ``-``, and ``.`` are
    removed from each extracted value.

    Args:
        fname: Path to an MGF text file.
        maxn: Maximum number of sequence records to read.

    Returns:
        A NumPy array of unique cleaned peptide sequences.
    """
    seqs = []
    count = 0
    print("Parsing:", fname)
    for line in fileinput.input(fname):
        if "SEQ=" in line:
            seq = line.split("SEQ=")[1].split("\n")[0]
            # Drop numbers and + or - values
            seq = ''.join([i for i in seq if not i.isdigit()])
            seq = seq.replace("+", "")
            seq = seq.replace("-", "")
            seq = seq.replace(".", "")

            seqs.append(seq)

            count += 1
        if count >= maxn:
            print("Max Reached:", maxn)
            break
    return np.unique(seqs)


def seqs_to_vectors(seqs):
    """Calculate FFT distributions, vectors, and masses in parallel.

    Args:
        seqs: Iterable of peptide sequence strings.

    Returns:
        A tuple containing accepted sequences, isotope distributions,
        amino-acid-count vectors, and masses. Four ``None`` values are returned
        when no sequence produces a valid distribution and vector.
    """
    goodseqs = []
    dists = []
    vecs = []
    masses = []

    with multiprocessing.Pool(processes=8) as pool:
        results = pool.map(seq_to_dist_vecs_seqs, seqs)

    for result in results:
        dist, vec, seq, mass = result
        if dist is not None and vec is not None:
            dists.append(dist)
            vecs.append(vec)
            goodseqs.append(seq)
            masses.append(mass)


    dists = np.array(dists)
    vecs = np.array(vecs)

    goodseqs = np.array(goodseqs)
    if len(goodseqs) < 1:
        return None, None, None, None
    return goodseqs, dists, vecs, masses


def parse_file(fname, maxn=1000000, maxlen=200):
    """Convert peptide records from a TSV or MGF file into an NPZ dataset.

    The output archive is written in the current directory using the input
    file's basename and contains ``dists``, ``vecs``, ``seqs``, and ``masses``.

    Args:
        fname: Input ``.tsv`` or ``.mgf`` path.
        maxn: Maximum number of unique sequences to process, or ``None`` for
            no count limit.
        maxlen: Maximum accepted peptide length.

    Returns:
        ``None``. Results are written to a compressed NumPy archive.

    Raises:
        Exception: If the input extension is not ``.tsv`` or ``.mgf``.
    """
    if fname.endswith(".tsv"):
        seqs = parse_tsv_file(fname)

    elif fname.endswith(".mgf"):
        seqs = parse_mgf_file(fname, maxn=maxn)
    else:
        raise Exception("Unknown File Type")
    if maxn is not None:
        seqs = seqs[:maxn]

    seqs = [s for s in seqs if len(s) <= maxlen]
    print("Retreived Sequences:", len(seqs))
    goodseqs, dists, vecs, masses = seqs_to_vectors(seqs)

    #Get the filename without the extension
    filename = os.path.splitext(os.path.basename(fname))[0]
    np.savez_compressed(filename + ".npz", dists=dists, vecs=vecs, seqs=goodseqs, masses=masses)
    #np.savez_compressed("peptidedists_" + str(len(dists)) + ".npz", dists=dists, vecs=vecs, seqs=goodseqs)


def gen_random_prots(all_seqs, organism, iterations=10):
    """Generate shuffled protein assessment datasets.

    Amino acids from the input proteins are pooled and repartitioned using the
    original sequence lengths for each iteration. The combined results are
    saved as ``assessment_random_prots_<organism>.npz``.

    Args:
        all_seqs: Iterable of source protein sequences.
        organism: Label embedded in the output filename.
        iterations: Number of shuffle-and-repartition passes.

    Returns:
        ``None``.
    """
    #First concatenate all_seqs into a single list
    all_aas = []
    lengths = []
    for seq in all_seqs:
        all_aas.extend(list(seq))
        lengths.append(len(seq))

    goodseqs = np.array([])
    dists = np.array([])
    vectors = np.array([])

    for i in range(iterations):
        #Generate random protein sequences by shuffling all_aas and lengths
        # and then concatenating them
        current_all_aas = all_aas.copy()
        random_seqs = []
        random.shuffle(all_aas)
        for length in lengths:
            if length > len(all_aas):
                continue
            seq = ''.join(all_aas[:length])
            random_seqs.append(seq)
            #now remove the sequence from all_aas
            current_all_aas = current_all_aas[length:]

        curr_goodseqs, curr_dists, curr_vectors = seqs_to_vectors(random_seqs)
        if curr_goodseqs is None:
            print("No good sequences found in random sequences")
            continue

        curr_goodseqs = np.array(curr_goodseqs)
        curr_dists = np.array(curr_dists)
        curr_vectors = np.array(curr_vectors)
        if i == 0:
            goodseqs = curr_goodseqs
            dists = curr_dists
            vectors = curr_vectors
        else:
            #Add the random sequences, dists, and vectors to the goodseqs, dists, and vectors lists
            goodseqs = np.concatenate((goodseqs, curr_goodseqs))
            dists = np.concatenate((dists, curr_dists))
            vectors = np.concatenate((vectors, curr_vectors))
        print("Completed iteration", i + 1)

    #Now save the random sequences, dists, and vectors to a file
    goodseqs = np.array(goodseqs)
    dists = np.array(dists)
    vectors = np.array(vectors)
    np.savez_compressed("assessment_random_prots_" + str(organism) + ".npz", dists=dists, vecs=vectors, seqs=goodseqs)
    return


def mod_to_chemicalformula(mod):
    """Convert a modification composition mapping to a formula string.

    Args:
        mod: Object exposing a ``composition`` mapping of element names to
            integer counts.

    Returns:
        A concatenated elemental formula string.
    """
    return ''.join(f"{key}{value}" for key, value in mod.composition.items())


def gen_random_seqs_even_length(all_seqs, organism, n=10, min_length=1, max_length=200):
    """Generate and save a length-balanced random protein dataset.

    For each integer length in the inclusive range, the function creates
    ``n`` sequences from a shuffled source amino-acid pool, calculates their
    FFT distributions, vectors, and masses in a multiprocessing pool, and
    writes a compressed NumPy archive in the current directory.

    Args:
        all_seqs: Source peptide or protein sequences used as the residue pool.
        organism: Label embedded in the output filename.
        n: Number of randomized sequences generated at each length.
        min_length: Inclusive minimum generated sequence length.
        max_length: Inclusive maximum generated sequence length.

    Returns:
        ``None``. Data are written to an
        ``assessment_random_<organism>_proteins_<n>_min_<min>_max_<max>.npz``
        file.
    """
    big_seq = get_big_seq(all_seqs)
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

    np.savez_compressed("assessment_random_"+ str(organism)+ "_proteins_"+str(n)+ "_min_" + str(min_length) + "_max_" +
                        str(max_length) + ".npz",
                        dists=dists, vecs=vectors, seqs=good_seqs,masses=masses)

if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    os.chdir(r"C:\Users\Admin\Documents\martylab\Protein\IntactProtein\OtherTesting")

    min_length = 5
    max_length = 1000

    human_prot_data = np.load("human_protein_seqs.npz")
    gen_random_seqs_even_length(human_prot_data["seqs"], "human", min_length=min_length, max_length=max_length)

    yeast_prot_data = np.load("yeast_protein_seqs.npz")
    gen_random_seqs_even_length(yeast_prot_data["seqs"], "yeast", min_length=min_length, max_length=max_length)

    ecoli_prot_data = np.load("ecoli_protein_seqs.npz")
    gen_random_seqs_even_length(ecoli_prot_data["seqs"], "ecoli", min_length=min_length, max_length=max_length)

    mouse_prot_data = np.load("mouse_protein_seqs.npz")
    gen_random_seqs_even_length(mouse_prot_data["seqs"], "mouse", min_length=min_length, max_length=max_length)











