import numpy as np
import os

if __package__:
    from .isogen_tools import aas
    from .isogenpep_trainingdata import seqs_to_vectors
else:
    from isogen_tools import aas
    from isogenpep_trainingdata import seqs_to_vectors

# Create all possible combinations for peptides of lengths 1 to 5
def create_peptides():
    peptides = []
    for a in aas:
        peptides.append(a)

    for a in aas:
        for b in aas:
            peptides.append(a + b)

    for a in aas:
        for b in aas:
            for c in aas:
                peptides.append(a + b + c)

    for a in aas:
        for b in aas:
            for c in aas:
                for d in aas:
                    peptides.append(a + b + c + d)

    for a in aas:
        for b in aas:
            for c in aas:
                for d in aas:
                    for e in aas:
                        peptides.append(a + b + c + d + e)

    for a in aas:
        peptides.append(a + a + a + a + a + a)
        peptides.append(a + a + a + a + a + a + a)
        peptides.append(a + a + a + a + a + a + a + a)
        peptides.append(a + a + a + a + a + a + a + a + a)
        peptides.append(a + a + a + a + a + a + a + a + a + a)
        peptides.append(a + a + a + a + a + a + a + a + a + a + a)
        peptides.append(a + a + a + a + a + a + a + a + a + a + a + a)
        peptides.append(a + a + a + a + a + a + a + a + a + a + a + a + a)
        peptides.append(a + a + a + a + a + a + a + a + a + a + a + a + a + a)
        peptides.append(a + a + a + a + a + a + a + a + a + a + a + a + a + a + a)
        peptides.append(a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a)
        peptides.append(
            a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a
        )
        peptides.append(
            a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a
        )
        peptides.append(
            a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a
        )
        peptides.append(
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

    return peptides

if __name__ == "__main__":
    peptides = create_peptides()
    print("Number of peptides:", len(peptides))
    goodseqs, dists, vecs = seqs_to_vectors(peptides)
    print("Vecs Created")

    np.savez_compressed(
        "peptidedists_synthetic_" + str(len(dists)) + ".npz",
        dists=dists,
        vecs=vecs,
        seqs=goodseqs,
    )
