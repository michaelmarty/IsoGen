"""Run and train neural-network isotope models for RNA sequences.

The RNA engine encodes sequences as A/C/G/U count vectors and selects a
pretrained model according to sequence length. Models produce either 64 or 128
relative isotope intensities. This training-oriented module is not part of
IsoGen's top-level public API.
"""

import numpy as np
import os

if __package__:
    from .isogen_base import IsoGenEngineBase, IsoGenModelBase
    from .isogen_tools import rnaseq_to_vector
else:
    from isogen_base import IsoGenEngineBase, IsoGenModelBase
    from isogen_tools import rnaseq_to_vector


class IsoGenRNAEngine(IsoGenEngineBase):
    """Manage RNA sequence models for supported sequence-length ranges.

    Sequences containing 1–200 residues use the 64-point model, while
    sequences containing 201–500 residues use the 128-point model.
    """

    def __init__(self, isolen=64):
        """Initialize the RNA model collection.

        Args:
            isolen: Output length assigned to the primary ``model`` attribute.
                Supported pretrained values are 64 and 128.
        """
        super().__init__()
        self.isolen = isolen
        self.seqlengthranges = np.array([[1, 200], [201, 500]])
        self.lengths = np.array([64, 128])
        self.inputname = "seqs"
        self.models = []
        for l in self.lengths:
            if l > 128:
                modelid = 1
            else:
                modelid = 0
            model = IsoGenModelBase(isolen=l, savename="isogenrna_model_", vectorlen=4, modelid=modelid)
            self.models.append(model)

        for i in range(len(self.lengths)):
            if self.lengths[i] == self.isolen:
                self.model = self.models[i]

    def inputs_to_vectors(self, inputs):
        """Encode RNA sequences as nucleotide-count vectors.

        Args:
            inputs: Iterable of RNA sequence strings.

        Returns:
            A two-dimensional NumPy array with A/C/G/U counts for each input
            sequence.
        """
        return np.array([rnaseq_to_vector(m) for m in inputs])

    def get_model_index(self, seqlength):
        """Select a model from an RNA sequence length.

        Args:
            seqlength: Number of residues in the RNA sequence.

        Returns:
            Zero for lengths 1–200 or one for lengths 201–500.

        Raises:
            ValueError: If the length is outside the supported ranges.
        """
        for i, range in enumerate(self.seqlengthranges):
            if seqlength >= range[0] and seqlength <= range[1]:
                return i

        print("Sequence length out of range.")
        raise ValueError("Sequence length out of range, must be under 500.")

    def predict(self, seq):
        """Predict a relative isotope-intensity vector for an RNA sequence.

        The sequence length automatically selects the 64- or 128-point model.

        Args:
            seq: RNA sequence using the A/C/G/U alphabet.

        Returns:
            The selected neural network's predicted isotope-intensity vector.

        Raises:
            ValueError: If the sequence length is outside 1–500 residues.
        """
        self.check(seq)
        vec = rnaseq_to_vector(seq)
        modelindex = self.get_model_index(len(seq))
        model = self.models[modelindex]
        return model.predict(vec)

    def check(self, seq):
        """Warn when an RNA sequence exceeds the trained length range.

        Args:
            seq: RNA sequence to inspect.

        Returns:
            ``None``.
        """
        if len(seq) > 500:
            print("Warning: Sequence is too long. Behavior may be unpredictable.")




if __name__ == "__main__":
    os.chdir(r"C:\Users\Admin\Documents\martylab\RNA_SeqData\Training")

    trainfile0 = "synthetic_RNAs_min2_max35.npz"
    trainfile1 = "synthetic_RNAs_min36_max220.npz"
    trainfile2 = "synthetic_RNAs_min30_max120.npz"

    trainfile3 = "synthetic_RNAs_min180_max520.npz"
    trainfile4 = "synthetic_RNAs_min200_max400.npz"

    eng1 = IsoGenRNAEngine(isolen=64)
    eng2 = IsoGenRNAEngine(isolen=128)

    #eng1.train_multiple([trainfile0, trainfile0, trainfile1, trainfile2], epochs=20, forcenew=True)
    eng2.train_multiple([trainfile3, trainfile4], epochs=20, forcenew=True)