"""Train and run neural-network isotope models driven by neutral mass.

The mass models approximate peptide averagine isotope distributions generated
by :mod:`isogen.isogenwrapper`'s native FFT engine. Separate pretrained models
cover several mass ranges and output lengths. This is a development and model
training module rather than part of IsoGen's top-level public API.
"""

import numpy as np

if __package__:
    from .isogen_base import IsoGenEngineBase, IsoGenModelBase
    from .isogen_tools import mass_to_vector
    from .isogenwrapper import fft_gen_isodist
else:
    from isogen_base import IsoGenEngineBase, IsoGenModelBase
    from isogen_tools import mass_to_vector
    from isogenwrapper import fft_gen_isodist


class IsoGenMassEngine(IsoGenEngineBase):
    """Manage mass-input neural networks for multiple isotope-vector lengths.

    The engine maintains models producing 8, 32, 64, 128, and 1024 isotope
    values. During automatic prediction, the neutral mass determines which
    model and output length are selected.
    """

    def __init__(self, isolen=128):
        """Initialize the requested model and the complete model collection.

        Args:
            isolen: Output length for the primary ``model`` attribute. Values
                above 128 use the high-mass neural-network architecture.
        """
        super().__init__()
        self.veclen = 5
        if isolen > 128:
            modelid = 1
        else:
            modelid = 0
        self.model = IsoGenModelBase(isolen=isolen, savename="isogenmass_model_", vectorlen=self.veclen, modelid=modelid)
        self.isolen = isolen

        self.massranges = np.array(
            [[10, 1200], [400, 12000], [5000, 60000], [8000, 120000], [80000, 1100000]]
        )
        self.lengths = np.array([8, 32, 64, 128, 1024])

        self.models = []
        for l in self.lengths:
            if l > 128:
                modelid = 1
            else:
                modelid = 0
            model = IsoGenModelBase(isolen=l, savename="isogenmass_model_", vectorlen=self.veclen, modelid=modelid)
            self.models.append(model)

    def train(self, n=100000, epochs=10, length=128, forcenew=False):
        """Train the model associated with one isotope-vector length.

        Training and validation targets are generated from random masses in
        the configured range for ``length`` using the peptide FFT averagine
        implementation.

        Args:
            n: Number of training masses to generate.
            epochs: Number of training epochs.
            length: Desired isotope-distribution output length.
            forcenew: Start from newly initialized weights instead of loading
                an existing model.

        Returns:
            ``None``. If ``length`` is unsupported, a message is printed and
            the method returns without training.
        """
        # find length in self.lengths
        try:
            modelindex = np.argwhere(self.lengths == length)[0][0]
        except:
            print("Length not found in list")
            return
        massrange = self.massranges[modelindex]
        model = self.models[modelindex]
        print("Generating Training Data for isolen:" + str(length))
        input, target = gen_training_data(n, isolen=length, massrange=massrange)
        testinput, testtarget = gen_training_data(
            int(n * 0.1), isolen=length, massrange=massrange
        )
        print("Created Data:", length, n)
        trd, ted = self.create_data_loaders([input, target], [testinput, testtarget])
        model.run_training(trd, ted, epochs=epochs, forcenew=forcenew)

    def transfer_train(self, trainfile, epochs=10, length=128):
        """Continue training from distributions stored in an NPZ archive.

        Args:
            trainfile: Path to an NPZ file containing ``dists`` and ``masses``
                arrays.
            epochs: Number of additional training epochs.
            length: Output length of the model to update.

        Returns:
            ``None``.

        Raises:
            ValueError: If no model exists for the requested output length.
        """
        tdata = np.load(trainfile)
        dists = tdata["dists"]
        masses = tdata["masses"]
        trd, ted = self.setup_data(dists, masses)

        indices = np.where(self.lengths == length)
        if len(indices) > 0:
            model = self.models[indices[0][0]]
            model.run_training(trd, ted, epochs=epochs)
        else:
            print("No model with the requested isolen (" + str(length) + ") exists")
            raise ValueError("No model with selected isolen exists.")



    def train_all(self, n=100000, epochs=10, forcenew=False, ignore=[]):
        """Train every configured output-length model except ignored values.

        Args:
            n: Number of training masses generated for each model.
            epochs: Number of epochs used for each model.
            forcenew: Request newly initialized weights for model training.
            ignore: Iterable of isotope-vector lengths to skip.

        Returns:
            ``None``.
        """
        for i, l in enumerate(self.lengths):
            if l in ignore:
                continue
            self.train(n, epochs, l)

    def get_model_index(self, mass):
        """Select the first configured model range containing a mass.

        Args:
            mass: Neutral peptide mass in daltons.

        Returns:
            Integer index into ``massranges``, ``lengths``, and ``models``.

        Raises:
            ValueError: If the mass is outside every supported range.
        """
        for i in range(len(self.massranges)):
            if self.massranges[i][0] <= mass <= self.massranges[i][1]:
                return i

        print("Mass out of range", mass)
        raise ValueError("Mass out of range. Must be less than 1 MDa.")

    def inputs_to_vectors(self, inputs):
        """Encode neutral masses as neural-network input vectors.

        Args:
            inputs: Iterable of neutral masses in daltons.

        Returns:
            A two-dimensional NumPy array with one encoded vector per mass.
        """
        return np.array([mass_to_vector(m) for m in inputs])

    def predict(self, mass, isolen=None):
        """Predict a relative isotope-intensity vector for a neutral mass.

        Args:
            mass: Neutral peptide mass in daltons.
            isolen: Optional explicit output length. When omitted, the mass
                selects the model automatically.

        Returns:
            The selected neural network's predicted isotope-intensity vector.

        Raises:
            ValueError: If the mass is outside the automatic selection ranges
                or no model exists for an explicitly requested output length.
        """
        vec = mass_to_vector(mass)
        if isolen is None:
            index = self.get_model_index(mass)
            model = self.models[index]
            return model.predict(vec)
        else:
            indices = np.where(self.lengths == isolen)
            if len(indices) > 0:
                model = self.models[indices[0][0]]
                return model.predict(vec)
            else:
                print("No model with the requested isolen (" + str(isolen) + ") exists")
                raise ValueError("No model with selected isolen exists.")



def gen_training_data(n, isolen=128, massrange=[100, 1000000], log=False):
    """Generate random masses and their peptide FFT training targets.

    Args:
        n: Number of mass/distribution pairs to generate.
        isolen: Number of isotope intensities in each target.
        massrange: Two-element inclusive sampling range in daltons.
        log: Sample uniformly in log10 mass when true; otherwise sample
            uniformly in linear mass.

    Returns:
        A tuple containing a one-dimensional mass array and a two-dimensional
        NumPy array of normalized FFT isotope-intensity targets.
    """
    if log:
        rlog = np.random.uniform(np.log10(massrange[0]), np.log10(massrange[1]), n)
        randmasses = np.power(10, rlog)
    else:
        randmasses = np.random.uniform(massrange[0], massrange[1], n)
    dists = [
        fft_gen_isodist(m, type="PEPTIDE", isolen=isolen)
        for m in randmasses
    ]
    return randmasses, np.array(dists)

if __name__ == "__main__":

    eng = IsoGenMassEngine()
    #os.chdir("Z:\\Group Share\\JGP\\PeptideTraining")
    #masses, dists = gen_training_data(1000000, isolen=1024, massrange=[80000, 1100000])
    #np.savez_compressed("massdata_"+str(len(masses))+".npz", masses=masses, dists=dists)

    if False:
        n = 60000
        # eng.train_all(n, 20, forcenew=False, ignore=[])
        eng.train(n, 10, 1024, forcenew=False)
        # eng.train_multiple(["massdata_100000.npz"], epochs=40, forcenew=False, inputname="masses")
    # exit()

    if True:
        eng.train_all(n=1000000, epochs=10, forcenew=False)

