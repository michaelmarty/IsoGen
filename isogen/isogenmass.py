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

    def train(self, masses=none, dists=None, n=100000, epochs=10, length=128, forcenew=False):
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
        if masses is None or dists is None:
            input, target = gen_training_data(n, isolen=length, massrange=massrange)
            testinput, testtarget = gen_training_data(int(n * 0.1), isolen=length, massrange=massrange)
        else:
            if dists.shape[1] < length:
                print("Warning: Isolen does not match training data. Changing to", dists.shape[1])
                self.isolen = dists.shape[1]

            elif dists.shape[1] > length:
                print("Warning: Isolen does not match training data. Truncating to", self.isolen)
                dists = dists[:, :self.isolen]

            # Shuffle once and split
            n = len(masses)
            indices = np.random.permutation(n)
            split = int(n * 0.9)

            train_idx = indices[:split]
            test_idx = indices[split:]

            input = masses[train_idx]
            target = dists[train_idx]
            testinput = masses[test_idx]
            testtarget = dists[test_idx]
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

    if False:
        eng.train_all(n=1000000, epochs=10, forcenew=False)

    if True:
        topdir = r"C:\Users\Admin\Documents\martylab\Protein\IntactProtein\Training"
        os.chdir(topdir)

        import fnmatch

        def match_files(directory, string, exclude=None):
            files = []
            for file in os.listdir(directory):
                if fnmatch.fnmatch(file, string):
                    if exclude is None or exclude not in file:
                        files.append(file)
            return np.array(files)

        def correct_dist_lengths(dists, length):
            corrected_dists = []
            for d in dists:
                if len(d) == length:
                    corrected_dists.append(d)
                if len(d) > length:
                    corrected = [d[i] for i in range(length)]
                    corrected_dists.append(corrected)
                if len(d) < length:
                    diff = length - len(d)
                    corrected = np.append(d, np.zeros(diff))
                    corrected_dists.append(corrected)

            return corrected_dists

        polyXdata = np.load(r"C:\Users\Admin\Documents\martylab\Protein\IntactProtein\Training\poly_X_peptides_min_6_max_25.npz")
        polyXdists = polyXdata['dists']
        polyXmasses = polyXdata['masses']

        smallpepdata = np.load(r"C:\Users\Admin\Documents\martylab\Protein\IntactProtein\Training\all_peptides_min_1_max_5.npz")
        smallpepdists = smallpepdata['dists']
        smallpepmasses = smallpepdata['masses']

        trainfiles = match_files(topdir, "*1100.npz")

        dists = []
        masses = []
        for file in trainfiles:
            data = np.load(file)
            dists.extend(data['dists'])
            masses.extend(data['masses'])

        data_indices = [[] for massrange in eng.massranges]

        for i in range(len(masses)):
            for j in range(len(eng.massranges)):
                if eng.massranges[j][0] <= masses[i] <= eng.massranges[j][1]:
                    data_indices[j].append(i)

        for i in range(len(eng.lengths)):
            curr_masses = np.array([masses[j] for j in data_indices[i]])
            # curr_masses = np.append(curr_masses, smallpepmasses)
            # curr_masses = np.append(curr_masses, polyXmasses)

            curr_dists = np.array([dists[j] for j in data_indices[i]])
            # curr_dists = np.concatenate((curr_dists, smallpepdists))
            # curr_dists = np.concatenate((curr_dists, polyXdists))

            curr_dists = np.array(correct_dist_lengths(curr_dists, eng.lengths[i]))

            eng.isolen = eng.lengths[i]
            eng.train(masses=curr_masses, dists=curr_dists, length=eng.lengths[i], epochs=10, forcenew=True)


