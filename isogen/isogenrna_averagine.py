import numpy as np

if __package__:
    from .isogen_base import IsoGenEngineBase, IsoGenModelBase
    from .isogen_tools import mass_to_vector, rna_dict_sep, rnamass_to_isolen
    from .isogenwrapper import fft_gen_isodist, fft_gen_seq_isodist
else:
    from isogen_base import IsoGenEngineBase, IsoGenModelBase
    from isogen_tools import mass_to_vector, rna_dict_sep, rnamass_to_isolen
    from isogenwrapper import fft_gen_isodist, fft_gen_seq_isodist


def process_sequence(sequence_isolen_tuple, length=None):
    """Generate an FFT RNA distribution for a ``(sequence, isolen)`` pair.

    ``length`` is retained for compatibility with existing multiprocessing
    callers; the tuple's isotope length remains authoritative.
    """
    sequence, isolen = sequence_isolen_tuple
    return fft_gen_seq_isodist(
        sequence,
        type="RNA",
        isolen=isolen,
    )

class IsoGenRNAveragineEngine(IsoGenEngineBase):
    def __init__(self, isolen=64):
        super().__init__()
        self.veclen = 5
        self.model = IsoGenModelBase(isolen=isolen, savename="isogen_rnaveragine_model", vectorlen=self.veclen, modelid=0)
        self.lengths = np.array([32, 64, 128])
        self.massranges = np.array(
            [[200, 25000], [23000, 70000], [60000, 165000]])
        self.isolen = isolen
        self.modelindex = 0
        try:
            self.modelindex = np.argwhere(self.lengths == isolen)[0][0]
        except:
            print("Length not found in list")
            return
        self.models = []
        for l in self.lengths:
            model = IsoGenModelBase(isolen=l, savename="isogen_rnaveragine_model", vectorlen=self.veclen, modelid=0)
            self.models.append(model)


    def train(self, n=100000, epochs=10, length=128, forcenew=False):
        print("Training RNA Averagine Model")
        # find length in self.lengths
        try:
            modelindex = np.argwhere(self.lengths == length)[0][0]
            model = self.models[modelindex]
            print("Generating Training Data")
            massrange = self.massranges[modelindex]
            print("Mass Range:", massrange)

            input, target = gen_training_data(n, isolen=length, massrange=massrange)
            testinput, testtarget = gen_training_data(int(n * 0.1), isolen=length, massrange=massrange)
            print("Created Data:", length, n)
            trd, ted = self.create_data_loaders([input, target], [testinput, testtarget])

            model.run_training(trd, ted, epochs=epochs, forcenew=forcenew)

        except:
            print("Length not found in list")
            return


    def train_all(self, n=100000, epochs=10, forcenew=False, ignore=[]):
        for i, l in enumerate(self.lengths):
            if l in ignore:
                continue
            self.train(n, epochs, length=l, forcenew=forcenew)
            print("Trained RNA Averagine Model for length", l)

    def inputs_to_vectors(self, inputs):
        return np.array([mass_to_vector(m) for m in inputs])


    def predict(self, mass):
        vec = mass_to_vector(mass)
        for i in range(len(vec)):
            print("Py vec", i, ":", vec[i])
        isolen = rnamass_to_isolen(mass)
        model = None
        for m in self.models:
            if m.isolen == isolen:
                model = m
                break
        if model == None:
            print("Unsupported mass.")
            return None
        else:
            return model.predict(vec)

def rna_mass_to_dist(input):
    """Generate an FFT RNA distribution for a ``(mass, isolen)`` pair."""
    return fft_gen_isodist(input[0], type="RNA", isolen=input[1])

def gen_training_data(n, isolen=128, massrange=[100, 1000000]):
    masses = np.random.uniform(massrange[0], massrange[1], n)
    inputs = [(mass, isolen) for mass in masses]
    print("Calculating Isotope Dists...")
    dists = []
    for i in range(len(masses)):
        dist = fft_gen_isodist(masses[i], type="RNA", isolen=isolen)
        dists.append(dist)
    return np.array(masses), np.array(dists)

def rna_seq_to_formula(seq):
    #Get the molecular formula of the sequence
    formula = np.zeros(5)
    for base in seq:
        formula += rna_dict_sep[base]
    return formula

def mass_to_rnaveragine_dist(mass_len_tuple):
    """Generate an RNA averagine FFT distribution for ``(mass, isolen)``."""
    mass, isolen = mass_len_tuple
    return fft_gen_isodist(mass, type="RNA", isolen=isolen)



if __name__ == "__main__":
    eng = IsoGenRNAveragineEngine()
    eng.train_all(1000000, forcenew=True)









