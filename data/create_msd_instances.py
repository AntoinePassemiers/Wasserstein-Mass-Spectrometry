import os
import shutil
import IsoSpecPy
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


class Spectrum:

    def __init__(self, formula, charge=1., cutoff=.001, resolution=.01):
        self.formulas = { formula }
        self.charge = charge
        self.cutoff = cutoff
        self.data = dict()
        masses, logprobs, _ = IsoSpecPy.IsoSpecPy.IsoSpec.IsoFromFormula(
                formula=formula, cutoff=cutoff, method='threshold_relative').getConfs()
        for mass, logprob in zip(masses, logprobs):
            mz = mass / charge
            self.data[mz] = float(np.exp(logprob))
        self.set_resolution(resolution)

    def with_charge(self, charge):
        result = deepcopy(self)
        result.set_charge(charge)
        return result

    def set_charge(self, charge):
        data = dict()
        for key in self.data.keys():
            mz = (self.charge * key) / charge
            data[mz] = self.data[key]
        self.charge = charge

    def with_resolution(self, resolution):
        result = deepcopy(self)
        result.set_resolution(resolution)
        return result

    def set_resolution(self, resolution):
        if resolution:
            data = dict()
            for key in self.data.keys():
                rounded_key = resolution * int(np.round(key / resolution))
                if rounded_key in data.keys():
                    data[rounded_key] += self.data[key]
                else:
                    data[rounded_key] = self.data[key]
            self.data = data
        self.resolution = resolution

    def keys(self):
        return np.asarray(list(sorted(self.data.keys())))

    def values(self):
        return np.asarray([self.data[mz] for mz in sorted(self.data.keys())])

    def __iadd__(self, other):
        if type(other) in [int, float]:
            for key in self.data.keys():
                self.data[key] += other
        else:
            self.formulas.union(other.formulas)
            for key in other.data.keys():
                if not (key in self.data.keys()):
                    self.data[key] = other.data[key]
                else:
                    self.data[key] += other.data[key]
        return self

    def __add__(self, other):
        return deepcopy(self).__iadd__(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __imul__(self, other):
        assert(type(other) in [int, float])
        for key in self.data.keys():
            self.data[key] *= other
        return self

    def __mul__(self, other):
        return deepcopy(self).__imul__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def name(self):
        if len(self.formulas) == 1:
            return list(self.formulas)[0]
        else:
            return 'Mixture of ' + ', '.join(self.formulas)

    def save(self, filepath):
        max_intensity = float(max(list(self.data.values())))
        with open(filepath, 'w') as f:
            f.write('CH$FORMULA: %s\n' % self.name())
            f.write('PK$NUM_PEAK: %i\n' % len(self.data))
            f.write('PK$PEAK: m/z int. rel.int.\n')

            for key in sorted(self.data.keys()):
                rel = int(np.round(999. * (self.data[key] / max_intensity)))
                if rel > 0:
                    f.write('  %f %i %i\n' % (key, rel, rel))
            f.write('//\n')

    @staticmethod
    def random(target_mass, **kwargs):
        counts = [0, 0, 0, 0, 0]
        symbols = ['C', 'H', 'O', 'N', 'S']
        masses = [12, 1, 16, 14, 32]
        weights = [0.28, 0.28, 0.14, 0.03, 0.27]
        mass = 0
        while mass < target_mass:
            idx = np.random.choice(np.arange(5), p=weights)
            if mass + masses[idx] <= target_mass:
                counts[idx] += 1
                mass += masses[idx]
            else:
                pass # TODO

        formula = ''
        for idx in range(5):
            if counts[idx] > 0:
                formula += (symbols[idx] + str(counts[idx]))
        return Spectrum(formula, **kwargs)


def create_msd_folder(root, folder_name, mus, nu):
    folder_path = os.path.join(root, folder_name)
    if os.path.isdir(folder_path):
        shutil.rmtree(folder_path)
    os.makedirs(folder_path)

    # Create txt file containing the list of theoretical spectra filenames
    with open(os.path.join(folder_path, 'molecule_list.txt'), 'w') as f:
        for mu in mus:
            f.write('%s.txt\n' % mu.name())

    # Save empirical spectrum
    nu.save(os.path.join(folder_path, 'mixture.txt'))

    # Save theoretical spectra
    os.makedirs(os.path.join(folder_path, 'molecules'))
    for mu in mus:
        mu.save(os.path.join(folder_path, 'molecules', '%s.txt' % mu.name()))



if __name__ == '__main__':


    mus = list()
    mus.append(Spectrum.random(30000))
    mus.append(Spectrum.random(30000))

    plt.plot(mus[0].keys(), mus[0].values(), label=list(mus[0].formulas)[0])
    plt.plot(mus[1].keys(), mus[1].values(), label=list(mus[1].formulas)[0])
    alpha = .8
    nu = alpha * mus[0] + (1. - alpha) * mus[1]
    plt.plot(nu.keys(), nu.values(), label='mixture')
    plt.legend()
    plt.show()

    create_msd_folder('msd', 'test1', mus, nu)

    """
    mus = list()
    mus.append(Spectrum('C1482H2328N408O444S12', cutoff=.1))
    mus.append(Spectrum('C1482H2329N408O444S12', cutoff=.1))
    mus.append(Spectrum('C1482H2330N408O444S12', cutoff=.1))
    mus.append(Spectrum('C1481H2341N408O444S12', cutoff=.1))
    nu = .3 * mus[0] + .5 * mus[1] + .1 * mus[2] + .1 * mus[3]
    create_msd_folder('msd', 'test3', mus, nu)
    """