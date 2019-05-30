import os
import subprocess
import shutil
import IsoSpecPy
import numpy as np
import scipy.stats
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

    def noisy(self):
        result = deepcopy(self)
        result.add_noise()
        return result

    def add_noise(self):
        resolution = self.resolution
        noisy_data = dict()

        # Add noise in the intensity domain
        for mz in self.data.keys():
            intensity = self.data[mz]
            if intensity > 0:
                log_intensity = np.log(intensity) + np.random.normal(0, .01)
                self.data[mz] = np.exp(log_intensity)

        # Add noise in the mass domain
        for mz in self.data.keys():
            intensity = self.data[mz]
            rounded_mz = resolution * int(np.round(mz / resolution))
            xs = np.arange(rounded_mz - 5*resolution, rounded_mz+5*resolution, resolution)
            cdf = scipy.stats.norm.cdf(xs, loc=mz, scale=.01)
            ys = (cdf[1:] - cdf[:-1])
            ys *= intensity
            noisy_mzs = (xs + 0.5*resolution)[:-1]
            for noisy_mz, y in zip(noisy_mzs, ys):
                noisy_mz = resolution * int(np.round(noisy_mz / resolution))
                if not noisy_mz in noisy_data.keys():
                    noisy_data[noisy_mz] = y
                else:
                    noisy_data[noisy_mz] += y
        self.data = noisy_data

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


def formula_to_latex(formula):
    latex = '$' + formula[0]
    for i in range(1, len(formula)):
        c_old = formula[i-1]
        c = formula[i]
        if c.isdigit():
            if not c_old.isdigit():
                latex = latex + '_{' + c
            else:
                latex = latex + c
        else:
            if c_old.isdigit():
                latex = latex + '} ' + c
            else:
                latex = latex + c
    if formula[-1].isdigit():
        latex = latex + '}'
    latex = latex + '$'
    return latex


DECONVMS_PATH = '../deconvms'


def solve_msd(mus, weights):
    k = len(mus)
    nu = sum(float(weights[i]) * mus[i] for i in range(k))
    nu = nu.noisy()
    create_msd_folder('msd', 'test3', mus, nu)
    folder = os.path.join('msd', 'test3')

    output = subprocess.check_output(
        [DECONVMS_PATH,
         os.path.join(folder, 'mixture.txt'),
         os.path.join(folder, 'molecule_list.txt'),
         os.path.join(folder, 'molecules'),
         '--niter',
         '10000'])

    predicted_weights = list()
    for line in output.decode('utf-8').split('\n'):
        line = line.replace('\r', '')
        print(line)
        if line.startswith('Weight'):
            predicted_weights.append(float(line.split(':')[1]))
    return np.sqrt(np.mean((predicted_weights - weights) ** 2.))


def test1():
    rmsd = list()
    std = list()
    for n_molecules in range(11, 18):
        s = list()
        for _ in range(2):
            mus = list()
            for _ in range(n_molecules):
                mus.append(Spectrum.random(30000))

            k = len(mus)

            alpha = np.ones(k)
            weights = np.random.dirichlet(alpha)
            print(weights)
            
            s.append(solve_msd(mus, weights))
        rmsd.append(np.mean(s))
        std.append(np.std(s))
        print(rmsd)
        print(std)
    return rmsd, std



def test3():

    rmsd = list()
    std = list()
    for charge in range(1, 31):
        s = list()
        for _ in range(10):
            mus = list()
            mus.append(Spectrum('C1482H2328N408O444S12', cutoff=.1, charge=charge))
            mus.append(Spectrum('C1482H2329N408O444S12', cutoff=.1, charge=charge))
            mus.append(Spectrum('C1482H2330N408O444S12', cutoff=.1, charge=charge))
            mus.append(Spectrum('C1481H2341N408O444S12', cutoff=.1, charge=charge))
            weights = np.asarray([.3, .5, .1, .1])
            nu = .3 * mus[0] + .5 * mus[1] + .1 * mus[2] + .1 * mus[3]
            s.append(solve_msd(mus, weights))
        rmsd.append(np.mean(s))
        std.append(np.std(s))
        print(rmsd)
        print(std)
    return rmsd, std


if __name__ == '__main__':

    test1()

