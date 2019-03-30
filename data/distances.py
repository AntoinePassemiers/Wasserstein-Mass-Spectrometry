import os
import re
import sys
import pickle
import numpy as np
import subprocess


METRIC = 'E'

if METRIC == 'E':
    out_filename, dataset, method = 'euclidean1', 'dataset1', 'E'
elif METRIC == 'W':
    out_filename, dataset, method = 'wasserstein1', 'dataset1', 'W'
elif METRIC == 'J':
    out_filename, dataset, method = 'jaccard1', 'dataset1', 'J'
elif METRIC == 'RE':
    out_filename, dataset, method = 'euclidean2', 'dataset2', 'E'
elif METRIC == 'RW':
    out_filename, dataset, method = 'wasserstein2', 'dataset2', 'W'
elif METRIC == 'J2':
    out_filename, dataset, method = 'jaccard2', 'dataset2', 'J'
n_molecules = 619 if dataset == 'dataset1' else 462


def get_mass(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[:13] == 'CH$EXACT_MASS':
                return float(line.split(':')[1].strip())
    assert(0)


massdiffs = np.zeros((n_molecules, n_molecules), dtype=np.float)
distances = np.zeros((n_molecules, n_molecules), dtype=np.float)
masses = list()
for i in range(n_molecules):
    filepath = '%s/%i.txt' % (dataset, i + 1)
    masses.append(get_mass(filepath))

for i in range(n_molecules):
    filepath1 = '%s/%i.txt' % (dataset, i + 1)
    mass1 = masses[i]
    for j in range(i):
        filepath2 = '%s/%i.txt' % (dataset, j + 1)
        mass2 = masses[j]
        output = subprocess.check_output(['../wassms', filepath1, filepath2, '--m', method])
        distance = float(re.findall(r"[-+]?\d*\.\d+|\d+", str(output))[0])
        if 'R' in METRIC:
            distance /= (mass1 * mass2)
        distances[i, j] = distances[j, i] = distance
        massdiffs[i, j] = massdiffs[j, i] = abs(mass1 - mass2)
    print(i)

with open(out_filename, 'wb') as f:
    pickle.dump((distances, massdiffs), f)
