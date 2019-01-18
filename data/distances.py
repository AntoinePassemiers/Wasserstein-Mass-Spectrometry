import os
import re
import sys
import pickle
import numpy as np
import subprocess
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt


def get_mass(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[:13] == 'CH$EXACT_MASS':
                return float(line.split(':')[1].strip())
    assert(0)


if not os.path.isfile('distances'):
    n_molecules = 619
    massdiffs = np.zeros((n_molecules, n_molecules), dtype=np.float)
    distances = np.zeros((n_molecules, n_molecules), dtype=np.float)
    for i in range(n_molecules):
        filepath1 = 'dataset1/%i.txt' % (i + 1)
        mass1 = get_mass(filepath1)
        for j in range(i):
            filepath2 = 'dataset1/%i.txt' % (j + 1)
            mass2 = get_mass(filepath2)
            output = subprocess.check_output(['../wassms', filepath1, filepath2])
            distances[i, j] = distances[j, i] = float(re.findall(r"[-+]?\d*\.\d+|\d+", str(output))[0])
            massdiffs[i, j] = massdiffs[j, i] = abs(mass1 - mass2)
        print(i)
            
    with open('distances', 'wb') as f:
        pickle.dump((distances, massdiffs), f)
else:
    n_molecules = 619
    with open('distances', 'rb') as f:
        distances, massdiffs = pickle.load(f)
    masses = np.asarray([get_mass('dataset1/%i.txt' % (i + 1)) for i in range(n_molecules)])

    X_embedded = TSNE(n_components=2, metric='precomputed').fit_transform(distances)
    plt.scatter(X_embedded[:, 0], X_embedded[:, 1], c=masses, cmap='plasma')
    plt.show()

    """
    indices = np.triu_indices(massdiffs.shape[0])
    massdiffs, distances = massdiffs[indices], distances[indices]
    indices = np.arange(massdiffs.shape[0])
    np.random.shuffle(indices)
    indices = indices[:1000]
    plt.scatter(massdiffs[indices], distances[indices])
    plt.xlabel('Absolute difference in mass')
    plt.ylabel('First Wasserstein distance')
    plt.title('Measuring the distance between molecules')
    plt.show()
    """
