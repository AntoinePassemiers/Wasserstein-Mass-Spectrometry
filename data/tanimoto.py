import os
import re
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


def get_smiles(filepath):
    with open(filepath, 'r') as f:
        for line in f.readlines():
            if line[:9] == 'CH$SMILES':
                return line.split(':')[1].strip()
    assert(0)


n_molecules = 462
distances = np.zeros((n_molecules, n_molecules), dtype=np.float)
fingerprints = list()
for i in range(n_molecules):
    filepath1 = 'dataset2/%i.txt' % (i + 1)
    m = Chem.MolFromSmiles(get_smiles(filepath1))
    fingerprints.append(AllChem.GetMorganFingerprintAsBitVect(m, 2))

for i in range(n_molecules):
    fp1 = fingerprints[i]
    for j in range(i):
        fp2 = fingerprints[j]
        # Compute Tanimoto similarity
        distances[i, j] = distances[j, i] = DataStructs.FingerprintSimilarity(fp1, fp2)
    print(i)

print(distances)

with open('tanimoto2', 'wb') as f:
    pickle.dump(distances, f)