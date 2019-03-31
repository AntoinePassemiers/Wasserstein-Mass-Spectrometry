import os
import pickle
import numpy as np
import scipy.stats
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt


def get_mass(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[:13] == 'CH$EXACT_MASS':
                return float(line.split(':')[1].strip())
    assert(0)


W, M1 = pickle.load(open('wasserstein1', 'rb'))
E, _ = pickle.load(open('euclidean1', 'rb'))
J, _ = pickle.load(open('jaccard1', 'rb'))

RW, M2 = pickle.load(open('wasserstein2', 'rb'))
RE, _ = pickle.load(open('euclidean2', 'rb'))
RJ, _ = pickle.load(open('jaccard2', 'rb'))
T = pickle.load(open('tanimoto2', 'rb'))

print(W.flatten())
print(E.flatten())

correlation = lambda x, y: scipy.stats.spearmanr(x.flatten(), y.flatten())[0]
print('M-W correlation for dataset 1: ', correlation(W, M1))
print('M-E correlation for dataset 1: ', correlation(E, M1))
print('M-J correlation for dataset 1: ', correlation(J, M1))
print('W-E correlation for dataset 1: ', correlation(W, E))
print('W-J correlation for dataset 1: ', correlation(W, J))
print('J-E correlation for dataset 1: ', correlation(J, E))
print('')
print('T-J correlation for dataset 2: ', correlation(T, RJ))
print('T-W correlation for dataset 2: ', correlation(T, RW))
print('J-W correlation for dataset 2: ', correlation(RJ, RW))
print('T-E correlation for dataset 2: ', correlation(T, RE))
print('W-E correlation for dataset 2: ', correlation(RW, RE))
print('J-E correlation for dataset 2: ', correlation(RJ, RE))


masses = np.asarray([get_mass('dataset1/%i.txt' % (i + 1)) for i in range(619)])

X_embedded = TSNE(n_components=2, metric='precomputed').fit_transform(W)
plt.scatter(X_embedded[:, 0], X_embedded[:, 1], c=masses, cmap='plasma')
#plt.xlabel()
plt.axis('off')
plt.colorbar()
plt.savefig('tsne.png')
plt.clf()

plt.axis('on')
plt.scatter(X_embedded[:, 0], X_embedded[:, 1], c=masses, cmap='plasma')
plt.colorbar()
plt.savefig('tsne_with_axes.png')

with open('embedding.txt', 'w') as f:
    for i in range(len(X_embedded)):
        f.write('%f, %f\n' % (X_embedded[i, 0], X_embedded[i, 1]))



"""
indices = np.triu_indices(M1.shape[0])
massdiffs, distances = M1[indices], W[indices]
indices = np.arange(massdiffs.shape[0])
np.random.shuffle(indices)
indices = indices[:1000]
plt.scatter(massdiffs[indices], distances[indices])
plt.xlabel('Absolute difference in mass')
plt.ylabel('First Wasserstein distance')
plt.title('Measuring the distance between molecules')
plt.show()
"""
