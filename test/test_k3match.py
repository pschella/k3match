import numpy as np
from k3match import spherical

N_catalog = 1e6
N_search = 1e6
ds = np.pi / (180 * 60)

theta_a = np.pi * np.random.rand(N_catalog)
phi_a = 2 * np.pi * np.random.rand(N_catalog)

theta_b = np.pi * np.random.rand(N_search)
phi_b = 2 * np.pi * np.random.rand(N_search)

idx, d = spherical(theta_a, phi_a, theta_b, phi_b, ds)

for i in range(idx.shape[0]):
    print idx[i][0], idx[i][1], d[i]

