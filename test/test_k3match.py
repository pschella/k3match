# This file is part of K3Match.
# Copyright (C) 2016 Pim Schellart <P.Schellart@astro.ru.nl>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy as np
import k3match

na = 1e6
nb = 1e6
ds = 1./60.

np.random.seed(325423)

print "ds =", ds

ra_a = np.random.rand(na) * 360
dec_a = np.random.rand(na) * 180 - 90

ra_b = np.random.rand(nb) * 360
dec_b = np.random.rand(nb) * 180 - 90

idx_a, idx_b, d = k3match.celestial(ra_a, dec_a, ra_b, dec_b, ds)

for i in range(idx_a.shape[0]):
    print idx_a[i], idx_b[i], d[i]

