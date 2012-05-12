#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np

setup(name='k3match',
      version='1.1',
      description='K3Match',
      author='Pim Schellart',
      author_email='P.Schellart@astro.ru.nl',
      url='',
      ext_modules=[Extension('k3match',
                   sources=['k3match/3dtree.c', 'k3match/brute.c',  'k3match/median.c', 'k3match/point.c', 'python/k3match.c'],
                   include_dirs=['.', 'k3match', np.get_include()],
                   libraries=['m'])]
     )

