# coding:utf-8

import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(currentdir, '..'))

import maskgen

seed = 123456
L = 8
nsamples = 10
for c in maskgen.percolation(seed, L, nsamples):
    print(c)
