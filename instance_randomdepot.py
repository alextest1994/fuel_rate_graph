import sys

import math

import os

import json

import random
import matplotlib.pylab as plt
import seaborn as sns
import matplotlib
from itertools import permutations

import gurobipy as gp
import numpy as np
import pandas as pd
from mpmath.libmp.backend import xrange
from numpy import sin, cos
from scipy.spatial import distance_matrix


from gurobipy import GRB
from scipy.spatial.distance import pdist

import numpy as np

rnd = np.random
np.random.seed(3)

rows = 10  # number of points
columns = 2  # number of dimensions - 2=2D, 3=3D etc.
samples = np.empty((rows, columns))
for i in xrange(0, rows):
    for j in xrange(0, columns):
        samples[i][j] = rnd.randint(0, 100)
n  = 5
m = 10

from numpy.random import multinomial

b = rnd.multinomial(m, np.ones(n)/n)
a = rnd.dirichlet(np.ones(n))*m

print(a)
print(b)