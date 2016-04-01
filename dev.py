import unittest
import numpy
from _tgo import *

from UQToolbox.sobol_lib import i4_sobol_generate
#from sobol import sobol_points
 # TODO: Replace with latinhypercube sampling used in differentialevolution.py
import numpy
import scipy.spatial
import scipy.optimize

from scipy.optimize._differentialevolution import DifferentialEvolutionSolver, _make_random_gen

#DES = DifferentialEvolutionSolver()

#print DifferentialEvolutionSolver.init_population_lhs()
#print _make_random_gen(None)

m =  2
n = int(1e3)
skip = 0
bounds = [(0, 6),
          (0, 6)]

B = i4_sobol_generate(m, n, skip)
C = numpy.column_stack([B[i] for i in range(m)])

# Distribute over bounds
# TODO: Find a better way to do this
for i in range(len(bounds)):
    C[:, i] = (C[:, i] *
               (bounds[i][1] - bounds[i][0])
               + bounds[i][0] )

sobol_points(N, D, f)

from matplotlib import pyplot as plot

plot.figure(1)
plot.plot(C[:, 0], C[:, 1], 'x')
#plot.show()


