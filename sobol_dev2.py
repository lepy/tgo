import numpy
import sys

m =  2
n = int(1e3)
skip = 0
bounds = [(0, 6),
          (0, 6)]


def sobol_points(N, D):
    """ sobol.cc translated to Python 3 by Carl Sandrock 2016-03-31

    The original program is available and described at
    http://web.maths.unsw.edu.au/~fkuo/sobol/

    """
    with open('new-joe-kuo-6.21201') as f:
        unsigned = "uint64"
        # swallow header
        buffer = next(f)

        L = int(numpy.log(N)//numpy.log(2.0)) + 1

        C = numpy.ones(N, dtype=unsigned)
        for i in range(1, N):
            value = i
            while value & 1:
                value >>= 1
                C[i] += 1

        points = numpy.zeros((N, D), dtype='double')

        # XXX: This appears not to set the first element of V
        V = numpy.empty(L+1, dtype=unsigned)
        for i in range(1, L+1):
            V[i] = 1 << (32 - i)

        X = numpy.empty(N, dtype=unsigned)
        X[0] = 0
        for i in range(1, N):
            X[i] = X[i-1] ^ V[C[i-1]]
            points[i, 0] = X[i]/2**32

        for j in range(1, D):
            F_int = [int(item) for item in next(f).strip().split()]
            (d, s, a), m = F_int[:3], [0] + F_int[3:]

            if L <= s:
                for i in range(1, L+1): V[i] = m[i] << (32 - i)
            else:
                for i in range(1, s+1): V[i] = m[i] << (32 - i)
                for i in range(s+1, L+1):
                    V[i] = V[i-s] ^ (V[i-s] >> numpy.array(s, dtype=unsigned))
                    for k in range(1, s):
                        V[i] ^= numpy.array((((a >> (s-1-k)) & 1) * V[i-k]), dtype=unsigned)

            X[0] = 0
            for i in range(1, N):
                X[i] = X[i-1] ^ V[C[i-1]]
                points[i, j] = X[i]/2**32  # *** the actual points

        return points

# Distribute over bounds
B = sobol_points(n, m)

#C = numpy.column_stack([B[i] for i in range(m)])
C = B
for i in range(len(bounds)):
    C[:, i] = (C[:, i] *
               (bounds[i][1] - bounds[i][0])
               + bounds[i][0] )

from matplotlib import pyplot as plot

plot.figure(1)
plot.plot(C[:, 0], C[:, 1], 'x')
plot.show()