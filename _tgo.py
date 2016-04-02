#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" execfile('tgo.py')
"""
from __future__ import division, print_function, absolute_import
import numpy
import scipy.spatial
import scipy.optimize


def tgo(func, bounds, args=(), g_cons=None, g_args=(), n=100,
        k_t=None, callback=None, minimizer_kwargs=None, disp=False):
    """
    Finds the global minima of a function using topograhphical global
    optimisation.

    Parameters
    ----------
    func : callable
        The objective function to be minimized.  Must be in the form
        ``f(x, *args)``, where ``x`` is the argument in the form of a 1-D array
        and ``args`` is a  tuple of any additional fixed parameters needed to
        completely specify the function.

    bounds : sequence
        Bounds for variables.  ``(min, max)`` pairs for each element in ``x``,
        defining the lower and upper bounds for the optimizing argument of
        `func`. It is required to have ``len(bounds) == len(x)``.
        ``len(bounds)`` is used to determine the number of parameters in ``x``.
        Use ``None`` for one of min or max when there is no bound in that
        direction. By default bounds are ``(None, None)``.

    args : tuple, optional
        Any additional fixed parameters needed to completely specify the
        objective function.

    g_cons : sequence of callable functions, optional
        Function used to define a limited subset to defining the feasible set
        of solutions in R^n in the form g(x) <= 0 applied as g : R^n -> R^m

        NOTE: If the ``constraints`` sequence used in the local optimization
              problem is not defined in ``minimizer_kwargs`` and a constrained
              method is used then the ``g_cons`` will be used.
              (Defining a ``constraints`` sequence in ``minimizer_kwargs``
               means that ``g_cons`` will not be added so if equality
               constraints and so forth need to be added then the inequality
               functions in g_cons need to be added again).

    g_args : sequence of tuples, optional
        Any additional fixed parameters needed to completely specify the
        feasible set functions ``g_cons``.
        ex. g_cons = (f1(x, *args1), f2(x, *args2))
        then
            g_args = (args1, args2)

    n : int, optional
        Number of sampling points used in the construction of the topography
        matrix.

    k_t : int, optional
        Defines the number of columns constructed in the k-t matrix. The higher
        k is the lower the amount of minimisers will be used for local search
        routines. If None the empirical model of Henderson et. al. (2015) will
        be used. (Note: Lower ``k_t`` values decrease performance depending
        on the local solver used, but could potentially be more robust due
        to testing more local minimisers in the function hypersuface.

    minimizer_kwargs : dict, optional
        Extra keyword arguments to be passed to the minimizer
        ``scipy.optimize.minimize`` Some important options could be:

            method : str
                The minimization method (e.g. ``SLSQP``)
            args : tuple
                Extra arguments passed to the objective function (``func``) and
                its derivatives (Jacobian, Hessian).

            options : {ftol: 1e-12}

    disp : bool, optional # (TODO)
        Display status messages

    callback : callable, `callback(xk, convergence=val)`, optional: # (TODO)
        A function to follow the progress of the minimization. ``xk`` is
        the current value of ``x0``. ``val`` represents the fractional
        value of the population convergence.  When ``val`` is greater than one
        the function halts. If callback returns `True`, then the minimization
        is halted (any polishing is still carried out).


    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a `OptimizeResult` object.
        Important attributes are:
        ``x`` the solution array corresponding to the global minimum,
        ``fun`` the function output at the global solution,
        ``xl`` an ordered list of local minima solutions,
        ``funl`` the function output at the corresponding local solutions,
        ``success`` a Boolean flag indicating if the optimizer exited
        successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes. If `polish`
        was employed, then OptimizeResult also contains the `jac` attribute.

    Notes
    -----
    Global optimization using the Topographical Global Optimization (TGO)
    method first proposed by Törn (1990) [1] with the the semi-empirical
    correlation by Hendorson et. al. (2015) [2] for k integer defining the
    k-t matrix.

    The TGO is a clustering method that uses graph theory to generate good
    starting points for local search methods from points distributed uniformly
    in the interior of the feasible set. These points are generated using the
    Sobol (1967) [3] sequence.

    The local search method may be specified using the ``minimizer_kwargs``
    parameter which is inputted to ``scipy.optimize.minimize``. By default
    the ``SLSQP`` method is used. In general it is recommended to use the
    ``SLSQP`` or ``COBYLA`` local minimization if inequality constraints
    are defined for the problem since the other methods do not use constraints.

    Performance can sometimes be improved by increasing ``n``

    Forcing a low ``k_t`` value will

    The primitive polynomials and various sets of initial direction numbers for
    generating Sobol sequences is provided by [4] by Frances Kuo and
    Stephen Joe. The original program sobol.cc is available and described at
    http://web.maths.unsw.edu.au/~fkuo/sobol/ translated to Python 3 by
    Carl Sandrock 2016-03-31

    Examples
    --------
    Consider first the Rosenbrock function

    >>> from scipy.optimize import rosen, tgo
    >>> bounds = [(0,2), (0, 2), (0, 2), (0, 2), (0, 2)]
    >>> result = tgo(rosen, bounds)
    >>> result.x, result.fun
    (array([ 1.,  1.,  1.,  1.,  1.]), 2.9203923741900809e-18)

     Note that bounds sequence determine the dimensionality of the objective
     function and is therefore not optional, however it's possible to specify
     empty bounds example (None, None) which will be converted to large
     float numbers.

    Now consider the Eggholder function
    (https://en.wikipedia.org/wiki/Test_functions_for_optimization)
    >>> from scipy.optimize import tgo
    >>> import numpy as np
    >>> def eggholder(x):
    ...     return (-(x[1] + 47.0)
    ...             * np.sin(np.sqrt(abs(x[0]/2.0 + (x[1] + 47.0))))
    ...             - x[0] * np.sin(np.sqrt(abs(x[0] - (x[1] + 47.0))))
    ...             )
    ...
    >>> bounds = [(-512, 512), (-512, 512)]
    >>> result = tgo(eggholder, bounds)
    >>> result.x, result.fun
    (array([ 512.        ,  404.23180542]), -959.64066272085051)

    ``tgo`` also has a return for any other local minima that was found, these
    can be called using example:
    >>> result.xl, result.funl
    (array([[ 512.        ,  404.23180542],
           [-456.88574619, -382.6233161 ],
           [ 283.07593402, -487.12566542],
           [ 324.99187533,  216.0475439 ],
           [-105.87688985,  423.15324143],
           [-242.97923629,  274.38032063],
           [-414.8157022 ,   98.73012628],
           [ 150.2320956 ,  301.31377513],
           [  91.00922754, -391.28375925],
           [ 361.66626134, -106.96489228]]),
           array([-959.64066272, -786.52599408, -718.16745962, -582.30628005,
           -565.99778097, -559.78685655, -557.85777903, -493.9605115 ,
           -426.48799655, -419.31194957]))

    Now suppose we want to find a larger amount of local minima, this can be
    accomplished by either increasing the amount of sampling points ``n`` or
    defining a lower value for ``k_t`` than the empirical correlation, ex.:

    >>> result2 = tgo(eggholder, bounds, n=1000)
    >>> len(result.xl), len(result2.xl)
    (10, 60)

    >>> from scipy.optimize import tgo
    >>> import numpy as np
    >>> def ackley(x):
    ...     arg1 = -0.2 * np.sqrt(0.5 * (x[0] ** 2 + x[1] ** 2))
    ...     arg2 = 0.5 * (np.cos(2. * np.pi * x[0]) + np.cos(2. * np.pi * x[1]))
    ...     return -20. * np.exp(arg1) - np.exp(arg2) + 20. + np.e
    >>> bounds = [(-5, 5), (-5, 5)]
    >>> result = tgo(ackley, bounds)
    >>> result.x, result.fun
    (array([ -9.02894984e-11,  -9.02900391e-11]), 3.6116221124871117e-10)


    To demonstrate solving problems with non-linear constraints consider the
    following example from [5] (Hock and Schittkowski problem 18):

    Minimize: f = 0.01 * (x_1)**2 + (x_2)**2

    Subject to: x_1 * x_2 - 25.0 >= 0,
                (x_1)**2 + (x_2)**2 - 25.0 >= 0,
                2 <= x_1 <= 50,
                0 <= x_2 <= 50.

    Approx. Answer:
        f([(250)**0.5 , (2.5)**0.5]) = 5.0

    >>> from scipy.optimize import tgo
    >>> def f(x):
    ...     return 0.01 * (x[0])**2 + (x[1])**2
    ...
    >>> def g1(x):
    ...     return x[0] * x[1] - 25.0
    ...
    >>> def g2(x):
    ...     return x[0]**2 + x[1]**2 - 25.0
    ...
    >>> g = (g1, g2)
    >>> bounds = [(2, 50), (0, 50)]
    >>> result = tgo(f, bounds, g_cons=g)
    >>> result.x, result.fun
    (array([ 15.81138847,   1.58113881]), 4.9999999999996252)


    References
    ----------
    .. [1] Törn, A (1990) "Topographical global optimization", Reports on
           Computer Science and Mathematics Ser. A, No 199, 8p. Abo Akademi
           University, Sweden
    .. [2] Henderson, N, de Sá Rêgo, M, Sacco, WF, Rodrigues, RA Jr. (2015) "A
           new look at the topographical global optimization method and its
           application to the phase stability analysis of mixtures",
           Chemical Engineering Science, 127, 151-174
    .. [3] Sobol, IM (1967) "The distribution of points in a cube and the
           approximate evaluation of integrals. USSR Comput. Math. Math. Phys.
           7, 86-112.
    .. [4] S. Joe and F. Y. Kuo (2008) "Constructing Sobol sequences with
           better  two-dimensional projections", SIAM J. Sci. Comput. 30,
           2635-2654
    .. [5] Hoch, W and Schittkowski, K (1981) "Test examples for nonlinear
           programming codes." Lecture Notes in Economics and mathematical
           Systems, 187. Springer-Verlag, New York.
           http://www.ai7.uni-bayreuth.de/test_problem_coll.pdf

    """
    # Initiate TGO class
    TGOc = TGO(func, bounds, args=args, g_funcs=g_cons, g_args=g_args, n=n,
               k_t=k_t, callback=callback, minimizer_kwargs=minimizer_kwargs,
               disp=disp)

    # Generate sampling points
    TGOc.sampling()

    # Find subspace of feasible points
    if g_cons is not None:
        TGOc.subspace()

    # Find topograph
    TGOc.topograph()

    ## Find the optimal k+ topograph
    # Find epsilon_i parameter for current system
    if k_t is None:
        TGOc.K_opt = TGOc.K_optimal()

    # %% Local Search: Find the minimiser float values and func vals.
    TGOc.l_minima()

    # Confirm the routine ran succesfully
    TGOc.res.message = 'Optimization terminated successfully.'
    TGOc.res.succes = True

    # Add local func evals to sampling func evals
    TGOc.res.nfev += TGOc.res.nlfev

    return TGOc.res

# %% Define tgo class
class TGO(object):
    """
    This class implements the tgo routine
    """

    def __init__(self, func, bounds, args=(), g_funcs=None, g_args=(), n=100,
                 k_t=None, callback=None, minimizer_kwargs=None,
                 disp=False):

        self.func = func
        self.bounds = bounds
        self.args = args
        if type(g_funcs) is not tuple and type(g_funcs) is not list:
            self.g_func = (g_funcs,)
        else:
            self.g_func = g_funcs

        self.g_args = g_args
        self.n = n
        self.k_t = k_t
        if k_t is not None:
            self.K_opt = k_t

        self.callback = callback
        self.disp = disp

        # set bounds
        if len(bounds) == 0:
            pass
            #raise ValueError('Error: bounds are required to de')

        abound = numpy.array(bounds, float)
        # Check if obunds are correctly specified
        if abound.ndim > 1:
            bnderr = numpy.where(abound[:, 0] > abound[:, 1])[0]
            # Set none finite values to large floats
            infind = ~numpy.isfinite(abound)
            abound[infind[:, 0], 0] = -1e50#e308
            abound[infind[:, 1], 1] = 1e50#e308
        else:
            bnderr = numpy.where(abound[0] > abound[1])[0]
        if bnderr.any():
            raise ValueError('Error: lb > ub in bounds %s.' %
                             ', '.join(str(b) for b in bnderr))

            # Set none finite values to large floats
            infind = ~numpy.isfinite(abound)
            abound[infind[0], 0] = -1e50#e308
            abound[infind[1], 1] = 1e50#e308

        self.bounds = abound

        # Define constraint function used in local minimisation
        if g_funcs is not None:
             self.min_cons = []
             for g in self.g_func:
                 self.min_cons.append({'type': 'ineq',
                                       'fun' : g})

        # Define local minimization keyword arguments
        if minimizer_kwargs is not None:
            self.minimizer_kwargs = minimizer_kwargs
            if 'args' not in minimizer_kwargs:
                self.minimizer_kwargs['args'] = self.args

            if 'method' not in minimizer_kwargs:
                self.minimizer_kwargs['method'] = 'SLSQP'

            if 'bounds' not in minimizer_kwargs:
                self.minimizer_kwargs['bounds'] = self.bounds

            if 'options' not in minimizer_kwargs:
                minimizer_kwargs['options'] = {'ftol': 1e-12}

            if self.minimizer_kwargs['method'] == 'SLSQP' or \
               self.minimizer_kwargs['method'] == 'COBYLA':
                if 'constraints' not in minimizer_kwargs:
                    minimizer_kwargs['constraints'] = self.min_cons
        else:
            self.minimizer_kwargs = {'args': self.args,
                                     'method': 'SLSQP',
                                     'bounds': self.bounds,
                                     'options': {'ftol': 1e-12}
                                     }
            if g_funcs is not None:
                self.minimizer_kwargs['constraints'] = self.min_cons


            # self.minimizer_kwargs['constraints'] = self.cons

        # Initialize return object
        self.res = scipy.optimize.OptimizeResult()
        self.res.nfev = n  # Include each sampling point as func evaluation
        self.res.nlfev = 0  # Local function evals for all minimisers
        self.res.nljev = 0  # Local jacobian evals for all minimisers



    def sobol_points(self, N, D):
        """
        sobol.cc by Frances Kuo and Stephen Joe translated to Python 3 by
        Carl Sandrock 2016-03-31

        The original program is available and described at
        http://web.maths.unsw.edu.au/~fkuo/sobol/
        """
        with open('new-joe-kuo-6.21201') as f:
            unsigned = "uint64"
            # swallow header
            buffer = next(f)

            L = int(numpy.log(N) // numpy.log(2.0)) + 1

            C = numpy.ones(N, dtype=unsigned)
            for i in range(1, N):
                value = i
                while value & 1:
                    value >>= 1
                    C[i] += 1

            points = numpy.zeros((N, D), dtype='double')

            # XXX: This appears not to set the first element of V
            V = numpy.empty(L + 1, dtype=unsigned)
            for i in range(1, L + 1):
                V[i] = 1 << (32 - i)

            X = numpy.empty(N, dtype=unsigned)
            X[0] = 0
            for i in range(1, N):
                X[i] = X[i - 1] ^ V[C[i - 1]]
                points[i, 0] = X[i] / 2 ** 32

            for j in range(1, D):
                F_int = [int(item) for item in next(f).strip().split()]
                (d, s, a), m = F_int[:3], [0] + F_int[3:]

                if L <= s:
                    for i in range(1, L + 1): V[i] = m[i] << (32 - i)
                else:
                    for i in range(1, s + 1): V[i] = m[i] << (32 - i)
                    for i in range(s + 1, L + 1):
                        V[i] = V[i - s] ^ (
                        V[i - s] >> numpy.array(s, dtype=unsigned))
                        for k in range(1, s):
                            V[i] ^= numpy.array(
                                (((a >> (s - 1 - k)) & 1) * V[i - k]),
                                dtype=unsigned)

                X[0] = 0
                for i in range(1, N):
                    X[i] = X[i - 1] ^ V[C[i - 1]]
                    points[i, j] = X[i] / 2 ** 32  # *** the actual points

            return points

    def sampling(self):
        """
        Generates uniform sampling points in a hypercube and scales the points
        to the bound limits.
        """
        # Generate sampling points.
        #  TODO Assert if func output matches dims. found from bounds
        self.m = len(self.bounds)  # Dimensions

        # Generate uniform sample points in R^m
        self.C = self.sobol_points(self.n, self.m)

        # Distribute over bounds
        # TODO: Find a better way to do this
        for i in range(len(self.bounds)):
            self.C[:, i] = (self.C[:, i] *
                            (self.bounds[i][1] - self.bounds[i][0])
                            + self.bounds[i][0])

        return self.C

    def subspace(self):
        """Find subspace of feasible points from g_func definition"""
        # Subspace of feasible points.
        for g in self.g_func:
            self.C = self.C[g(self.C.T, *self.g_args) >= 0.0]

        #TODO: Check if container is empty fail test or increase n

    def topograph(self):
        """
        Returns the topographical matrix with True boolean values indicating
        positive entries and False ref. values indicating negative values.
        """
        self.Y = scipy.spatial.distance.cdist(self.C, self.C, 'euclidean')
        self.Z = numpy.argsort(self.Y, axis=-1)
        self.A = numpy.delete(self.Z, 0, axis=-1)  # Topographical matrix
                                                   #  without signs
        # Obj. function returns to be used as reference table.:
        self.F = numpy.zeros(numpy.shape(self.C)[0])
        for i in range(numpy.shape(self.C)[0]):
            self.F[i] = self.func(self.C[i,:], *self.args)
        # TODO: see scipy.spatial.KDTree for F lookup?

        # %% Create float value and bool topograph:
        self.H = self.F[self.A] # This replaces all index values in A with the
                           # function result

        self.T = (self.H.T > self.F.T).T  # Topograph with Boolean entries
        return self.T, self.H, self.F

    def k_t_matrix(self, T, k):
        """Returns the k-t topograph matrix"""
        # TODO: Replace delete with simpler array access
        return numpy.delete(T, numpy.s_[k:numpy.shape(T)[1]], axis=-1)

    def minimizers(self, K):
        """Returns the minimizer indexes of a k-t matrix"""
        Minimizers = numpy.all(K, axis=-1)
        # Find data point indexes of minimizers:
        return numpy.where(Minimizers)[0]

    def K_optimal(self):
        """
        Returns the optimal k-t topograph with the semi-empirical correlation
        proposed by Henderson et. al. (2015)
        """
        # TODO: Recheck correct implementation, compare with HS19
        K_1 = self.k_t_matrix(self.T, 1)  # 1-t topograph
        k_1 = len(self.minimizers(K_1))
        k_i = k_1
        i = 2
        while k_1 == k_i:
            K_i = self.k_t_matrix(self.T, i)
            k_i = len(self.minimizers(K_i))
            i += 1

        ep = i * k_i / (k_1 - k_i)
        k_c = numpy.floor((-(ep - 1) + numpy.sqrt((ep - 1.0)**2 + 80.0 * ep))
                          / 2.0)

        k_opt = int(k_c + 1)
        if k_opt > numpy.shape(self.T)[1]:
            # If size of k_opt exceeds t-graph size.
            k_opt = int(numpy.shape(self.T)[1])

        self.K_opt = self.k_t_matrix(self.T, k_opt)
        return self.K_opt

    def l_minima(self):
        """
        Find the local minima using the chosen local minimisation method with
        the minimisers as starting points.
        """
        Min_ind = self.minimizers(self.K_opt)
        self.x_vals = []
        self.Func_min = numpy.zeros_like(Min_ind, dtype=float)

        for i, ind in zip(range(len(Min_ind)), Min_ind):
            # Find minimum x vals
            lres = scipy.optimize.minimize(self.func, self.C[ind, :],
                                            **self.minimizer_kwargs)
            self.x_vals.append(lres.x)
            self.Func_min[i] = lres.fun

            # Local function evals for all minimisers
            self.res.nlfev += lres.nfev

        self.x_vals = numpy.array(self.x_vals)
        # Sort and save
        ind_sorted = numpy.argsort(self.Func_min)  # Sorted indexes in Func_min

        # Save ordered list of minima
        self.res.xl = self.x_vals[ind_sorted]  # Ordered x vals
        self.res.funl = self.Func_min[ind_sorted]  # Ordered fun values

        # Find global of all minimisers
        self.res.x = self.x_vals[ind_sorted[0]]  # Save global minima
        x_global_min = self.x_vals[ind_sorted[0]][0]
        self.res.fun = self.Func_min[ind_sorted[0]]  # Save global fun value
        return x_global_min


    
if __name__ == '__main__':
    pass
    
    
    
 