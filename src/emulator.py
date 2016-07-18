#@Author Sean McLaughlin
#Several functions that can train an emulator and make predictions from a trained emulator

from os import path
from itertools import izip
from glob import glob
import numpy as np
import scipy.optimize as op
import george
from george.kernels import *

DIRECTORY = '/u/ki/swmclau2/des/EmulatorData/'
NBINS = 19
#In general, the params don't need to be ordered
#However, at this final step consistancy is necessary.
#This list defines that order.
PARAMS = ['logMmin', 'sigma_logM', 'logM0', 'logM1', 'alpha', 'f_c']

def file_reader(corr_file, cov_file):
    '''Data is stored in two formats. Numpy writes r and xi, and the comments hold what parameters they describe.
    I'll use numpy to load the r and xi and read the file manually to get the params.'''

    #TODO check filename exists, etc.
    r, xi = np.loadtxt(corr_file)# not sure if this will work, might nead to transpose
    cov = np.loadtxt(cov_file)
    params = {}
    with open(corr_file) as f:
        for i, line in enumerate(f):
            if line[0] != '#' or i < 2:
                continue #only looking at comments, and first two lines don't have params. Note: Does have cosmo!
            splitLine = line.strip('# \n').split(':')#split into key val pair
            params[splitLine[0]] = float(splitLine[1])

    return params, r, xi, cov

def get_training_data(fixed_params = {}):
    '''load the GP's x,y, and yerr from the output of paramCube'''
    corr_files = sorted(glob(path.join(DIRECTORY, '*corr*.npy')))
    cov_files = sorted(glob(path.join(DIRECTORY, '*cov*.npy')))  # since they're sorted, they'll be paired up by params.
    npoints = len(corr_files) * NBINS  # each file contains NBINS points in r, and each file is a 6-d point

    varied_params = set(PARAMS)-set(fixed_params.keys())

    ndim = len(varied_params)

    x = np.zeros((npoints, ndim))
    y = np.zeros((npoints,))
    yerr = np.zeros((npoints,))

    for idx, (corr_file, cov_file) in enumerate(izip(corr_files, cov_files)):
        params, r, xi, cov = file_reader(corr_file, cov_file)

        # skip values that aren't where we've fixed them to be.
        # It'd be nice to do this before the file I/O. Not possible without putting all info in the filename.
        #or, a more nuanced file structure
        if any(params[key] != val for key, val in fixed_params.iteritems()):
            continue

        #doing some shuffling and stacking
        file_params = [np.ones((NBINS,)) * params[p] for p in varied_params]
        file_params.append(np.log10(r))

        x[idx * NBINS:(idx + 1) * NBINS, :] = np.stack(file_params).T
        y[idx * NBINS:(idx + 1) * NBINS] = np.log10(xi)
        #Approximately true, may need to revisit
        yerr[idx * NBINS:(idx + 1) * NBINS] = np.sqrt(np.diag(cov)) / (xi * np.log(10))

    #remove rows that were skipped due to the fixed thing
    #NOTE: HACK
    #a reshape may be faster.
    zeros_slice = np.all(x != 0.0, axis=1)

    return x[zeros_slice],y[zeros_slice],yerr[zeros_slice]


def build_emulator(fixed_params = {}):
    ndim = len(PARAMS) - len(fixed_params)

    x,y,yerr = get_training_data(fixed_params)

    metric = (1.0 for i in xrange(ndim)) #could make better guesses:
    a = 1e5
    kernel = a * ExpSquaredKernel(metric, ndim=2)
    gp = george.GP(kernel)
    gp.compute(x, yerr)

    def nll(p):
        # Update the kernel parameters and compute the likelihood.
        # params are log(a) and log(m)
        gp.kernel[:] = p
        ll = gp.lnlikelihood(y, quiet=True)

        # The scipy optimizer doesn't play well with infinities.
        return -ll if np.isfinite(ll) else 1e25

    # And the gradient of the objective function.
    def grad_nll(p):
        # Update the kernel parameters and compute the likelihood.
        gp.kernel[:] = p
        return -gp.grad_lnlikelihood(y, quiet=True)

    p0 = gp.kernel.vector
    results = op.minimize(nll, p0, jac=grad_nll)

    if not results.success:
        print 'WARNING: GP Optimization failed!'
        #TODO a real warning

    gp.kernel[:] = results.x

    return gp , y
