# @Author Sean McLaughlin
# Several functions that can train an emulator and make predictions from a trained emulator

from os import path
from itertools import izip
from glob import glob
import numpy as np
import scipy.optimize as op
import george
from george.kernels import *

from allCorrFunc import RBINS

DIRECTORY = '/u/ki/swmclau2/des/EmulatorData/'
NBINS = len(RBINS) -1
# In general, the params don't need to be ordered
# However, at this final step consistancy is necessary.
# This list defines that order.
PARAMS = ['logMmin', 'sigma_logM', 'logM0', 'logM1', 'alpha', 'f_c']


def file_reader(corr_file, cov_file):
    '''Data is stored in two formats. Numpy writes r and xi, and the comments hold what parameters they describe.
    I'll use numpy to load the r and xi and read the file manually to get the params.'''

    assert path.exists(corr_file) and path.exists(cov_file)
    r, xi = np.loadtxt(corr_file)  # not sure if this will work, might nead to transpose
    cov = np.loadtxt(cov_file)
    params = {}
    with open(corr_file) as f:
        for i, line in enumerate(f):
            if line[0] != '#' or i < 2:
                continue  # only looking at comments, and first two lines don't have params. Note: Does have cosmo!
            splitLine = line.strip('# \n').split(':')  # split into key val pair
            params[splitLine[0]] = float(splitLine[1])

    return params, r, xi, cov


def get_training_data(fixed_params={}, directory=DIRECTORY):
    '''load the GP's x,y, and yerr from the output of paramCube'''
    corr_files = sorted(glob(path.join(DIRECTORY, '*corr*.npy')))
    cov_files = sorted(glob(path.join(DIRECTORY, '*cov*.npy')))  # since they're sorted, they'll be paired up by params.
    npoints = len(corr_files) * NBINS  # each file contains NBINS points in r, and each file is a 6-d point

    varied_params = set(PARAMS) - set(fixed_params.keys())

    ndim = len(varied_params) +1 #lest we forget r

    x = np.zeros((npoints, ndim))
    y = np.zeros((npoints,))
    yerr = np.zeros((npoints,))

    for idx, (corr_file, cov_file) in enumerate(izip(corr_files, cov_files)):
        params, r, xi, cov = file_reader(corr_file, cov_file)

        # skip values that aren't where we've fixed them to be.
        # It'd be nice to do this before the file I/O. Not possible without putting all info in the filename.
        # or, a more nuanced file structure
        #TODO check if a fixed_param is not one of the options
        if any(params[key] != val for key, val in fixed_params.iteritems()):
            continue

        assert NBINS == len(r)#at least it'll throw an error if there's an issue.

        # doing some shuffling and stacking
        file_params = []
        for p in PARAMS:
            if p in fixed_params:
                continue
            file_params.append(np.ones((NBINS,)) * params[p])

        file_params.append(np.log10(r))

        x[idx * NBINS:(idx + 1) * NBINS, :] = np.stack(file_params).T
        y[idx * NBINS:(idx + 1) * NBINS] = np.log10(xi)
        # Approximately true, may need to revisit
        yerr[idx * NBINS:(idx + 1) * NBINS] = np.sqrt(np.diag(cov)) / (xi * np.log(10))

    # remove rows that were skipped due to the fixed thing
    # NOTE: HACK
    # a reshape may be faster.
    zeros_slice = np.all(x != 0.0, axis=1)

    return x[zeros_slice], y[zeros_slice], yerr[zeros_slice]


def build_emulator(fixed_params={}, directory=DIRECTORY):
    '''Actually build the emulator. '''

    ndim = len(PARAMS) - len(fixed_params) +1 #include r
    x, xi, xi_err = get_training_data(fixed_params,directory)
    print x.shape 

    metric = [1.0 for i in xrange(ndim)]  # could make better guesses:
    a = 1e5
    kernel = a * ExpSquaredKernel(metric, ndim=ndim)
    gp = george.GP(kernel)
    #In the test module some of the errors are NaNs
    #TODO remove this in the main implementation
    xi_err[np.isnan(xi_err)] = 1.0
    gp.compute(x, xi_err)

    def nll(p):
        # Update the kernel parameters and compute the likelihood.
        # params are log(a) and log(m)
        gp.kernel[:] = p
        ll = gp.lnlikelihood(xi, quiet=True)

        # The scipy optimizer doesn't play well with infinities.
        return -ll if np.isfinite(ll) else 1e25

    # And the gradient of the objective function.
    def grad_nll(p):
        # Update the kernel parameters and compute the likelihood.
        gp.kernel[:] = p
        return -gp.grad_lnlikelihood(xi, quiet=True)

    p0 = gp.kernel.vector
    results = op.minimize(nll, p0, jac=grad_nll)

    if not results.success:
        import warnings
        warnings.warn('WARNING: GP Optimization failed!')

    gp.kernel[:] = results.x

    return gp, xi


# unsure on the design here. I'd like to combin the x,y* into one thing each, but idk how that's easy
# just having them be len(x)==1 dicts seems silly.
#TODO clarity of log_xi, xi
#TODO fixed params is different than the above; clarify
def emulate(gp, xi, fixed_params, x_param, x_points, y_param=None, y_points=None):
    # check that all params have been accounted for!
    # could wrap this in a try block to make it more informative
    input_params = set(fixed_params) | set([x_param])
    if y_param is not None:
        input_params.add(y_param)
    assert len(input_params) == gp._x.shape[1]  # check dimenstionality
    assert all([i in PARAMS or i == 'r' for i in input_params])
    # used to have a param check but forgot about the whole dropping params thing.

    assert all(y is None for y in [y_param, y_points]) or all(y is not None for y in [y_param, y_points])

    if y_param is None:
        # We only have to vary one parameter, as given by x_points
        t_list = []
        for p in PARAMS:
            if p in fixed_params:
                t_list.append(np.ones_like(x_points) * fixed_params[p])
            elif p == x_param:
                t_list.append(x_points)
            else:
                continue
        #adding 'r' in as a special case
        if 'r' in fixed_params:
            t_list.append(np.ones_like(x_points)*fixed_params['r'])
        elif 'r'==x_param:
            t_list.append(x_points)

        t = np.stack(t_list).T

        # TODO mean subtraction?
        print xi.shape, t.shape, gp._x.shape
        mu, cov = gp.predict(xi, t)

        # TODO return std or cov?
        return mu, np.diag(cov)
    else:
        output = []
        assert len(y_points) <= 20 #y_points has a limit, otherwise this'd be crazy

        for y in y_points: # I thought this had a little too mcuh copying, but
            #this is the best wayt ensure the ordering is consistent.
            t_list = []
            for p in PARAMS:
                if p in fixed_params:
                    t_list.append(np.ones_like(x_points) * fixed_params[p])
                elif p == x_param:
                    t_list.append(x_points)
                elif p == y_param:
                    t_list.append(np.ones_like(x_points)*y)
                else:
                    continue

            if 'r' in fixed_params:
                t_list.append(np.ones_like(x_points) * fixed_params['r'])
            elif 'r' == x_param:
                t_list.append(x_points)
            elif 'r' == y_param:
                t_list.append(np.ones_like(x_points) * y)

            t = np.stack(t_list)

            mu, cov = gp.predict(xi, t)
            output.append((mu, np.sqrt(np.diag(cov))))
        return output

def emulate_wrt_r(gp,xi,fixed_params, rbins):
    '''simplified version of the above that implements the most common case, 1-D emulation in r'''
    assert 'r' not in fixed_params
    return emulate(gp, xi, fixed_params, x_param='r', x_points=np.log10((rbins[1:]+rbins[:-1])/2))
