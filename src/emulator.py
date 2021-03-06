# @Author Sean McLaughlin
# Several functions that can train an emulator and make predictions from a trained emulator

from os import path
from itertools import izip
import warnings
from glob import glob
import numpy as np
import scipy.optimize as op
from scipy.linalg import block_diag, det, eigh
import george
from george.kernels import *

from allCorrFunc import RBINS
from doBatchCalls import BOUNDS, N_PER_DIM

#DIRECTORY = '/u/ki/swmclau2/des/EmulatorData/'
DIRECTORY = '/u/ki/swmclau2/des/EmulatorLHC/'
NBINS = len(RBINS) - 1
# In general, the params don't need to be ordered
# However, at this final step consistancy is necessary.
# This list defines that order.
PARAMS = ['logMmin', 'sigma_logM', 'logM0', 'logM1', 'alpha', 'f_c']
#The optimization struggles in higher dimensions
#These are values gleaned from lower dims, good guesses.
INITIAL_GUESSES = {'amp': 0.481, 'logMmin':0.1349,'sigma_logM':0.089,
                   'logM0': 2.0, 'logM1':0.204, 'alpha':0.039,
                   'f_c':0.041, 'r':0.040}

INITIAL_GUESSES_BIAS = {'amp':6, 'r':0.03}
for p in PARAMS:
    INITIAL_GUESSES_BIAS[p] = 0.03

#TODO Harcoding in bias for now, which is very lazy
#I should load it in somewhere, presumably by cosmology in the future.
xi_mm = np.array([1.447577e02,8.069545e01,4.991154e01,2.921271e01,
         1.618499e01,8.596292e00,4.510719e00,2.434982e00,
         1.389996e00,8.395150e-1,5.251286e-1,3.297156e-1,
         2.020197e-1,1.181109e-1,6.463543e-2,3.244891e-2,
         1.454078e-2,5.575480e-3])

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

#Fixing params is not possible with the LHC
def get_training_data(directory=DIRECTORY, independent_variable='xi',nbins=NBINS):
    '''load the GP's x,y, and yerr from the output of paramCube'''
    assert independent_variable in set(['xi', 'bias', 'r2xi'])

    corr_files = sorted(glob(path.join(directory, '*corr*.npy')))
    cov_files = sorted(glob(path.join(directory, '*cov*.npy')))  # since they're sorted, they'll be paired up by params.
    npoints = len(corr_files) * nbins# each file contains NBINS points in r, and each file is a 6-d point

    varied_params = set(PARAMS)# - set(fixed_params.keys())

    ndim = len(varied_params) + 1  # lest we forget r

    x = np.zeros((npoints, ndim))
    y = np.zeros((npoints,))
    # yerr = np.zeros((npoints,))
    ycovs = []

    warned = False
    num_skipped = 0
    num_used = 0
    for idx, (corr_file, cov_file) in enumerate(izip(corr_files, cov_files)):
        params, r, xi, cov = file_reader(corr_file, cov_file)

        if np.any(np.isnan(cov)) or np.any(np.isnan(xi)):
            #TODO I think some of the small bins might all be NaNs, but i'm throwing them out anyway.
            #May have to change the procedure here.
            if not warned:
                warnings.warn('WARNING: NaN detected. Skipping point in %s' % cov_file)
                warned = True
            num_skipped += 1
            continue

        assert nbins == len(r)  # at least it'll throw an error if there's an issue.

        num_used+=1

        # doing some shuffling and stacking
        file_params = []
        #NOTE could do a param ordering here
        for p in PARAMS:
            file_params.append(np.ones((nbins,)) * params[p])

        file_params.append(np.log10(r))

        x[idx * nbins:(idx + 1) *nbins, :] = np.stack(file_params).T
        if independent_variable == 'bias':
            y[idx * nbins:(idx + 1) * nbins] = xi/xi_mm
            ycovs.append(cov/np.outer(xi_mm, xi_mm))
        elif independent_variable == 'xi':
            y[idx * nbins:(idx + 1) * nbins] = np.log10(xi)
        # Approximately true, may need to revisit
        # yerr[idx * NBINS:(idx + 1) * NBINS] = np.sqrt(np.diag(cov)) / (xi * np.log(10))
            ycovs.append(
                cov / (np.outer(xi, xi)*np.log(10)**2))  # I think this is right, extrapolating from the above.
            #doing some voodoo on the cov mat.
            '''
            norm = ycovs[-1][0,0]
            mat = ycovs[-1]/norm
            w,v = eigh(mat)
            w[w**2< np.sqrt(2.0/125)] = 0.5
            new_mat = np.dot(v, np.dot(np.diag(w), v.T))
            ycovs[-1] = new_mat*norm
            '''
        else: #r2xi
            y[idx * nbins:(idx + 1) * nbins] = xi*r*r
            ycovs.append(cov*np.outer(r,r))#

    ycov = block_diag(*ycovs)
    #ycov = np.sqrt(np.diag(ycov))
    if num_used==0:
        raise RuntimeError('One of your parameters wasn\'t exact!')
    else:
        print '%d Points used for training.'%num_used

    #print '\nSkipped %.2f %% of points due to NaNs.'%(100.0*(num_skipped)/(num_used+num_skipped))

    # remove rows that were skipped due to the fixed thing
    # NOTE: HACK
    # a reshape may be faster.
    zeros_slice = np.all(x != 0.0, axis=1)

    return x[zeros_slice], y[zeros_slice], ycov

#Not sure this will work at all in an LHC scheme.
def get_plot_data(em_params,fixed_params, directory=DIRECTORY, independent_variable = 'xi',nbins=NBINS):
    '''Load truths from the output of paramCube for plotting alongside the GP emulations.'''
    assert independent_variable in set(['xi', 'bias', 'r2xi'])

    assert len(em_params)+len(fixed_params) +1 == len(PARAMS)

    corr_files = sorted(glob(path.join(directory, '*corr*.npy')))
    cov_files = sorted(glob(path.join(directory, '*cov*.npy')))  # since they're sorted, they'll be paired up by params.
    npoints = len(corr_files)# each file contains NBINS points in r, and each file is a 6-d point

    log_r = np.zeros((npoints,nbins ))
    y = np.zeros((npoints,nbins))
    y_err = np.zeros((npoints,nbins))

    warned = False
    num_skipped = 0
    num_used = 0
    for idx, (corr_file, cov_file) in enumerate(izip(corr_files, cov_files)):
        params, r, xi, cov = file_reader(corr_file, cov_file)

        # skip values that aren't where we've fixed them to be.
        # It'd be nice to do this before the file I/O. Not possible without putting all info in the filename.
        # or, a more nuanced file structure
        #TODO check if a fixed_param is not one of the options
        if any(params[key] != val for key, val in fixed_params.iteritems()):
            continue

        if any(params[key] != val for key, val in em_params.iteritems()):
            continue

        if np.any(np.isnan(cov)) or np.any(np.isnan(xi)):
            if not warned:
                warnings.warn('WARNING: NaN detected. Skipping point in %s' % cov_file)
                warned = True
            num_skipped += 1
            continue

        assert nbins == len(r)  # at least it'll throw an error if there's an issue.

        num_used+=1

        log_r[idx] = np.log10(r) 
        if independent_variable == 'bias':
            y[idx] = xi/xi_mm
            y_err[idx] = np.sqrt(np.diag(cov))/(xi_mm)  # I think this is right, extrapolating from the above.
        elif independent_variable == 'xi':
            y[idx] = np.log10(xi)
            # Approximately true, may need to revisit
            # yerr[idx * NBINS:(idx + 1) * NBINS] = np.sqrt(np.diag(cov)) / (xi * np.log(10))
            y_err[idx] = np.sqrt(np.diag(cov))/(xi*np.log(10))  # I think this is right, extrapolating from the above.
        else: #r2xi
            y[idx] = xi*r*r
            y_err[idx] = np.sqrt(np.diag(cov))*r  # I think this is right, extrapolating from the above.


    # remove rows that were skipped due to the fixed thing
    # NOTE: HACK
    # a reshape may be faster.
    zeros_slice = np.all(y != 0.0, axis=1)

    return log_r[zeros_slice], y[zeros_slice], y_err[zeros_slice] 

#TODO figure out what to do about the nbins thing; i'm not happy with it.
def build_emulator(directory=DIRECTORY,independent_variable = 'xi',nbins=NBINS):
    '''Actually build the emulator. '''
    from time import time
    ndim = len(PARAMS) + 1  # include r
    t0 = time() 
    x, y, y_cov = get_training_data(directory, independent_variable,nbins)
    print 'Get Data: %.2f'%(time()-t0)

    #TODO have to make initial guesses dict now
    ig = INITIAL_GUESSES_BIAS if independent_variable=='bias' else INITIAL_GUESSES
    
    metric = []
    for p in PARAMS:
        metric.append(ig[p])

    metric.append(ig['r'])

    a = ig['amp'] 
    kernel = a * ExpSquaredKernel(metric, ndim=ndim)
    gp = george.GP(kernel)
    #gp = george.GP(kernel, solver=george.HODLRSolver, nleaf=x.shape[0]+1,tol=1e-18)

    t0 = time()
    gp.compute(x, np.sqrt(np.diag(y_cov)))  # NOTE I'm using a modified version of george!
    print 'Compute Time: %.2f'%(time()-t0)

    print 'Initial Params: '
    for p in gp.kernel:
        print '%.6f'%np.exp(p)
    print 

    return gp, y, y_cov

def train_emulator(gp, y):
    '''Attempt to optimize the emulator!'''
    from time import time
    y_hat = y.mean() # NOTE consider doing a copy here rather than modifying y directly
    y -= y_hat

    ndim = gp._x.shape[1]

    def nll(p):
        t0 = time()
        # Update the kernel parameters and compute the likelihood.
        # params are log(a) and log(m)
        gp.kernel[:] = p
        ll = gp.lnlikelihood(y, quiet=True)

        # The scipy optimizer doesn't play well with infinities.
        print 'Nll\t%f\t%.2f s'%(ll, time()-t0)
        return -ll if np.isfinite(ll) else 1e25

    # And the gradient of the objective function.
    def grad_nll(p):
        # Update the kernel parameters and compute the likelihood.
        from time import time
        t0 = time()
        gp.kernel[:] = p
        output = -gp.grad_lnlikelihood(y, quiet=True)
        print 'Grad Nll\t%f\t%.2f s'%(output.mean(), time()-t0)
        #return -gp.grad_lnlikelihood(y, quiet=True)
        return output

    p0 = gp.kernel.vector
    t0 = time()
    results = op.minimize(nll, p0, jac=grad_nll, method='TNC', bounds = [(np.log(0.01), np.log(10)) for i in xrange(ndim+1)],options={'maxiter':50})
    print 'Training time: %.2f s'%(time()-t0)

    if not results.success:
        warnings.warn('WARNING: GP Optimization failed!')
        print 'GP Optimization Failed.'#The warning doesn't warn more than once. 

    gp.kernel[:] = results.x
    gp.recompute()
    print 'GP Params: '
    for p in results.x:
        print '%.6f'%np.exp(p)
    print
    print 'Computed: %s'%gp.computed
    y+=y_hat
    return #don't need to return anything!

# unsure on the design here. I'd like to combin the x,y* into one thing each, but idk how that's easy
# just having them be len(x)==1 dicts seems silly.
def emulate(gp, em_y, em_params, x_param, x_points, y_param=None, y_points=None):
    # check that all params have been accounted for!
    # could wrap this in a try block to make it more informative
    input_params = set(em_params) | set([x_param])
    if y_param is not None:
        input_params.add(y_param)
    assert len(input_params) == gp._x.shape[1]  # check dimenstionality
    assert all([i in PARAMS or i == 'r' for i in input_params])
    # used to have a param check but forgot about the whole dropping params thing.

    assert all(y is None for y in [y_param, y_points]) or all(y is not None for y in [y_param, y_points])

    em_y_hat = em_y.mean()
    em_y-=em_y_hat

    if y_param is None:
        # We only have to vary one parameter, as given by x_points
        t_list = []
        for p in PARAMS:
            if p in em_params:
                t_list.append(np.ones_like(x_points) * em_params[p])
            elif p == x_param:
                t_list.append(x_points)
            else:
                continue
        # adding 'r' in as a special case
        if 'r' in em_params:
            t_list.append(np.ones_like(x_points) * em_params['r'])
        elif 'r' == x_param:
            t_list.append(x_points)

        t = np.stack(t_list).T

        # TODO mean subtraction?
        mu, cov = gp.predict(em_y, t)

        # TODO return std or cov?
        # TODO return r's too? Just have those be passed in?
        em_y+=em_y_hat
        return mu+em_y_hat, np.diag(cov)
    else:
        output = []
        assert len(y_points) <= 20  # y_points has a limit, otherwise this'd be crazy


        for y in y_points:  # I thought this had a little too mcuh copying, but
            # this is the best wayt ensure the ordering is consistent.
            t_list = []
            for p in PARAMS:
                if p in em_params:
                    t_list.append(np.ones_like(x_points) * em_params[p])
                elif p == x_param:
                    t_list.append(x_points)
                elif p == y_param:
                    t_list.append(np.ones_like(x_points) * y)
                else:
                    continue

            if 'r' in em_params:
                t_list.append(np.ones_like(x_points) * em_params['r'])
            elif 'r' == x_param:
                t_list.append(x_points)
            elif 'r' == y_param:
                t_list.append(np.ones_like(x_points) * y)

            t = np.stack(t_list).T

            mu, cov = gp.predict(em_y, t)
            output.append((mu+em_y_hat, np.sqrt(np.diag(cov))))
        em_y+=em_y_hat
        return output


def emulate_wrt_r(gp, em_y, em_params, rpoints, y_param=None, y_points=None):
    '''simplified version of the above that implements the most common case, 1-D emulation in r'''
    assert 'r' not in em_params
    return emulate(gp,em_y,em_params, x_param='r', x_points=rpoints,
                    y_param=y_param,y_points=y_points)

if __name__ == '__main__':
    from time import time

    y_param = 'logMmin'
    ep = ['sigma_logM', 'logM0', 'logM1', 'alpha']#, 'f_c'] 

    emulation_point = [('f_c', 0.233),
                       ('logM0', 12.0), ('sigma_logM', 0.533), ('alpha', 1.083),
                       ('logM1', 13.5), ('logMmin', 12.233)]

    fiducial_point = {'logM0': 12.20, 'logM1': 13.7, 'alpha': 1.02,
                      'logMmin': 12.1, 'f_c': 0.19, 'sigma_logM': 0.46}

    yp = np.linspace(BOUNDS[y_param][0],BOUNDS[y_param][1], num=N_PER_DIM)

    rpoints = np.linspace(-1, 1.5, 200)

    RBIN_CENTERS = ((RBINS[:-1] + RBINS[1:]) / 2)
    flip_rpoints = np.log10(RBIN_CENTERS)
    flip_yp = np.linspace(BOUNDS[y_param][0], BOUNDS[y_param][1], num=200)
    '''
    for i in xrange(len(ep)):
        fixed_params = {key: val for key, val in emulation_point}

        em_params = {}
        _ep = ep[:]
        skipped = _ep.pop(i)
        print i, skipped

        for param in _ep:
            em_params[param] = fixed_params[param]
            del fixed_params[param]

        del fixed_params[y_param]

        t0 = time()
        gp,xi,xi_cov = build_emulator(fixed_params)
        print 'Build time: %.2f seconds'%(time()-t0)

        t1 = time()
        train_emulator(gp, xi)
        print 'Train time: %.2f seconds' % (time() - t1)

        outputs = emulate_wrt_r(gp,xi,em_params,rpoints,y_param=y_param,y_points=yp)
        outputs = np.stack(outputs)
        print 'Total time: %.2f seconds'%(time()-t0)

        flip_outputs = emulate(gp, xi, em_params, x_param=y_param, x_points=flip_yp, y_param='r', y_points=flip_rpoints)
        flip_outputs = np.stack(flip_outputs)
        
        plot_outputs = get_plot_data(em_params, fixed_params)
        plot_outputs = np.stack(plot_outputs)#.reshape((-1,3))

        output_dir = '/u/ki/swmclau2/des/EmulatorTest'
        np.save(path.join(output_dir, 'output_%d.npy'%i), outputs)
        np.save(path.join(output_dir, 'flip_output_%d.npy'%i), flip_outputs)
        np.save(path.join(output_dir, 'plot_output_%d.npy'%i), plot_outputs)
        print '*-'*30
    '''

    fixed_params = {key: val for key, val in emulation_point}

    em_params = {}

    for param in ep:
        em_params[param] = fixed_params[param]
        del fixed_params[param]

    del fixed_params[y_param]
    independent_variable = 'xi'

    t0 = time()
    gp, xi, xi_cov = build_emulator(independent_variable = independent_variable)
    print 'Build time: %.2f seconds' % (time() - t0)

    t1 = time()
    #train_emulator(gp, xi)
    print 'Train time: %.2f seconds' % (time() - t1)

    outputs = emulate_wrt_r(gp, xi, em_params, rpoints, y_param=y_param, y_points=yp)
    outputs = np.stack(outputs)
    print 'Total time: %.2f seconds' % (time() - t0)

    flip_outputs = emulate(gp, xi, em_params, x_param=y_param, x_points=flip_yp, y_param='r', y_points=flip_rpoints)
    flip_outputs = np.stack(flip_outputs)

    plot_outputs = get_plot_data(em_params, fixed_params,independent_variable = independent_variable)
    plot_outputs = np.stack(plot_outputs)  # .reshape((-1,3))

    output_dir = '/u/ki/swmclau2/des/EmulatorTest'
    if not independent_variable == 'bias':
        np.save(path.join(output_dir, 'output_all.npy'), outputs)
        np.save(path.join(output_dir, 'flip_output_all.npy'), flip_outputs)
        np.save(path.join(output_dir, 'plot_output_all.npy'), plot_outputs)
    else:
        np.save(path.join(output_dir, 'bias_output_all.npy'), outputs)
        np.save(path.join(output_dir, 'bias_flip_output_all.npy'), flip_outputs)
        np.save(path.join(output_dir, 'bias_plot_output_all.npy'), plot_outputs)
