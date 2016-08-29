#!/bin/bash
# @Author Sean McLaughlin
# This module will use the emulator to run MCMC to recover HOD parameters.

from multiprocessing import cpu_count
import numpy as np
import emcee as mc
from emulator import build_emulator, emulate_wrt_r, PARAMS
from myCats import cat_dict
from allCorrFunc import loadHaloAndModel, popAndCorr, RBINS, MIN_PTCL
from doBatchCalls import BOUNDS


# TODO clarity of log_xi v. xi
def mock_truth(simname, scale_factor, HOD='redMagic', true_params={}, min_ptcl=MIN_PTCL, rbins=RBINS, **kwargs):
    '''Generates mock "true" data to recover a given set of params.'''
    cat = cat_dict[simname](**kwargs)
    halocat, model = loadHaloAndModel(cat, HOD, scale_factor)
    xi, cov = popAndCorr(halocat, model, cat, true_params, True, min_ptcl, rbins)

    return xi, cov


def lnprior(theta):
    # for now, just a uniform over the params. Something more intelligent later
    for key, val in zip(PARAMS, theta):
        if val < BOUNDS[key][0] or val > BOUNDS[key][1]:
            return -np.inf
    return 0

#TODO change xi to something more generic
def lnlike(theta, xi, rpoints, cov, gp, training_xi):
    em_params = dict(zip(PARAMS, theta))

    log_xi_em, err_em = emulate_wrt_r(gp, training_xi, em_params, rpoints)
    # TODO figure out how to use err_em
    xi_em = 10 ** log_xi_em  # log to normal
    invcov = np.lingalg.inv(cov)
    delta_xi = xi - xi_em
    chi2 = -0.5 * np.sum(np.dot(delta_xi, np.dot(invcov, delta_xi)))
    return chi2


def lnprob(theta, **args):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, **args)


#TODO add bias, etc. options
def run_mcmc(xi,rpoints, cov, **kwargs):
    nwalkers = 2000
    nsteps = 10
    nburn = 0#50  # TODO convergence measures
    ndim = len(PARAMS) +1
    num_threads = cpu_count()

    assert xi.shape[0] == cov.shape[0] and cov.shape[1] == cov.shape[0]
    assert xi.shape[0] == rpoints.shape[0]

    gp, training_xi = build_emulator(**kwargs)
    # I'm not happy about how this works, but allows me to carry information into the liklihood

    sampler = mc.EnsembleSampler(nwalkers, ndim, lnprob,
                                 threads=num_threads, args=(xi,rpoints, cov, gp, training_xi))

    pos0 = np.zeros((nwalkers, ndim))
    for idx, p in enumerate(PARAMS):
        pos0[:, idx] = np.random.uniform(BOUNDS[p][0], BOUNDS[p][1], size = nwalkers)

    sampler.run_mcmc(pos0, nsteps)

    #Note, still an issue of param label ordering here.
    chain = sampler.chain[:, nburn:, :].reshape((-1, ndim))

    return chain

if __name__ == '__main__':
    #TODO has this be argv
    simname = 'chinchilla'
    scale_factor = 1/(1+0.5) #z = 0.5
    fiducial_point = {'logM0': 12.20, 'logM1': 13.7, 'alpha': 1.02,
                      'logMmin': 12.1, 'f_c': 0.19, 'sigma_logM': 0.46}
    xi, cov = mock_truth(simname, scale_factor,true_params=fiducial_point)
    rpoints = (RBINS[:1]+RBINS[:-1])/2
    chain = run_mcmc(xi, rpoints, cov)
    print chain.mean(axis =1)
