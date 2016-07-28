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


def lnlike(theta, xi, cov, gp, training_xi, fp):
    idx = 0
    # a bit of a hack
    # a way to combine the fixed and varying params in one spot.
    for p in PARAMS:
        if fp[p] is None:
            fp[p] = theta[idx]
        idx += 1

    log_xi_em, err_em = emulate_wrt_r(gp, training_xi, fp, RBINS)
    # TODO figure out how to use err_am
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


def run_mcmc(xi, cov, fixed_params={}):
    nwalkers = 2000
    nsteps = 200
    nburn = 50  # TODO convergence measures
    ndim = len(PARAMS) - len(fixed_params)
    num_threads = cpu_count()

    assert all(x is None for x in (xi, cov)) or all(x is not None for x in (xi, cov))
    assert xi.shape[0] == cov.shape[0] and cov.shape[1] == cov.shape[0]

    gp, training_xi = build_emulator(fixed_params)
    # I'm not happy about how this works, but allows me to carry information into the liklihood
    fp = {}
    for p in PARAMS:
        if p not in fp:
            fp[p] = None

    sampler = mc.EnsembleSampler(nwalkers, ndim, lnprob, threads=num_threads, args=(xi, cov, gp, training_xi, fp))

    pos0 = np.zeros((nwalkers, ndim))
    idx = 0
    for p in PARAMS:
        if p in fixed_params:
            continue
        pos0[:, idx] = np.random.uniform(BOUNDS[p][0], BOUNDS[p][1], size = nwalkers)
        idx+=1

    sampler.run_mcmc(pos0, nsteps)

    #Note, still an issue of param label ordering here.
    chain = sampler.chain[:, nburn:, :].reshape((-1, ndim))

    return chain