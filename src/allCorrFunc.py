#!bin/bash/
#@Author: Sean McLaughlin
#This file takes cached halo catalogs, populates them according to some HOD, and calculates the gg cross correlations

import numpy as np
import argparse
from multiprocessing import cpu_count
from os.path import isdir
from halotools.empirical_models import HodModelFactory, TrivialPhaseSpace, NFWPhaseSpace
from halotools.empirical_models import Zheng07Cens, Zheng07Sats
from halotools.sim_manager import CachedHaloCatalog
from halotools.mock_observables import return_xyz_formatted_array, tpcf, tpcf_one_two_halo_decomp, wp
from .redMagicHOD import RedMagicCens, RedMagicSats, StepFuncCens, StepFuncSats
from .myCats import *

N_PTCL = 0 
PI_MAX = 40

RBINS = np.logspace(-1, 1.25, 15)
#RBIN_CENTERS = (RBINS[1:]+RBINS[:-1])/2 #just for plotting

def corrFunc(simname, scale_factor, outputdir,params = {}, n_ptcl = N_PTCL, rbins = RBINS, **kwargs):
    'Calculate the cross correlation for a single catalog at a single scale factor'

    if not isdir(outputdir):
        raise IOError("%s is not a directory."%outputdir)

    cat = cat_dict[simname](**kwargs)
    print str(cat)
    halocat, model = loadHaloAndModel(cat, scale_factor)
    data = popAndCorr(halocat,model, cat, params, n_ptcl, rbins)

    for d, name in zip(data, ['xi_all', 'xi_1h', 'xi_2h', 'wp_all']):
        np.savetxt(outputdir + '%s_%.3f.npy' %(name, scale_factor), d)

def allCorrFunc(simname, outputdir,params = {}, n_ptcl = N_PTCL, rbins = RBINS, **kwargs):
    'Calculates cross correlations for all scale factors cached for one halocatalog'
    if not isdir(outputdir):
        raise IOError("%s is not a directory."%outputdir)

    cat = cat_dict[simname](**kwargs)
    print str(cat)
    for a in cat.scale_factors:
        halocat, model = loadHaloAndModel(cat, a)
        data = popAndCorr(halocat, model, cat, params, n_ptcl, rbins)

        for d, name in zip(data, ['xi_all', 'xi_1h', 'xi_2h', 'wp_all']):
            np.savetxt(outputdir + '%s_%.3f.npy' % (name, a), d)

def loadHaloAndModel(cat, scale_factor):
    '''Return the cached halo catalog and the appropriate HOD model'''
    try:
        idx = cat.scale_factors.index(scale_factor)
    except:
        print 'Provided scale_factor %.3f not cached for %s.'%(scale_factor, cat.simname)
        raise

    halocat = CachedHaloCatalog(simname=cat.simname, halo_finder=cat.halo_finder, version_name=cat.version_name,
                                redshift=cat.redshifts[idx])

    model = HodModelFactory(
        centrals_occupation=RedMagicCens(redshift=cat.redshifts[idx]),
        centrals_profile=TrivialPhaseSpace(redshift=cat.redshifts[idx]),
        satellites_occupation=RedMagicSats(redshift=cat.redshifts[idx]),
        satellites_profile=NFWPhaseSpace(redshift=cat.redshifts[idx]))

    return halocat, model

def popAndCorr(halocat,model, cat, params = {}, n_ptcl = N_PTCL, rbins = RBINS):
    '''Populate a halocat with a model and calculate the tpcf, tpcf_1h, tpcf_2h, and projected corr fun'''
    print 'Min Num Particles: %d\t%d bins'%(n_ptcl, len(rbins))
    model.update(params)#insert new params into model
    # Note: slow
    model.populate_mock(halocat,Num_ptcl_requirement=n_ptcl)

    # Now, calculate with Halotools builtin
    # TODO include the fast version
    x, y, z = [model.mock.galaxy_table[c] for c in ['x', 'y', 'z']]
    pos = return_xyz_formatted_array(x, y, z)*cat.h
    Lbox = model.mock.Lbox*cat.h #remember damn little h's

    xi_all = tpcf(pos, rbins, period=Lbox, num_threads=cpu_count())#TODO change numthreads

    halo_hostid = model.mock.galaxy_table['halo_id']

    xi_1h, xi_2h = tpcf_one_two_halo_decomp(pos,
                                            halo_hostid, rbins,
                                            period=Lbox, num_threads=cpu_count(),
                                            max_sample_size=1e7)

    wp_all = wp(pos, rbins, PI_MAX, period=Lbox, num_threads = cpu_count())

    return xi_all, xi_1h, xi_2h, wp_all #TODO I don't need all these, use kwargs to determine which?

if __name__ == '__main__':
    desc = 'Populate a particular halo catalog and calculate cross correlations. '
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('simname', type=str,
                        help='The name of the simulation to populate. Defaults are stored in the myCats module.')
    parser.add_argument('outputdir', type = str,
                        help='The directory to store the outputs of the calculations. ')
    parser.add_argument('--plot', action = 'store_true',
                        help = 'Determine if plots should be saved along with the calculations.')
    # TODO do I want to have an advanced CLI? Connect to kwargs at all?
    args = parser.parse_args()

    allCorrFunc(args.simname, args.outputdir, args.plot)
