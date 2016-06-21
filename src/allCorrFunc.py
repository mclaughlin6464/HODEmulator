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
from halotools.mock_observables import return_xyz_formatted_array,tpcf, tpcf_jackknife, tpcf_one_two_halo_decomp, wp
from redMagicHOD import RedMagicCens, RedMagicSats, StepFuncCens, StepFuncSats
from myCats import *

from memory_profiler import profile

N_PTCL = 0 
PI_MAX = 40

RBINS = np.logspace(-1, 1.7, 20)
RBIN_CENTERS = (RBINS[1:]+RBINS[:-1])/2

#TODO will need ways to pass params into the model when populating. Could just use kwargs, but how to separate cat kwargs?
#See other branch
def corrFunc(simname, scale_factor, outputdir, plot = False,logMmin = 12.1,  **kwargs):
    'Calculate the cross correlation for a single catalog at a single scale factor'
    cat = cat_dict[simname](**kwargs)
    _corrFunc(cat, scale_factor, outputdir, plot, logMmin = logMmin)

def allCorrFunc(simname, outputdir, plot = False, **kwargs):
    'Calculates cross correlations for all scale factors cached for one halocatalog'
    cat = cat_dict[simname](**kwargs)
    for a in cat.scale_factors:
        _corrFunc(cat, a, outputdir, plot)

@profile
def _corrFunc(cat, scale_factor, outputdir, plot = False, logMmin = 12.1):
    'Helper function that uses the built in cat object'

    print str(cat)
    print 'Min Num Particles: %d'%N_PTCL

    if not isdir(outputdir):
        raise IOError("%s is not a directory"%outputdir)

    try:
        idx = cat.scale_factors.index(scale_factor)
    except:
        print 'Provided scale_factor %.3f not cached for %s.'%(scale_factor, cat.simname)
        raise

    #Note: Confusing name between cat and halocat. Consider changing.
    halocat = CachedHaloCatalog(simname = cat.simname, halo_finder = cat.halo_finder,version_name = cat.version_name, redshift = cat.redshifts[idx])

    model = HodModelFactory(
        centrals_occupation=RedMagicCens(redshift=cat.redshifts[idx]),
        #centrals_occupation=StepFuncCens(redshift=cat.redshifts[idx]),
        centrals_profile=TrivialPhaseSpace(redshift=cat.redshifts[idx]),
        satellites_occupation=RedMagicSats(redshift=cat.redshifts[idx],modulate_with_cenocc = True),#TODO consider putting this in redMagicSats intself
        #satellites_occupation=StepFuncSats(redshift=cat.redshifts[idx]),
        satellites_profile=NFWPhaseSpace(redshift=cat.redshifts[idx]))

    model.param_dict['logMmin'] = logMmin# - np.log10(cat.h) #this is gonna have to obviously be changed when merged with the other branch.

    print model.param_dict

    #Note: slow
    model.populate_mock(halocat, Num_ptcl_requirement = N_PTCL) 

    #Now, calculate with Halotools builtin
    #TODO include the fast version
    print 1
    x, y, z = [model.mock.galaxy_table[c] for c in ['x','y','z'] ]
    #mask = model.mock.galaxy_table['halo_mvir'] < 1e15/cat.h
    pos = return_xyz_formatted_array(x,y,z)#, mask = mask)
    print 2
    #TODO N procs
    print cpu_count()
    xi_all = tpcf(pos*cat.h, RBINS, period = model.mock.Lbox*cat.h, num_threads =  cpu_count())

    #TODO ways to decide which of these to call
    #randoms = np.random.random(pos.shape)*model.mock.Lbox*cat.h
    #xi_all, xi_cov = tpcf_jackknife(pos*cat.h,randoms, RBINS, period = model.mock.Lbox*cat.h, num_threads =  cpu_count())
    print 3
    halo_hostid = model.mock.galaxy_table['halo_id']
    print 4

    xi_1h, xi_2h = tpcf_one_two_halo_decomp(pos*cat.h,
                    halo_hostid, RBINS,
                    period = cat.h*model.mock.Lbox, num_threads =  cpu_count(),
                    max_sample_size = 1e7)
    print 5

    #wp_all = wp(pos*cat.h, RBINS, PI_MAX, period=model.mock.Lbox*cat.h, num_threads = cpu_count())

    output = np.stack([RBIN_CENTERS, xi_all, xi_1h, xi_2h])
    print 6

    np.savetxt(outputdir + 'corr_%.3f_default_mm_%.2f.npy' %(scale_factor, logMmin), output)
    #np.savetxt(outputdir + 'xi_cov_%.3f_default_125_2048.npy' %(scale_factor), xi_cov)

    #np.savetxt(outputdir + 'wp_all_%.3f_default.npy' %(scale_factor), wp_all)

if __name__ == '__main__':
    desc = 'Populate a particular halo catalog and calculate cross correlations. '
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('simname', type=str,
                        help='The name of the simulation to populate. Defaults are stored in the myCats module.')
    parser.add_argument('outputdir', type = str,
                        help='The directory to store the outputs of the calculations. ')
    parser.add_argument('--plot', action = 'store_true',
                        help = 'Determine if plots should be saved along with the calculations.')

    args = parser.parse_args()

    allCorrFunc(args.simname, args.outputdir, args.plot)
