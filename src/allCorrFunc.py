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
from redMagicHOD import RedMagicCens, RedMagicSats, StepFuncCens, StepFuncSats
from myCats import *

N_PTCL = 0 
PI_MAX = 40

RBINS = np.logspace(-1, 1.7, 20)
#RBIN_CENTERS = (RBINS[1:]+RBINS[:-1])/2 #just for plotting

#TODO will need ways to pass params into the model when populating. Could just use kwargs, but how to separate cat kwargs?
#Could pass in a dict and update the param dict
def corrFunc(simname, scale_factor, outputdir, plot = False,mMin = 7e12,  **kwargs):
    'Calculate the cross correlation for a single catalog at a single scale factor'
    cat = cat_dict[simname](**kwargs) #TODO better handling of arguements
    _corrFunc(cat, scale_factor, outputdir, plot, mMin = mMin)

def allCorrFunc(simname, outputdir, plot = False, **kwargs):
    'Calculates cross correlations for all scale factors cached for one halocatalog'
    cat = cat_dict[simname](**kwargs) #TODO better handling or arguements
    for a in cat.scale_factors:
        _corrFunc(cat, a, outputdir, plot)

def _corrFunc(cat, scale_factor, outputdir, plot = False, mMin = 7e12):
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
        #centrals_occupation=RedMagicCens(redshift=cat.redshifts[idx]),
        centrals_occupation=StepFuncCens(redshift=cat.redshifts[idx], mMin=mMin, h = cat.h),
        centrals_profile=TrivialPhaseSpace(redshift=cat.redshifts[idx]),
        #satellites_occupation=RedMagicSats(redshift=cat.redshifts[idx]),
        satellites_occupation=StepFuncSats(redshift=cat.redshifts[idx]),
        satellites_profile=NFWPhaseSpace(redshift=cat.redshifts[idx]))

    #model.param_dict['logMmin'] = 13.1/cat.h

    #Note: slow
    model.populate_mock(halocat, Num_ptcl_requirement = N_PTCL) #TODO try again with 300 or a larger number for more robustness

    #Now, calculate with Halotools builtin
    #TODO include the fast version
    x, y, z = [model.mock.galaxy_table[c] for c in ['x','y','z'] ]
    #mask = model.mock.galaxy_table['halo_mvir'] < 1e15/cat.h
    pos = return_xyz_formatted_array(x,y,z)#, mask = mask)
    #TODO N procs
    xi_all = tpcf(pos*cat.h, RBINS, period = model.mock.Lbox*cat.h, num_threads =  cpu_count())

    #xi_all, xi_cov = tpcf_jackknife(pos*cat.h, RBINS, period = model.mock.Lbox*cat.h, num_threads =  cpu_count())

    #halo_hostid = model.mock.galaxy_table['halo_id']
    '''
    xi_1h, xi_2h = tpcf_one_two_halo_decomp(pos*cat.h,
                    halo_hostid, RBINS,
                    period = cat.h*model.mock.Lbox, num_threads =  cpu_count(),
                    max_sample_size = 1e7)
    '''
    #wp_all = wp(pos*cat.h, RBINS, PI_MAX, period=model.mock.Lbox*cat.h, num_threads = cpu_count())

    np.savetxt(outputdir + 'xi_all_%.3f_stepFunc_400_2048_mm_%2f.npy' %(scale_factor, mMin), xi_all)
    #np.savetxt(outputdir + 'xi_cov_%.3f_default_125_2048.npy' %(scale_factor), xi_cov)

    #np.savetxt(outputdir + 'xi_1h_%.3f_stepFunc.npy' %(scale_factor), xi_1h)
    #np.savetxt(outputdir + 'xi_2h_%.3f_stepFunc.npy' %(scale_factor), xi_2h)
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
