#!/bin/bash/
#@Author Sean McLaughlin
#There's been problems with the correlation function. I'm comparing the halo corr func to the gal corr func 
#I'm using two different functions; one is from halotools and the other is from Yao.

import numpy as np
import argparse
from multiprocessing import cpu_count

from halotools.sim_manager import CachedHaloCatalog
from halotools.mock_observables import return_xyz_formatted_array, tpcf
from halotools.empirical_models import HodModelFactory, TrivialPhaseSpace, NFWPhaseSpace
from halotools.empirical_models import Zheng07Cens, Zheng07Sats

from redMagicHOD import *

from helpers.CorrelationFunction import correlation3d

from myCats import *

def galCorr(simname, scale_factor, outputdir,  **kwargs):
    'Calculate the cross correlation for a single catalog at a single scale factor'
    cat = cat_dict[simname](**kwargs) #TODO better handling of arguements
    _galCorr(cat, scale_factor, outputdir)

def _galCorr(cat, scale_factor, outputdir):
    'Helper function that uses the built in cat object'
    h = 0.7
    RBINS = np.logspace(-1, 1.25, 15)
    redshift = 1.0/scale_factor - 1.0

    if outputdir[-1] != '/':
        outputdir+='/'

    #Note: Confusing name between cat and halocat. Consider changing.
    halocat = CachedHaloCatalog(simname = cat.simname, halo_finder = cat.halo_finder,version_name = cat.version_name, redshift = redshift)

    model = HodModelFactory(
            centrals_occupation=RedMagicCens(redshift=redshift),
            centrals_profile=TrivialPhaseSpace(redshift=redshift),
            satellites_occupation=RedMagicSats(redshift=redshift),
            satellites_profile=NFWPhaseSpace(redshift=redshift))

    model.populate_mock(halocat) #default NPTCL

    #Now, calculate with Halotools builtin
    #TODO include the fast version
    x, y, z = [model.mock.galaxy_table[c]*h for c in ['x','y','z'] ]
    pos = return_xyz_formatted_array(x,y,z)
    #TODO N procs
    xi_all = tpcf(pos, RBINS, period = model.mock.Lbox*h, num_threads =  cpu_count())

    np.savetxt(outputdir + 'xi_all_gal_%.3f_ht_h.npy' %(scale_factor), xi_all)

    #xi_all = correlation3d(pos, RBINS, model.mock.Lbox)

    #np.savetxt(outputdir + 'xi_all_gal_%.3f_yao.npy' %(scale_factor), xi_all)

if __name__ == '__main__':
    desc = 'Calculate the correlation function for galaxies.  '
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('simname', type=str,
                        help='The name of the simulation to populate. Defaults are stored in the myCats module.')
    parser.add_argument('outputdir', type = str,
                        help='The directory to store the outputs of the calculations. ')

    args = parser.parse_args()
    #TODO remove kwargs
    galCorr(args.simname, 1.0, args.outputdir, Lbox = 250.0, npart = 2560)



