#!/bin/bash/
#@Author Sean McLaughlin
#There's been problems with the correlation function. The dummy test now is to compare the halo corr funcs
#I'm using two different functions; one is from halotools and the other is from Yao.

import numpy as np
import argparse
from multiprocessing import cpu_count

from halotools.sim_manager import CachedHaloCatalog
from halotools.mock_observables import return_xyz_formatted_array, tpcf

from helpers.CorrelationFunction import correlation3d

from myCats import *

def haloCorr(simname, scale_factor, outputdir,  **kwargs):
    'Calculate the cross correlation for a single catalog at a single scale factor'
    cat = cat_dict[simname](**kwargs) #TODO better handling of arguements
    _haloCorr(cat, scale_factor, outputdir)

def _haloCorr(cat, scale_factor, outputdir):
    'Helper function that uses the built in cat object'
    RBINS = np.logspace(-1, 1.25, 15)
    redshift = 1.0/scale_factor - 1.0

    #Note: Confusing name between cat and halocat. Consider changing.
    halocat = CachedHaloCatalog(simname = cat.simname, halo_finder = cat.halo_finder,version_name = cat.version_name, redshift = redshift)

    #Now, calculate with Halotools builtin
    #TODO include the fast version
    sample_mask = halocat.halo_table['halo_mvir'] > 7e12
    x, y, z = [halocat.halo_table[c] for c in ['halo_x','halo_y','halo_z'] ]
    pos = return_xyz_formatted_array(x,y,z, mask = sample_mask)
    #TODO N procs
    xi_all = tpcf(pos, RBINS, period = halocat.Lbox, num_threads =  cpu_count())

    np.savetxt(outputdir + 'xi_all_halo_%.3f_ht.npy' %(scale_factor), xi_all)

    xi_all = correlation3d(pos, RBINS, halocat.Lbox)

    np.savetxt(outputdir + 'xi_all_halo_%.3f_yao.npy' %(scale_factor), xi_all)

if __name__ == '__main__':
    desc = 'Calculate the correlation function for halos.  '
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('simname', type=str,
                        help='The name of the simulation to populate. Defaults are stored in the myCats module.')
    parser.add_argument('outputdir', type = str,
                        help='The directory to store the outputs of the calculations. ')

    args = parser.parse_args()
    #TODO remove kwargs
    haloCorr(args.simname, 1.0, args.outputdir, Lbox = 250.0, npart = 2560)



