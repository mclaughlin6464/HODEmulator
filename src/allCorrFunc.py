#!bin/bash/
#@Author: Sean McLaughlin
#This file takes cached halo catalogs, populates them according to some HOD, and calculates the gg cross correlations

import numpy as np
import argparse
from multiprocessing import cpu_count
from os.path import isdir, dirname, abspath, join
from halotools.empirical_models import HodModelFactory, TrivialPhaseSpace, NFWPhaseSpace
from halotools.empirical_models import Zheng07Cens, Zheng07Sats
from halotools.sim_manager import CachedHaloCatalog
from halotools.mock_observables import return_xyz_formatted_array,tpcf, tpcf_jackknife, tpcf_one_two_halo_decomp, wp
try:
    CORRFUNC = True
    from Corrfunc._countpairs import coutpairs_xi


except ImportError:
    CORRFUNC = False

from redMagicHOD import RedMagicCens, RedMagicSats, StepFuncCens, StepFuncSats
from myCats import *

N_PTCL = 0 
PI_MAX = 40

RBINS = np.logspace(-1, 1.25, 15)
RBIN_CENTERS = (RBINS[1:]+RBINS[:-1])/2 #just for plotting

#TODO change name so as to not overlap with CorrFunc
def corrFunc(simname, scale_factor, outputdir, HOD = 'redMagic', params = {}, n_ptcl = N_PTCL, rbins = RBINS, **kwargs):
    'Calculate the cross correlation for a single catalog at a single scale factor'

    if not isdir(outputdir):
        raise IOError("%s is not a directory."%outputdir)

    cat = cat_dict[simname](**kwargs)
    print str(cat)
    halocat, model = loadHaloAndModel(cat, HOD, scale_factor)
    data = popAndCorr(halocat,model, cat, params, n_ptcl, rbins)

    np.savetxt(outputdir + 'corr_%.3f_%s_mm_%.2f.npy' % (HOD, scale_factor, params['logMmin']), data)


def allCorrFunc(simname, outputdir,HOD = 'redmagic', params = {}, n_ptcl = N_PTCL, rbins = RBINS, **kwargs):
    'Calculates cross correlations for all scale factors cached for one halocatalog'
    if not isdir(outputdir):
        raise IOError("%s is not a directory."%outputdir)

    cat = cat_dict[simname](**kwargs)
    print str(cat)
    for a in cat.scale_factors:
        halocat, model = loadHaloAndModel(cat,HOD, a)
        data = popAndCorr(halocat, model, cat, params, n_ptcl, rbins)

        np.savetxt(outputdir + 'corr_%.3f_%s_mm_%.2f.npy' %(HOD, a, params['logMmin']), data)


def loadHaloAndModel(cat, HOD, scale_factor):
    '''Return the cached halo catalog and the appropriate HOD model'''
    try:
        idx = cat.scale_factors.index(scale_factor)
    except:
        print 'Provided scale_factor %.3f not cached for %s.'%(scale_factor, cat.simname)
        raise

    if HOD == 'redMagic':
        cens_occ = RedMagicCens(redshift=cat.redshifts[idx])
        sats_occ = RedMagicSats(redshift=cat.redshifts[idx],modulate_with_cenocc = True)
    elif HOD == 'stepFunc':
        cens_occ = StepFuncCens(redshift=cat.redshifts[idx])
        sats_occ = StepFuncSats(redshift=cat.redshifts[idx])
    else:
        raise ValueError('%s invalid input for HOD'%HOD)

    #Note: Confusing name between cat and halocat. Consider changing.
    halocat = CachedHaloCatalog(simname = cat.simname, halo_finder = cat.halo_finder,version_name = cat.version_name, redshift = cat.redshifts[idx])

    model = HodModelFactory(
        centrals_occupation=cens_occ,
        centrals_profile=TrivialPhaseSpace(redshift=cat.redshifts[idx]),
        satellites_occupation=sats_occ,
        satellites_profile=NFWPhaseSpace(redshift=cat.redshifts[idx]))

    return halocat, model

def popAndCorr(halocat,model, cat, params = {}, n_ptcl = N_PTCL, rbins = RBINS):
    '''Populate a halocat with a model and calculate the tpcf, tpcf_1h, tpcf_2h, and projected corr fun'''
    print 'Min Num Particles: %d\t%d bins'%(n_ptcl, len(rbins))
    model.update(params)#insert new params into model
    # Note: slow
    model.populate_mock(halocat,Num_ptcl_requirement=n_ptcl)

    #Now, calculate with Halotools builtin
    x, y, z = [model.mock.galaxy_table[c] for c in ['x','y','z'] ]
    #mask = model.mock.galaxy_table['halo_mvir'] < 1e15/cat.h
    pos = return_xyz_formatted_array(x,y,z)#, mask = mask)

    #TODO N procs
    if CORRFUNC:
        #write bins to file
        BINDIR = dirname(abspath(__file__))  # location of files with bin edges
        with open(join(BINDIR , '/binfile'), 'w') as f:
            for low, high in zip(RBINS[:-1], RBINS[1:]):
                f.write('\t%f\t%f\n' % (low, high))

        xi_all = coutpairs_xi(model.mock.Lbox*cat.h, cpu_count(), join(BINDIR,'/binfile'), x*cat.h, y*cat.h, z*cat.h )
        xi_all = np.array(xi_all)

    else:

        xi_all = tpcf(pos*cat.h, RBINS, period = model.mock.Lbox*cat.h, num_threads =  cpu_count())

    #TODO ways to decide which of these to call
    #randoms = np.random.random(pos.shape)*model.mock.Lbox*cat.h
    #xi_all, xi_cov = tpcf_jackknife(pos*cat.h,randoms, RBINS, period = model.mock.Lbox*cat.h, num_threads =  cpu_count())

    halo_hostid = model.mock.galaxy_table['halo_id']

    xi_1h, xi_2h = tpcf_one_two_halo_decomp(pos*cat.h,
                                            halo_hostid, rbins,
                                            period=model.mock.Lbox*cat.h, num_threads=cpu_count(),
                                            max_sample_size=1e7)

    #wp_all = wp(pos*cat.h, RBINS, PI_MAX, period=model.mock.Lbox*cat.h, num_threads = cpu_count())

    output = np.stack([RBIN_CENTERS, xi_all, xi_1h, xi_2h])

    return output


    #np.savetxt(outputdir + 'corr_%.3f_default_mm_%.2f.npy' %(scale_factor, logMmin), output)
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
