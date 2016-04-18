#!bin/bash/
#@Author: Sean McLaughlin
#This file takes cached halo catalogs, populates them according to some HOD, and calculates the gg cross correlations

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set()

from itertools import cycle
from multiprocessing import cpu_count
from halotools.empirical_models import HodModelFactory, TrivialPhaseSpace, NFWPhaseSpace
from halotools.sim_manager import CachedHaloCatalog
from halotools.mock_observables import return_xyz_formatted_array, tpcf, tpcf_one_two_halo_decomp, wp
#TODO fix relative import
from redMagicHOD import RedMagicCens, RedMagicSats

#TODO argparse for simname
simname = 'fox'

outputdir = '/u/ki/swmclau2/des/HODOutput/%s/'%simname

scale_factors = [0.25,0.333,0.5,  0.540541, 0.588235, 0.645161, 0.714286, 0.8, 0.909091, 1.0 ] #sf of emu and fox

halocats = {}
models = {}
for sf in scale_factors:
    rz = 1.0/sf-1
    halocats[sf] = CachedHaloCatalog(simname = simname, halo_finder = 'rockstar',version_name = 'most_recent',redshift = rz)
    #halocats[sf] = CachedHaloCatalog(simname = 'multidark', halo_finder = 'rockstar',redshift = rz)


    models[sf] = HodModelFactory(
        centrals_occupation = RedMagicCens(redshift = rz),
        centrals_profile = TrivialPhaseSpace(redshift = rz),
        satellites_occupation = RedMagicSats(redshift = rz),
        satellites_profile = NFWPhaseSpace(redshift = rz))

#NOTE replace sf with a?
for sf in scale_factors:
    models[sf].populate_mock(halocats[sf], Num_ptcl_requirement = 30)

data = {}
pi_max = 40

rbins = np.logspace(-1, 1.25, 15)
rbin_centers = (rbins[1:]+rbins[:-1])/2

for sf in scale_factors:
    x, y, z = [models[sf].mock.galaxy_table[c] for c in ['x','y','z'] ]
    pos = return_xyz_formatted_array(x,y,z)

    xi_all = tpcf(pos, rbins, period = models[sf].mock.Lbox, num_threads =  cpu_count())

    halo_hostid = models[sf].mock.galaxy_table['halo_id']

    xi_1h, xi_2h = tpcf_one_two_halo_decomp(pos,
                    halo_hostid, rbins,
                    period = models[sf].mock.Lbox, num_threads =  cpu_count(),
                    max_sample_size = 1e7)


    wp_all = wp(pos, rbins, pi_max, period=models[sf].mock.Lbox, num_threads = cpu_count())

    data[sf] = (xi_all, xi_1h, xi_2h, wp_all)

colors = cycle(sns.color_palette())
fig = plt.figure(figsize = (15, 15))

for sf, color in zip(scale_factors, colors):
    plt.subplot(211)
    rz = 1.0/sf -1
    plt.plot(rbin_centers, data[sf][0],
             label='$z = %.2f$'%rz, color=color)
    plt.plot(rbin_centers, data[sf][1], ls = '--', color = color)
    plt.plot(rbin_centers, data[sf][2], ls = '-.', color = color)

    plt.subplot(221)
    plt.plot(rbin_centers, data[sf][3],
             label='$z = %.2f$'%rz,
             color= color )

plt.subplot(211)
plt.title('Cross Correlations')
plt.xlim(xmin = 0.1, xmax = 10)
plt.ylim(ymin = 1, ymax = 1e4)
plt.loglog()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel(r'$r $  $\rm{[Mpc]}$', fontsize=25)
plt.ylabel(r'$\xi_{\rm gg}(r)$', fontsize=25)
plt.legend(loc='best', fontsize=20)

plt.subplot(221)
plt.title('Projected Cross Correlations')
plt.xlim(xmin = 0.1, xmax = 10)
plt.ylim(ymin = 0.5, ymax = 5e3)
plt.loglog()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel(r'$r_{\rm p} $  $\rm{[Mpc]}$', fontsize=25)
plt.ylabel(r'$w_{\rm p}(r_{\rm p})$', fontsize=25)
plt.legend(loc='best', fontsize=20)

plt.savefig(outputdir+'emu_xi.png')

for sf in scale_factors:
    np.savetxt(outputdir+'xi_all_%.3f.npy'%sf, data[sf][0])
    np.savetxt(outputdir+'xi_1h_%.3f.npy'%sf, data[sf][1])
    np.savetxt(outputdir+'xi_2h_%.3f.npy'%sf, data[sf][2])
    np.savetxt(outputdir+'wp_all_%.3f.npy'%sf, data[sf][3])
