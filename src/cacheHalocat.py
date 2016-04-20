#!/bin/bash
#@Author: Sean McLaughlin
#This module caches a halocat at a given location.
#TODO make this handle arguements?
#Don't know if that's something I'd want this to be capable of.

from halotools.sim_manager import RockstarHlistReader

simname = 'emu'

if simname == 'emu':
    #TODO specify for which box since we have 40
    #loc = '/nfs/slac/g/ki/ki22/cosmo/beckermr/tinkers_emu/Box000/halos/m200b/'
    loc =  '/u/ki/swmclau2/des/emu/Box000/'

    columns_to_keep = {'halo_id': (0, 'i8'), 'halo_upid': (41, 'i8'), 'halo_x': (8, 'f4'), 'halo_y': (9, 'f4'),
                       'halo_z': (10, 'f4')
        , 'halo_vx': (11, 'f4'), 'halo_vy': (12, 'f4'), 'halo_vz': (13, 'f4'), 'halo_mvir': (2, 'f4'),
                       'halo_rvir': (5, 'f4'), 'halo_rs': (6, 'f4')}  # what else?

    halo_finder = 'rockstar'
    version_name = 'most_recent'
    Lbox = 1050.0
    pmass = 3.9876e10

    fnames = ['out_%d.list' % i for i in xrange(10)]  # redshift = 0
    scale_a = [0.25, 0.333, 0.5, 0.540541, 0.588235, 0.645161, 0.714286, 0.8, 0.909091, 1.0]

elif simname == 'fox':
    loc = '/nfs/slac/g/ki/ki23/des/BCCSims/Fox/Lb400/halos/rockstar/output/hlists/'
    columns_to_keep = {'halo_id': (1, 'i8'), 'halo_upid': (6, 'i8'), 'halo_x': (17, 'f4'), 'halo_y': (18, 'f4'),
                       'halo_z': (19, 'f4'), 'halo_vx': (20, 'f4'), 'halo_vy': (21, 'f4'), 'halo_vz': (22, 'f4'),
                       'halo_mvir': (10, 'f4'), 'halo_rvir': (11, 'f4'), 'halo_rs': (12, 'f4')}

    halo_finder = 'rockstar'
    version_name = 'most_recent'
    Lbox = 400.0
    pmass = 6.58298e8

    fnames = ['hlist_%d' % n for n in [46, 57, 73, 76, 79, 82, 86, 90, 95, 99]]
    scale_a = [0.25, 0.333, 0.5, 0.540541, 0.588235, 0.645161, 0.714286, 0.8, 0.909091, 1.0]  # not exact, but close.

elif simname == 'multidark_highres':
    loc = '/nfs/slac/g/ki/ki20/cosmo/behroozi/MultiDark/hlists/'
    columns_to_keep = {'halo_id': (1, 'i8'), 'halo_upid': (6, 'i8'), 'halo_x': (17, 'f4'), 'halo_y': (18, 'f4'),
                       'halo_z': (19, 'f4'), 'halo_vx': (20, 'f4'), 'halo_vy': (21, 'f4'), 'halo_vz': (22, 'f4'),
                       'halo_mvir': (10, 'f4'), 'halo_rvir': (11, 'f4'), 'halo_rs': (12, 'f4')}

    halo_finder = 'rockstar'
    version_name = 'most_recent'
    Lbox = 1e3
    pmass = 8.721e9  # foudn on MD site
    # Trying to be close to the ones i have for the other
    scale_a = [0.25690, 0.34800, 0.49990, 0.53030, 0.65180, 0.71250, 0.80370, 0.91000, 1.00110]
    fnames = ['hlist_%.5f.list' % sf for sf in scale_a]  # redshift = 0

#TODO else raise error

for fn, a in zip(fnames, scale_a):
    z = 1.0/a-1
    out = '/u/ki/swmclau2/des/halocats/hlist_%.2f.list.%s.hdf5'%(a, simname)
    reader = RockstarHlistReader(loc+fn, columns_to_keep, out, simname, halo_finder, z,
                           version_name,Lbox, pmass, overwrite = True)
    reader.read_halocat(["halo_rvir", "halo_rs"], write_to_disk = True, update_cache_log = True)