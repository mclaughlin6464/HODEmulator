#!/bin/bash
#@Author: Sean McLaughlin
#This module caches a halocat at a given location.
#Don't know if that's something I'd want this to be capable of.

import argparse
from halotools.sim_manager import RockstarHlistReader
from myCats import *

def cacheHalocat(simname,Lbox,  **kwargs):

    cat = cat_dict[simname](Lbox,**kwargs) #TODO better handling of arguements

    for i in xrange(len(cat)): #Pythonic iteration?
        reader = RockstarHlistReader(cat.filenames[i], cat.columns_to_keep, cat.cache_locs[i], cat.simname,
                                     cat.halo_finder, cat.redshifts[i],
                                    cat.version_name,cat.Lbox, cat.pmass,
                                     overwrite = True)
        reader.read_halocat(cat.columns_to_convert, write_to_disk = True, update_cache_log = True)

if __name__ == '__main__':
    desc = 'Cache a particular simulation. Defaults are in the myCats module.'
    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument('simname', metavar='simname', type=str,
                        help='The name of the simulation to cache. Defaults are stored in the myCats module.')
    parser.add_argument('Lbox', type=float,
                        help='Size of the box in Mpc. Necessaryto identify a unique box.')
    #TODO do I want to have an advanced CLI? Connect to kwargs at all?
    parser.add_argument('--scale_factor', type = float, help = 'Scale factor to cache. Default is all. ')
    args = parser.parse_args()
    if args.scale_factor is None:
        cacheHalocat(args.simname, args.Lbox)
    else:
        cacheHalocat(args.simname, args.Lbox, scale_factors = [args.scale_factor])
