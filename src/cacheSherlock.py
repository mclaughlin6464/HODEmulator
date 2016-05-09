#!/bin/bash/
#@Author Sean McLaughlin
#This script will cache the hlists I've moved to Sherlock.
from cacheHalocat import cacheHalocat

sims = {}
sims['emu'] = []#kwargs for each sim
sims['multidark_highres'] = []
sims['fox'] = []

boxsize_npart = [(125.0, 1024),(125.0, 2048), (250.0, 1024), (250.0, 2048),
                (250.0, 2560), (250.0, 512), (250.0, 768), (400.0, 1024),
                (400.0, 2048), (400.0, 768)]

sims['chinchilla'] = [{'Lbox':lb, 'npart':np} for lb,np in boxsize_npart ]

scale_factor = 1.0

for simname, kwargs in sims.iteritems():
    print simname
    try:
        if len(kwargs) == 0:
            cacheHalocat(simname)
        else:
            for kw in kwargs:
		#only do z = 0 for cinchillas
                print kw
                cacheHalocat(simname, scale_factors=[scale_factor], **kw)
        print

    except:
        print 'An error occured for %s'%simname
        #continue
        raise
