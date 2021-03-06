#@Author Sean McLaughlin
#Tmp helper script i'm writing to cache and cross-correlate all the chinchillas there are.

from os import mkdir, path
from cacheHalocat import cacheHalocat
from allCorrFunc import corrFunc

simname = 'chinchilla'
boxsize_npart = [(125.0, 1024),(125.0, 2048), (250.0, 1024), (250.0, 2048),
                (250.0, 2560), (250.0, 512), (250.0, 768), (400.0, 1024),
                (400.0, 2048), (400.0, 512), (400.0, 768)]
#boxsize_npart = [(250.0, 1024), (250.0, 2048),(250.0, 2560), (250.0, 512), (250.0, 768)]

scale_factor = 0.658
outputdir = '/u/ki/swmclau2/des/HODOutput/chinchilla/'

for boxsize, npart in boxsize_npart: #see what I did there?
    try:
        print boxsize, npart
        cacheHalocat(simname, Lbox = boxsize, npart = npart, scale_factors = [scale_factor])

    except:
        print 'An error occured for %.2f, %d'%(boxsize, npart)
        #continue
        raise
