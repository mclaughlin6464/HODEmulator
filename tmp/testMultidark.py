#@Author Sean McLaughlin
#Tmp helper script i'm writing to correlate the multidark with a slight tweak.

from os import mkdir, path
from allCrossCorr import crossCorr

simname = 'bolshoi'
scale_factor = 1.0
outputdir = '/u/ki/swmclau2/des/HODOutput/%s/'%simname

crossCorr(simname, scale_factor, outputdir)
