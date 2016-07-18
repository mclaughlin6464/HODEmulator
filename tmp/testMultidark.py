#@Author Sean McLaughlin
#Tmp helper script i'm writing to correlate the multidark with a slight tweak.

from os import mkdir, path
from ..src.allCorrFunc import corrFunc

simname = 'bolshoi'
scale_factor = 1.0
outputdir = '/u/ki/swmclau2/des/HODOutput/%s/'%simname

corrFunc(simname, scale_factor, outputdir)
