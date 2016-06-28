#Having trouble with f_c changing the correlation function; I'm not sure if it's a real effect or not. Gonna calculate.
#TODO stil realtive import problems
from ..src.allCorrFunc import corrFunc

simname = 'chinchilla'
for f_c in [.2, .4, .6, .8, 1.0]:
    corrFunc(simname, 1.0, '/u/ki/swmclau2/des/HODOutput/%s/'%simname, f_c = f_c, Lbox = 250.0, npart = 2560)
