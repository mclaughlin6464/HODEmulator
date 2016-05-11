#Having trouble with f_c changing the correlation function; I'm not sure if it's a real effect or not. Gonna calculate.
from ..src.allCorrFunc import corrFunc

simname = 'chinchilla'
corrFunc(simname, 1.0, '/u/ki/swmclau2/des/HODOutput/%s_tests/'%simname, Lbox = 125.0, npart = 2048)
