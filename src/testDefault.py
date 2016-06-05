#Having trouble with f_c changing the correlation function; I'm not sure if it's a real effect or not. Gonna calculate.
from allCorrFunc import corrFunc

simname = 'aardvark'
corrFunc(simname, 1.0, '/u/ki/swmclau2/des/HODOutput/%s_step_tests/'%simname,logMmin = 15.5, Lbox = 400.0, npart = 2048)
