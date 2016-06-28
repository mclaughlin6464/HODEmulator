#!/opt/rh/python27/root/usr/bin/
#Having trouble with f_c changing the correlation function; I'm not sure if it's a real effect or not. Gonna calculate.
from ..src.allCorrFunc import corrFunc

simname = 'chinchilla'
corrFunc(simname, 1.0, '/home/swmclau2/scratch/HODOutput/%s_tests/'%simname,logMmin = 12.5, Lbox = 400.0, npart = 2048)
