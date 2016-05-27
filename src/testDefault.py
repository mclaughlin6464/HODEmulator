#!/opt/rh/python27/root/usr/bin/
#Having trouble with f_c changing the correlation function; I'm not sure if it's a real effect or not. Gonna calculate.
from allCorrFunc import corrFunc

simname = 'multidark'
corrFunc(simname, 1.0, '/scratch/users/swmclau2/HODOutput/%s_tests/'%simname)
