#I'm going to be testing the correlation functions of the step Function with several different thresholds.
from allCorrFunc import corrFunc

simname = 'chinchilla'
for m in [12.25 + i*0.25 for i in xrange(7)]:
    print 'Min Mass: %e'%m
    corrFunc(simname, 1.0, '/u/ki/swmclau2/des/HODOutput/%s_step_tests/'%simname,mMin = 10**m, Lbox = 400.0, npart = 2048)
