#I'm going to be testing the correlation functions of the step Function with several different thresholds.
from allCorrFunc import corrFunc

simname = 'emu'
for m in [12.1 + i*0.2 for i in xrange(6)]:
    print 'Log Min Mass: %e'%m
    corrFunc(simname, 1.0, '/u/ki/swmclau2/des/HODOutput/%s_step_tests/'%simname,logMmin = m)
