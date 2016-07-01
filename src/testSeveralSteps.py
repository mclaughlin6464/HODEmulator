#I'm going to be testing the correlation functions of the step Function with several different thresholds.
from allCorrFunc import corrFunc
import argparse

def main(simname,outputdir, **kwargs):

    for m in [12.0 + i*0.25 for i in xrange(8)]:
        print 'Log Min Mass: %e'%m
        params = {'logMmin':m}
        corrFunc(simname, 1.0, outputdir+'/%s_step_tests/'%simname,params = params, **kwargs)

if __name__ == '__main__':
    desc = 'Run corrFunc on several mass steps.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('simname', type=str,
                        help='The name of the simulation to populate. Defaults are stored in the myCats module.')
    parser.add_argument('outputdir', type=str,
                        help='The directory to store the outputs of the calculations. Stored in <simname>_step_tests')
    parser.add_argument('--stepFunc', action='store_true',
                        help = 'Use a stepfunction HOD instead of the default RedMagic')
    parser.add_argument('--n_ptcl', type = int, default=0, help = 'Number of particles in a halo. Default is none.' )

    args = parser.parse_args()

    HOD = 'stepFunc' if args.stepFunc else 'redMagic'

    main(args.simname, args.outputdir, HOD = HOD, n_ptcl = args.n_ptcl)


