# This module is similar to testSeveralSteps, but with an increase in scale.
#@Author Sean McLaughlin
import numpy as np
from time import time
from os import path
from itertools import izip
import argparse

from myCats import cat_dict
from allCorrFunc import loadHaloAndModel, popAndCorr
from doBatchCalls import BOUNDS #i Need them in both places but it's smarter to ahve ti there.
# Given by Elisabeth from on high


# TODO not hardcoding some of these? Depends on my use i guess.
# Will have to see how i end up using this.
SIMNAME = 'chinchilla'  # hardcode for noew
REDSHIFT = 0.0#0.5
N_PTCL = 200
RBINS = np.logspace(-1, 1.7, 20)
RBIN_CENTERS = (RBINS[1:] + RBINS[:-1]) / 2

def paramCube(outputdir, fixed_params={}, n_per_dim=5):
    if type(n_per_dim) is int:
        n_per_dim = {key: n_per_dim for key in BOUNDS.iterkeys()}
    elif type(n_per_dim) is not dict:
        pass  # WOOPWOOP ERROR
        # TODO error here

    #TODO if user passed in fixed_params and n_per_dim but they have different keys!

    values = {}
    for param in BOUNDS.iterkeys():
        if param in fixed_params:
            n_per_dim[param] = 1
            values[param] = np.array([fixed_params[param]])
        else:  # param in varied_params
            values[param] = np.linspace(BOUNDS[param][0], BOUNDS[param][1], num=n_per_dim[param])

    n_total = np.prod(n_per_dim.values())
    points = [{} for i in xrange(n_total)]
    fixed_base = '_'.join('%s%.2f' % (key, val) for key, val in fixed_params.iteritems()) + '_'
    outbase = [fixed_base for i in xrange(n_total)]

    n_segment = n_total  # not necessary, but notaionally clearer
    for param in sorted(BOUNDS.iterkeys()):  # sorted to make deterministic, though it may already be.
        n_segment /= n_per_dim[param]
        for i, p in enumerate(points):
            idx = (i / n_segment) % n_per_dim[param]
            p[param] = values[param][idx]
            outbase[i] += str(idx)  # now each outbase has a unique combination of indexes
    # now each dictionary in values carries a unique combination of parameters for the emulator
    # if memory was an issue one could just run the model at each step instead of generating them all.
    # i don't think 1000 dictionaries is the worst of my memory issues.

    # now, send each fo these to my code.
    for p, out in izip(points, outbase):
        calc_galaxy_autocorr(SIMNAME, 1 / (1 + REDSHIFT), path.join(outputdir, out), params=p, Lbox=400, npart=2048)


# mostly copied from allCorrFunc. I don't wanan break backwards compatibaility yet
# but I need to make some changes here.
def calc_galaxy_autocorr(simname, scale_factor, outbase, params={}, **kwargs):
    'Calculate the cross correlation for a single catalog at a single scale factor'
    t0 = time()

    cat = cat_dict[simname](**kwargs)
    print str(cat)
    halocat, model = loadHaloAndModel(cat, 'redMagic', scale_factor)
    data, cov = popAndCorr(halocat, model, cat, params, N_PTCL, RBINS)
    header_start = ['Cosmology: %s'%simname, 'Params for HOD:' ]
    header_start.extend('%s:%.3f'%(key,val) for key, val in params.iteritems())
    header = '\n'.join(header_start)
    np.savetxt(outbase + '_corr_%.3f.npy' % (scale_factor), data,
            header = header)
    np.savetxt(outbase + '_cov_%.3f.npy' % (scale_factor), cov,
            header = header)

    print '\nTotal Time: %.3f\n' % (time() - t0)


if __name__ == '__main__':
    desc = 'Run my correlation function calculator over a large range of parameters'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('outputdir', type=str,
                        help='The directory to store the outputs of the calculations.')

    for param in BOUNDS.iterkeys():
        parser.add_argument(''.join(['--', param])) #no help scripts #YOLO

    args = vars(parser.parse_args())

    outputdir = args['outputdir']
    del args['outputdir']
    for key in args.keys():
        if args[key] is not None:
            args[key] = float(args[key])
        else:
            del args[key]

    #pretty smart if i say so myself
    #leave default nperdim for now..
    print args
    paramCube(outputdir, fixed_params=args)
