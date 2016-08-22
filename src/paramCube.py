# This module is similar to testSeveralSteps, but with an increase in scale.
#@Author Sean McLaughlin
import numpy as np
from time import time
from os import path
from itertools import izip
import argparse

from myCats import cat_dict
from allCorrFunc import loadHaloAndModel, popAndCorr, RBINS,  MIN_PTCL
from doBatchCalls import BOUNDS #i Need them in both places but it's smarter to ahve ti there.
from guppy import hpy
h = hpy()

# TODO not hardcoding some of these? Depends on my use i guess.
# Will have to see how i end up using this.
SIMNAME = 'chinchilla'  # hardcode for noew

REDSHIFT = 0.5#0.0
#N_PTCL = 200

RBIN_CENTERS = (RBINS[1:] + RBINS[:-1]) / 2

def paramCube(outputdir, fixed_params={}, n_per_dim=4, id_no=None):
    if type(n_per_dim) is int:
        n_per_dim = {key: n_per_dim for key in BOUNDS.iterkeys()}

    assert type(n_per_dim) is dict

    values = {}
    for param in BOUNDS.iterkeys():
        if param in fixed_params:
            n_per_dim[param] = 1
            values[param] = np.array([fixed_params[param]])
        else:  # param in varied_params
            values[param] = np.linspace(BOUNDS[param][0], BOUNDS[param][1], num=n_per_dim[param])

    n_total = np.prod(n_per_dim.values())
    if n_total == 1: #only one, we can skip all this stuff.
        calc_galaxy_autocorr(SIMNAME, 1 / (1 + REDSHIFT),
                             path.join(outputdir,'Emulator_lhc_'+ '%03d'%id_no if id_no is not None else 'Emulator'),
                             params=fixed_params, do_jackknife=True, Lbox=400, npart=2048)
        return

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
        #check if the file already exists
        if path.isfile(path.join(outputdir,out+'_corr_%.3f.npy'%(1/(1+REDSHIFT)) )):
            continue #file already exists!
            print 'Skipping %s'%out

        calc_galaxy_autocorr(SIMNAME, 1 / (1 + REDSHIFT), path.join(outputdir, out),
                             params=p,do_jackknife=True, Lbox=400, npart=2048)


# mostly copied from allCorrFunc. I don't wanan break backwards compatibaility yet
# but I need to make some changes here.

def calc_galaxy_autocorr(simname, scale_factor, outbase, params={},do_jackknife=True, **kwargs):
    'Calculate the cross correlation for a single catalog at a single scale factor'
    print h.heap()
    print '--'*25 
    t0 = time()

    cat = cat_dict[simname](**kwargs)
    print str(cat)
    halocat, model = loadHaloAndModel(cat, 'redMagic', scale_factor)
    if do_jackknife:
        data, cov = popAndCorr(halocat, model, cat, params,do_jackknife, MIN_PTCL, RBINS)
    else:
        data = popAndCorr(halocat, model, cat, params,do_jackknife, MIN_PTCL, RBINS)

    header_start = ['Cosmology: %s'%simname, 'Params for HOD:' ]
    header_start.extend('%s:%.3f'%(key,val) for key, val in params.iteritems())
    header = '\n'.join(header_start)
    np.savetxt(outbase + '_corr_%.3f.npy' % (scale_factor), data,
            header = header)
    if do_jackknife:
        np.savetxt(outbase + '_cov_%.3f.npy' % (scale_factor), cov,
            header = header)

    print '\nTotal Time: %.3f\n' % (time() - t0)

def testCube(outputdir, fixed_params={}, n_per_dim=4):
    '''Create fake data of the same structure as paramCube for testing. '''
    if type(n_per_dim) is int:
        n_per_dim = {key: n_per_dim for key in BOUNDS.iterkeys()}

    assert type(n_per_dim) is dict

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
    simname = SIMNAME
    scale_factor = 1/(1+REDSHIFT)

    for p, out in izip(points, outbase):

        ob = path.join(outputdir, out)

        #I could maybe do something more interesting than rands.

        data = np.stack( [(RBINS[1:] + RBINS[:-1]) / 2, np.random.rand(len(RBINS)-1)] )
        #cov = np.random.rand((len(RBINS), len(RBINS))) #could just do an eye matrix too.
        cov = np.eye((len(RBINS)-1))*np.random.rand()

        header_start = ['Cosmology: %s'%simname, 'Params for HOD:' ]
        header_start.extend('%s:%.3f'%(key,val) for key, val in p.iteritems())
        header = '\n'.join(header_start)
        np.savetxt(ob + '_corr_test_%.3f.npy' % (scale_factor), data,
                header = header)
        np.savetxt(ob + '_cov_test_%.3f.npy' % (scale_factor), cov,
                header = header)

if __name__ == '__main__':
    desc = 'Run my correlation function calculator over a large range of parameters'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('outputdir', type=str,
                        help='The directory to store the outputs of the calculations.')
    parser.add_argument('--test', action='store_true', help='Create fake data with a similar structure for testing.')
    parser.add_argument('--id',type=int,default=None, help='The job id for this call.')

    for param in BOUNDS.iterkeys():
        parser.add_argument(''.join(['--', param])) #no help scripts #YOLO

    args = vars(parser.parse_args())

    test = args['test']
    del args['test']
    outputdir = args['outputdir']
    del args['outputdir']
    id_no = args['id']
    del args['id']

    for key in args.keys():
        if args[key] is not None:
            args[key] = float(args[key])
        else:
            del args[key]

    #pretty smart if i say so myself
    #leave default nperdim for now..
    print args
    if not test:
        paramCube(outputdir, fixed_params=args, id_no=id_no)
    else:
        testCube(outputdir, fixed_params=args)
