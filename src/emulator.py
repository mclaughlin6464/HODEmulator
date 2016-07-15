#!bin/bash/
#Author: Sean McLaughlin
#This module utilizies George to construct an emulator for the correlation function.

from glob import glob
import os
import numpy as np
import george
from george.kernels import ExpSquaredKernel

def emulator(directory):
    '''Main emulator module.'''

    if not os.path.isdir(directory): 
        raise IOError('%s is not a directory!'%directory)

    files = sorted(glob(os.path.join(directory, 'corr_1.000_redMagic_mm_*.npy')))
    #cov_files = sort(glob(os.path.join(directory, 'cov_1.000_redMagic_mm_*.npy')))
    used_files = files[0::2]#only use half to train, for now. 
    #u_cov_files = cov_files[xrange(0,len(files), 2)]

    x,y,yerr = [],[],[]

    #for f,cf in zip(used_files, u_cov_files):
    for f in used_files:
        data = np.loadtxt(f)
        #cov = np.loadtxt(cf)
        basename = os.path.basename(f)
        mass = float(basename.split('_')[-1][:-4]) #mess of parsing to get the mass
        for r, xi in data.T:
            x.append([mass, r])
            y.append(xi)
        #yerr.extend(np.sqrt(np.diag(cov)))

    #x,y,yerr = np.array(x), np.array(y), np.array(yerr)
    x,y = np.array(x), np.array(y)

    metric = 1.0 #variance in each dimension of the kernel. Don't understand its use entirely, yet.
    kernel = ExpSquaredKernel(metric, ndim = 1)
    gp = george.GP(kernel)

    print x.shape

    gp.compute(x[:19,1])#,yerr)
    print gp.lnliklihood(y)

from sys import argv 
emulator(argv[1])
