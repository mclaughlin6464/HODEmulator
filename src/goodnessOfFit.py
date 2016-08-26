# /bin/bash
# @Author Sean McLaughlin
# A collection of functions that compare "true" calculations to a trained emulator.
from os import path
from time import time
from glob import glob
import numpy as np

from allCorrFunc import RBINS
from emulator import build_emulator,file_reader,\
                        train_emulator, emulate_wrt_r

TRUTH_DIR = '/u/ki/swmclau2/des/EmulatorData/'
NBINS = len(RBINS)-1

def root_mean_square(N = None,truth_dir=TRUTH_DIR, **kwargs):
    gp,em_y, _ = build_emulator(**kwargs)
    print 'Emulator built!'
    #train_emulator(gp, y)
    print 'Emulator trained!'

    corr_files = sorted(glob(path.join(truth_dir, '*corr*.npy')))
    cov_files = sorted(glob(path.join(truth_dir, '*cov*.npy')))

    np.random.seed(int(time()))

    if N is None:
        idxs = np.arange(len(corr_files))
    else:
        idxs = np.random.choice(len(corr_files), N, replace = False)

    errors = []
    print 'Starting...'
    t0 = time()
    for idx in idxs:
        if idx%100 == 0:
            print 'Index: %d\tTime Elapsed:%.2f s'%(idx, time()-t0)
        params, r_centers, true_xi ,_ = file_reader(corr_files[idx], cov_files[idx])
        pred_log_xi, _ = emulate_wrt_r(gp, em_y, params, np.log10(r_centers) )

        #print ((pred_log_xi-np.log10(true_xi))**2)
        #print (np.log10(true_xi)**2) 
        #print ((pred_log_xi-np.log10(true_xi))**2)/(np.log10(true_xi)**2) 
        #print np.mean(((pred_log_xi-np.log10(true_xi))**2)/(np.log10(true_xi)**2)) 

        RMSE = np.sqrt(np.mean(((pred_log_xi-np.log10(true_xi))**2)/(np.log10(true_xi)**2)))

        errors.append(RMSE)

    return np.array(errors)

def calcR2(N = None,truth_dir=TRUTH_DIR, **kwargs):
    gp,em_y, _ = build_emulator(**kwargs)
    print 'Emulator built!'
    #train_emulator(gp, y)
    print 'Emulator trained!'

    corr_files = sorted(glob(path.join(truth_dir, '*corr*.npy')))
    cov_files = sorted(glob(path.join(truth_dir, '*cov*.npy')))

    np.random.seed(int(time()))

    if N is None:
        idxs = np.arange(len(corr_files))
    else:
        idxs = np.random.choice(len(corr_files), N, replace = False)

    r2_vals= []
    print 'Starting...'
    t0 = time()
    for idx in idxs:
        if idx%100 == 0:
            print 'Index: %d\tTime Elapsed:%.2f s'%(idx, time()-t0)
        params, r_centers, true_xi ,_ = file_reader(corr_files[idx], cov_files[idx])
        pred_log_xi, _ = emulate_wrt_r(gp, em_y, params, np.log10(r_centers) )

        SSR = np.sum((pred_log_xi-np.log10(true_xi))**2)
        SST = np.sum((np.log10(true_xi)-np.log10(true_xi).mean())**2)

        r2_vals.append(1-SSR/SST)

    return np.array(r2_vals)


if __name__ == '__main__':
    from sys import argv
    N = None
    if len(argv) > 1:
        N = int(argv[1])
    errors = root_mean_square(N)*100
    print errors.mean(), errors.std()
    print errors.max(), errors.min()
    print 
    print np.nanmean(errors), np.nanstd(errors)
    print np.nanmax(errors), np.nanmin(errors)

