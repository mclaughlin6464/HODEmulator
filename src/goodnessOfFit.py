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
    #train_emulator(gp, y)

    #need to go to truth_dir
    #TODO pass in a different Truth dir, possibly by modifying kwargs
    if 'directory' in kwargs:
        del kwargs['directory']

    corr_files = sorted(glob(path.join(truth_dir, '*corr*.npy')))
    cov_files = sorted(glob(path.join(truth_dir, '*cov*.npy')))

    np.random.seed(int(time()))

    idxs = np.random.choice(len(corr_files), N, replace = False)

    errors = []

    for idx in idxs:
        params, r_centers, true_xi ,_ = file_reader(corr_files[idx], cov_files[idx])
        pred_log_xi, _ = emulate_wrt_r(gp, em_y, params, r_centers)

        RMSE = np.sqrt(np.mean((pred_log_xi-np.log10(true_xi))**2))
        NRMSE = RMSE/true_xi.mean()

        errors.append(NRMSE)

    return np.array(errors)

if __name__ == '__main__':
    from sys import argv
    N = None
    if len(argv) > 1:
        N = int(argv[1])
    errors = root_mean_square(N)
    print errors.mean(), errors.std()

