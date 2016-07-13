#@Author Sean McLaughlin
#Send a collection of calls to paramCube to ki-ls
import numpy as np
import pandas as pd
from subprocess import call
from os import path
import numpy as np

from paramCube import BOUNDS

outputdir='/u/ki/swmclau2/des/EmulatorData/'

for f_c in np.linspace(BOUNDS['f_c'][0], BOUNDS['f_c'][1], num = 5):
    for alpha in np.linspace(BOUNDS['alpha'][0], BOUNDS['alpha'][1], num=5):
        for logM1 in np.linspace(BOUNDS['logM1'][0], BOUNDS['logM'][1], num=5):

            jobname = '_'.join('Emulator_fc_%.2f_a_%.2f_M1_%.2f'%(f_c, alpha, logM1))
            print jobname
            logfile = jobname + '.out'
            command = ['bsub',
                       '-J', jobname,
                       '-o', logfile,
                       '-W', '12:00',
                       'python', path.join(path.dirname(__file__), 'paramCube.py'),
                       outputdir,
                       '--f_c', f_c,
                       '--alpha', alpha,
                       '--logM1', logM1]
            call(command)
