#@Author Sean McLaughlin
#Send a collection of calls to paramCube to ki-ls
import numpy as np
import pandas as pd
from subprocess import call
from os import path
import numpy as np

from paramCube import BOUNDS

outputdir='/u/ki/swmclau2/des/EmulatorData/'

queue = 'medium'

for f_c in np.linspace(BOUNDS['f_c'][0], BOUNDS['f_c'][1], num = 5):
    for alpha in np.linspace(BOUNDS['alpha'][0], BOUNDS['alpha'][1], num=5):
        for logM1 in np.linspace(BOUNDS['logM1'][0], BOUNDS['logM1'][1], num=5):

            jobname = 'Emulator_fc_%.2f_a_%.2f_M1_%.2f'%(f_c, alpha, logM1)
            print jobname

            logfile = jobname + '.out'
            command = ['bsub',
                       '-q', queue,
                       '-n', str(16),
                       '-J', jobname,
                       '-oo', path.join(outputdir, logfile),
                       '-W', '12:00',
                       'python', path.join(path.dirname(__file__), 'paramCube.py'),
                       outputdir,
                       '--f_c', str(f_c),
                       '--alpha', str(alpha),
                       '--logM1', str(logM1)]
            call(command)
            break
        break
    break
