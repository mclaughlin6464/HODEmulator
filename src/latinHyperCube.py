# /bin/bash
# @Author Sean McLaughlin
# Construct a Latin hyper cube for training the emulator.

from subprocess import call
from os import path
import argparse
from time import time
import numpy as np

# These imports I'd like moved somewhere else, or hosted here
# For the time being, I'll prefer this to a duplication.
from doBatchCalls import BOUNDS
from emulator import PARAMS

# TODO put all this into a config!
# not happy about copying this code
# need to put this into configs so this isn't an issue.
QUEUE = 'bulletmpi'
TIME = 24  # hours
system = 'ki-ls'
if system == 'ki-ls':
    outputdir = '/u/ki/swmclau2/des/EmulatorLHC/'
    # outputdir = '/u/ki/swmclau2/des/TestData/'

elif system == 'sherlock':
    outputdir = '/home/swmclau2/scratch/EmulatorLHC/'


def make_cube(M=500):
    '''Return an array detailing the values of each parameter.'''

    # NOTE I'm not modifying PARAMS at all. A more general case would
    # allow the user to pass in arbitrary dimensions. Probably not necessary.
    np.random.seed(int(time()))

    # ranges= []
    # for p in PARAMS:
    #   ranges.append(np.linspace(BOUNDS[p][0], BOUNDS[p][1], num=M)) )
    #   np.random.shuffle(ranges[-1])
    # return np.stack(ranges).T

    # Unreadable, but i'm a code golf champion.
    return np.stack([np.random.shuffle(np.linspace(BOUNDS[p][0], BOUNDS[p][1], num=M)) for p in PARAMS]).T


def make_kils_command(jobname, params, queue=QUEUE):
    logfile = jobname + '.out'
    command = ['bsub',
               '-q', queue,
               '-n', str(16),
               '-J', jobname,
               '-oo', path.join(outputdir, logfile),
               '-W', '%d:00' % TIME,
               'python', path.join(path.dirname(__file__), 'paramCube.py'),
               outputdir]

    param_list = []
    # param_list.append('--test')
    param_list.append('--id')
    param_list.append(jobname[-3:])
    for param, val in params.iteritems():
        param_list.append('--%s' % param)
        param_list.append(str(val))


    command.extend(param_list)

    return command


def make_sherlock_command(jobname, params):
    logfile = jobname + '.out'
    errfile = jobname + '.err'

    sbatch_header = ['#!/bin/bash',
                     '--job-name=%s' % jobname,
                     '-p iric',  # KIPAC queu
                     '--output=%s' % path.join(outputdir, logfile),
                     '--error=%s' % path.join(outputdir, errfile),
                     '--time=%d:00' % (TIME * 60),
                     '--qos=normal',
                     '--nodes=%d' % 1,
                     # '--exclusive',
                     '--mem-per-cpu=32000',
                     '--ntasks-per-node=%d' % 1,
                     '--cpus-per-task=%d' % 16]

    sbatch_header = '\n#SBATCH '.join(sbatch_header)

    call_str = ['python', path.join(path.dirname(__file__), 'paramCube.py'),
                outputdir]

    call_str.append('--id')
    call_str.append(jobname[-3:])

    for param, val in params.iteritems():
        call_str.append('--%s' % param)
        call_str.append(str(val))

    call_str = ' '.join(call_str)

    with open('./tmp.sbatch', 'w') as f:
        f.write(sbatch_header + '\n' + call_str)

    return 'sbatch ./tmp.sbatch'


make_command = make_kils_command if system == 'ki-ls' else make_sherlock_command


def send_calls(cube):
    for idx, point in enumerate(cube):
        params = dict(zip(PARAMS, point))
        jobname = 'Emulator_lhc_%3d' % idx

        command = make_command(jobname, params)
        call(command, shell=system == 'sherlock')
        break
