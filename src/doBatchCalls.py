# @Author Sean McLaughlin
# Send a collection of calls to paramCube to ki-ls
from subprocess import call
from os import path
import numpy as np

system = 'sherlock'
if system == 'ki-ls':
    outputdir = '/u/ki/swmclau2/des/EmulatorData/'
elif system == 'sherlock':
    outputdir = '/home/swmclau2/scratch/EmulatorData/'

QUEUE = 'bulletmpi'
N_PER_DIM = 4
TIME = 6 #hours

BOUNDS = {'logMmin': (11.7, 12.5), 'sigma_logM': (0.2, 0.7), 'logM0': (10, 13), 'logM1': (13.1, 14.3),
          'alpha': (0.75, 1.25), 'f_c': (0.1, 0.5)}


def make_kils_command(jobname, params, queue='bulletmpi'):
    logfile = jobname + '.out'
    command = ['bsub',
               '-q', queue,
               '-n', str(16),
               '-J', jobname,
               '-oo', path.join(outputdir, logfile),
               '-W', '%d:00'%TIME,
               'python', path.join(path.dirname(__file__), 'paramCube.py'),
               outputdir]

    param_list = []
    for param, val in params.iteritems():
        param_list.append('--%s' % param)
        param_list.append(str(val))

    command.extend(param_list)

    return command

#not sure if this will work the same way, the spaces may be an issue.
def make_sherlock_command(jobname, params):
    logfile = jobname+'.out'
    errfile = jobname+'.err'

    sbatch_header = ['#!/bin/bash',
               '--job-name=%s'%jobname,
               '-p iric', #KIPAC queu
               '--output=%s'%path.join(outputdir, logfile),
               '--error=%s'%path.join(outputdir, errfile),
               '--time=%d:00'%(TIME*60),
               '--qos=normal',
               '--nodes=%d'%1,
               #'--exclusive',
               '--mem-per-cpu=32000',
               '--ntasks-per-node=%d'%1,
               '--cpus-per-task=%d'%16]

    sbatch_header = '\n#SBATCH '.join(sbatch_header)

    call_str =['python', path.join(path.dirname(__file__), 'paramCube.py'),
               outputdir]

    for param, val in params.iteritems():
        call_str.append('--%s' % param)
        call_str.append(str(val))

    call_str = ' '.join(call_str)

    with open('./tmp.sbatch', 'w') as f:
        f.write(sbatch_header + '\n' + call_str)

    return 'sbatch ./tmp.sbatch'

make_command = make_kils_command if system=='ki-ls' else make_sherlock_command

if __name__ == "__main__":
    for f_c in np.linspace(BOUNDS['f_c'][0], BOUNDS['f_c'][1], num=N_PER_DIM):
        for alpha in np.linspace(BOUNDS['alpha'][0], BOUNDS['alpha'][1], num=N_PER_DIM):
            for logM1 in np.linspace(BOUNDS['logM1'][0], BOUNDS['logM1'][1], num=N_PER_DIM):
                jobname = 'Emulator_fc_%.2f_a_%.2f_M1_%.2f' % (f_c, alpha, logM1)
                print jobname
                params = {'f_c':f_c, 'alpha':alpha, 'logM1':logM1}

                command = make_command(jobname, params)
                call(command, shell = True)
                break
            break
        break
