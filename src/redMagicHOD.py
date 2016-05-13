#!bin/bash
#@Author Sean McLaughlin
#This module contains the components of the Red Magic HOD, a subclass of the Zheng07 HOD built into halotools.

from halotools.empirical_models import Zheng07Cens, Zheng07Sats
import numpy as np

class RedMagicCens(Zheng07Cens):
    """Slight tweak of the Zheng model to add a new parameter, f_c
    """
    #TODO what meaning does Luminosity threshold hold here?
    def __init__(self,**kwargs):

        super(RedMagicCens,self).__init__(**kwargs)

        if 'f_c' not in self.param_dict:
            self.param_dict['f_c'] = 0.19 #add in best fit of new param.

    def mean_occupation(self, **kwargs):
        "See Zheng07 for details"
        return self.param_dict['f_c']*super(RedMagicCens,self).mean_occupation(**kwargs)

class RedMagicSats(Zheng07Sats):
    "Slight tweak of Zheng model to add new parameter, f_c"

    def mean_occupation(self, **kwargs):
        "See Zheng07 for details"
        f_c = 1
        if 'f_c' in self.param_dict:
            f_c = self.param_dict['f_c']

        return super(RedMagicSats,self).mean_occupation(**kwargs)/f_c

class StepFuncCens(Zheng07Cens):
    "Testing HOD that is a step function for centrals"
    def __init__(self, mMin = 7e12, h = 0.7, **kwargs):
        upper_occupation_bound = 1.0
        super(StepFuncCens, self).__init__(**kwargs)
        print 'logMmin: %e'%np.log10(mMin/h)
        self.param_dict['logMmin'] = np.log10(mMin/h)#200 Chinchilla particles

    def mean_occupation(self, **kwargs):
        "See Zheng07 for details"
        if 'table' in kwargs.keys():
            mass = kwargs['table'][self.prim_haloprop_key]
        elif 'prim_haloprop' in kwargs.keys():
            mass = kwargs['prim_haloprop']
        else:
            msg = ("\nYou must pass either a ``table`` or ``prim_haloprop`` argument \n"
                   "to the ``mean_occupation`` function of the ``StepFuncCens`` class.\n")
            raise RuntimeError(msg)
        mass = np.array(mass)
        if np.shape(mass) == ():
            mass = np.array([mass])

        Mmin = 10**self.param_dict['logMmin']

        return np.array(mass > Mmin, dtype = int)

class StepFuncSats(Zheng07Sats):
    "Testing HOD that is 0 for satellites"

    def mean_occupation(self, **kwargs):
        "See Zheng07 for details"
        if 'table' in kwargs.keys():
            mass = kwargs['table'][self.prim_haloprop_key]
        elif 'prim_haloprop' in kwargs.keys():
            mass = kwargs['prim_haloprop']
        else:
            msg = ("\nYou must pass either a ``table`` or ``prim_haloprop`` argument \n"
                   "to the ``mean_occupation`` function of the ``StepFuncCens`` class.\n")
            raise RuntimeError(msg)
        mass = np.array(mass)
        if np.shape(mass) == ():
            mass = np.array([mass])

        return np.zeros_like(mass)
