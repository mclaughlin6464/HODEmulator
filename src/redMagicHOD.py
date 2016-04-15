#!bin/bash
#@Author Sean McLaughlin
#This module contains the components of the Red Magic HOD, a subclass of the Zheng07 HOD built into halotools.

from halotools.empirical_models import Zheng07Cens, Zheng07Sats

class RedMagicCens(Zheng07Cens):
    """Slight tweak of the Zheng model to add a new parameter, f_c
    """
    #TODO what meaning does Luminosity threshold hold here?
    def __init__(self,**kwargs):

        upper_occupation_bound = 1.0
        super(RedMagicCens,self).__init__(**kwargs)

        if 'f_c' not in self.param_dict:
            self.param_dict['f_c'] = 0.19 #add in best fit of new param.

    def mean_occupation(self, **kwargs):
        "See Zheng07 for details"
        return self.param_dict['f_c']*super(RedMagicCens,self).mean_occupation(**kwargs)

class RedMagicSats(Zheng07Sats):
    "Slight tweak of Zheng model to add new parameter, f_c"

    def __init__(self, **kwargs):

        upper_occupation_bound = float("inf")
        super(RedMagicSats,self).__init__(**kwargs)

        #TODO not sure if this will work for the whole model; will need to test
        #may need to 'modulate'
        #if 'f_c' not in self.param_dict:
        #    self.param_dict['f_c'] = 0.19 #add in best fit of new param.

    def mean_occupation(self, **kwargs):
        "See Zheng07 for details"
        f_c = 1
        if 'f_c' in self.param_dict:
            f_c = self.param_dict['f_c']

        return super(RedMagicSats,self).mean_occupation(**kwargs)/f_c