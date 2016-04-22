#!bin/bash
#@Author: Sean McLaughlin
#This file contains objects that are little more than holders for data I need to use over and over
#Introducting a little heirarchy to reduce copy-pasting seems like a good idea to me, if a little overkill
#important object at the end is the cat_dict, which links simnames to the objects here.

__all__ = ['Emu', 'Fox', 'MDHR', 'cat_dict']

class Cat(object):

    def __init__(self, simname = 'Cat',
                 loc = '', columns_to_keep = {},
                 halo_finder = 'rockstar', version_name = 'most_recent',
                 Lbox = 1.0, pmass = 1.0, scale_factors = [], filenames = []):

        self.simname = simname
        self.loc = loc
        self.columns_to_keep = columns_to_keep
        self.columns_to_convert = set(["halo_rvir", "halo_rs"])
        self.columns_to_convert = list(self.columns_to_convert & set(self.columns_to_keep.keys()) )
        self.halo_finder = halo_finder
        self.version_name = version_name
        self.Lbox = Lbox
        self.pmass = pmass
        self.scale_factors = scale_factors
        self.redshifts = [1.0/a-1 for a in self.scale_factors] #TODO let user pass in redshift and get a

        self.filenames = filenames
        for i, fname in enumerate(self.filenames):
            self.filenames[i] = self.loc+fname

        self.cache_locs =['/u/ki/swmclau2/des/halocats/hlist_%.2f.list.%s.hdf5'%(a, self.simname)
                                    for a in self.scale_factors]

    def __len__(self):
        return len(self.filenames)

class Hlist(Cat):

    def __init__(self,**kwargs):

        if 'columns_to_keep' not in kwargs or kwargs['columns_to_keep'] is None:
            kwargs['columns_to_keep'] =  {'halo_id': (1, 'i8'), 'halo_upid': (6, 'i8'),
                        'halo_x': (17, 'f4'), 'halo_y': (18, 'f4'),'halo_z': (19, 'f4'),
                        'halo_vx': (20, 'f4'), 'halo_vy': (21, 'f4'), 'halo_vz': (22, 'f4'),
                        'halo_mvir': (10, 'f4'), 'halo_rvir': (11, 'f4'), 'halo_rs': (12, 'f4')}

        if 'simname' not in kwargs or kwargs['simname'] is None:
            kwargs['simname'] = 'Hlist'

        super(self, Hlist).__init__(self, **kwargs)

class OutList(Cat):

    def __init__(self, **kwargs):
        #Default dict here?
        if 'columns_to_keep' not in kwargs or kwargs['columns_to_keep'] is None:
            kwargs['columns_to_keep'] = {'halo_id': (0, 'i8'), 'halo_upid': (36, 'i8'),
                                         'halo_x': (8, 'f4'), 'halo_y': (9, 'f4'),'halo_z': (10, 'f4'),
                                         'halo_vx': (11, 'f4'), 'halo_vy': (12, 'f4'), 'halo_vz': (13, 'f4'),
                                         'halo_mvir': (2, 'f4'),'halo_rvir': (5, 'f4'), 'halo_rs': (6, 'f4')}

        if 'simname' not in kwargs or kwargs['simname'] is None:
            kwargs['simname'] = 'Outlist'

        super(self, OutList).__init__(self, **kwargs)

class Emu(OutList):
    #TODO define as Box000
    #Actually could subclass boxes. Or with Chichilla, handle that as version info
    def __init__(self, **kwargs):

        defaults = {'simname': 'emu', 'loc':'/u/ki/swmclau2/des/emu/Box000/',
                    'Lbox':1050.0,'pmass':3.9876e10,
                    'filenames':['out_%d.list' % i for i in xrange(10)],
                    'scale_factors':[0.25, 0.333, 0.5, 0.540541, 0.588235, 0.645161, 0.714286, 0.8, 0.909091, 1.0] }
        #TODO move this step to the superclass
        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        super(self, Emu).__init__(self, **kwargs)

class Fox(Hlist):

    def __init__(self, **kwargs):
        defaults = {'simname':'fox', 'loc': '/nfs/slac/g/ki/ki23/des/BCCSims/Fox/Lb400/halos/rockstar/output/hlists/',
                    'Lbox': 400.0, 'pmass': 6.58298e8,
                    'filenames': ['hlist_%d' % n for n in [46, 57, 73, 76, 79, 82, 86, 90, 95, 99]],
                    'scale_factors': [0.25, 0.333, 0.5, 0.540541, 0.588235, 0.645161, 0.714286, 0.8, 0.909091, 1.0] }

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        super(self, Fox).__init__(self, **kwargs)

class MDHR(Hlist):

    def __init__(self, **kwargs):
        defaults = {'simname': 'multidark_highres', 'loc': '/nfs/slac/g/ki/ki20/cosmo/behroozi/MultiDark/hlists/',
                    'Lbox': 1e3, 'pmass': 8.721e9,
                    'scale_factors': [0.25690, 0.34800, 0.49990, 0.53030, 0.65180, 0.71250, 0.80370, 0.91000, 1.00110]}

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        if 'filenames' not in kwargs or kwargs['filenames'] is None:
            kwargs['filenames'] = ['hlist_%.5f.list' % a for a in kwargs['scale_factors']]

        super(self, MDHR).__init__(self, **kwargs)

class Chinchilla(Hlist):

    #Lbox and npart are required!
    def __init__(self,Lbox, npart, **kwargs):
        from glob import glob
        #NOTE not sure if loc should be in default, or pmass for that matter
        defaults = {'simname':'chinchilla', 'loc':'/nfs/slac/g/ki/ki21/cosmo/yymao/sham_test/resolution-test/'}

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        #TODO make a set of valid version_names

        kwargs['version_name'] = 'Lb%d-%d'%(int(Lbox), npart )
        kwargs['loc'] += 'c%d-%d/rockstar/hlists/'%(int(Lbox), npart )
        #TODO make it possible to do cuts on scale factor like "use only cats for z < 1".
        fnames =  glob(kwargs['loc']+ 'hlist_*.list') #snag all the hlists
        #TODO write code that makes it so when I pass in a particular sf the filenames are also cut down.
        #TODO also write it so that if filenames is not as long as sf throw an error
        kwargs['filenames'] = [fname[len(kwargs['loc']):] for fname in fnames] #just want the names in the dir
        kwargs['scale_factors'] = [float(fname[len(kwargs['loc'])+6:-5] for fname in fnames)] #pull out scale factors
        kwargs['pmass'] =  #TODO figure out general relationship of params and pmass


cat_dict = {'emu': Emu, 'fox': Fox, 'multidark_highres': MDHR}
