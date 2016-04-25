#!bin/bash
#@Author: Sean McLaughlin
#This file contains objects that are little more than holders for data I need to use over and over
#Introducting a little heirarchy to reduce copy-pasting seems like a good idea to me, if a little overkill
#important object at the end is the cat_dict, which links simnames to the objects here.

__all__ = ['Multidark','Emu', 'Fox', 'MDHR','Chinchilla', 'cat_dict']

class Cat(object):

    def __init__(self, simname = 'Cat',
                 loc = '', columns_to_keep = {},
                 halo_finder = 'rockstar', version_name = 'most_recent',
                 Lbox = 1.0, pmass = 1.0, scale_factors = [], filenames = [], **kwargs):

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

        assert len(self.filenames) == len(self.redshifts)
        assert len(self.scale_factors) == len(self.redshifts)

        self.cache_locs =['/u/ki/swmclau2/des/halocats/hlist_%.2f.list.%s.hdf5'%(a, self.simname)
                                    for a in self.scale_factors]

    def __len__(self):
        return len(self.filenames)

    def __str__(self):
        output = []
        output.append(self.simname)
        output.append('-'*25)
        output.append('Halo finder:\t %s'%self.halo_finder)
        output.append('Version name:\t%s'%self.version_name)
        output.append('Redshifts:\t%s'%str(self.redshifts))
        output.append('-'*25)
        output.append('Location:\t%s'%self.loc)
        output.append('Lbox:\t%.1f'%self.Lbox)
        output.append('Particle Mass:\t%f'%self.pmass)
        output.append('Columns to Keep:\n%s'%str(self.columns_to_keep))

        return '\n'.join(output)

class Hlist(Cat):

    def __init__(self,**kwargs):

        if 'columns_to_keep' not in kwargs or kwargs['columns_to_keep'] is None:
            kwargs['columns_to_keep'] =  {'halo_id': (1, 'i8'), 'halo_upid': (6, 'i8'),
                        'halo_x': (17, 'f4'), 'halo_y': (18, 'f4'),'halo_z': (19, 'f4'),
                        'halo_vx': (20, 'f4'), 'halo_vy': (21, 'f4'), 'halo_vz': (22, 'f4'),
                        'halo_mvir': (10, 'f4'), 'halo_rvir': (11, 'f4'), 'halo_rs': (12, 'f4')}

        if 'simname' not in kwargs or kwargs['simname'] is None:
            kwargs['simname'] = 'Hlist'

        super(Hlist,self).__init__(**kwargs)

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

        super(OutList,self).__init__(**kwargs)

class Multidark(Cat):
    'Builtin, so needs very little. Will never be cached!'
    def __init__(self, **kwargs):
        if 'simname' not in kwargs or kwargs['simname'] is None:
            kwargs['simname'] = 'multidark'

        if 'scale_factors' not in kwargs or kwargs['scale_factors'] is None:
            kwargs['scale_factors'] = [2.0/3, 1.0]

        super(Multidark, self).__init__(**kwargs)

class Emu(OutList):
    #TODO define as Box000
    #Actually could subclass boxes. Or with Chichilla, handle that as version info
    def __init__(self, **kwargs):

        defaults = {'simname': 'emu', 'loc':'/u/ki/swmclau2/des/emu/Box000/',
                    'Lbox':1050.0,'pmass':3.9876e10,
                    'filenames':['out_%d.list' % i for i in xrange(10)],
                    'scale_factors':[0.25, 0.333, 0.5, 0.540541, 0.588235, 0.645161, 0.714286, 0.8, 0.909091, 1.0] }

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        super(Emu,self).__init__(**kwargs)

class Fox(Hlist):

    def __init__(self, **kwargs):
        defaults = {'simname':'fox', 'loc': '/nfs/slac/g/ki/ki23/des/BCCSims/Fox/Lb400/halos/rockstar/output/hlists/',
                    'Lbox': 400.0, 'pmass': 6.58298e8,
                    'filenames': ['hlist_%d' % n for n in [46, 57, 73, 76, 79, 82, 86, 90, 95, 99]],
                    'scale_factors': [0.25, 0.333, 0.5, 0.540541, 0.588235, 0.645161, 0.714286, 0.8, 0.909091, 1.0] }

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        super(Fox, self).__init__(**kwargs)

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

        super(MDHR,self).__init__(**kwargs)

class Chinchilla(Hlist):

    #Lbox and npart are required!
    def __init__(self, **kwargs):

        assert 'Lbox' in kwargs and 'npart' in kwargs

        from glob import glob
        #NOTE not sure if loc should be in default, or pmass for that matter
        defaults = {'simname':'chinchilla', 'loc':'/nfs/slac/g/ki/ki21/cosmo/yymao/sham_test/resolution-test/',
                    'pmass':  1.44390e+08} #mass for 125-1024}
        #TODO make it possible to do cuts on scale factor like "use only cats for z < 1".

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        valid_version_names = set(['Lb125-1024','Lb125-2048','Lb250-1024','Lb250-128',
                                   'Lb250-196','Lb250-2048', 'Lb250-2560', 'Lb250-320',
                                   'Lb250-512', 'Lb250-768', 'Lb250-85', 'Lb400-1024',
                                   'Lb400-136', 'Lb400-2048', 'Lb400-210', 'Lb400-315',
                                   'Lb400-512', 'Lb400-768'])

        if 'version_name' not in kwargs:
            kwargs['version_name'] = 'Lb%d-%d'%(int(kwargs['Lbox']), kwargs['npart'] )
        else:# check version name is legit
            split_vname = kwargs['version_name'].split('-')
            assert int(split_vname[0][2:]) == kwargs['Lbox']
            assert int(split_vname[1]) == kwargs['npart']

        assert kwargs['version_name'] in valid_version_names
            #raise ValueError('%s is not a valid version of %s'%(kwargs['version_name'], kwargs['simname']))

        kwargs['loc'] += 'c%d-%d/rockstar/hlists/'%(int(kwargs['Lbox']), kwargs['npart'] )
        fnames =  glob(kwargs['loc']+ 'hlist_*.list') #snag all the hlists
        fnames = [fname[len(kwargs['loc']):] for fname in fnames] #just want the names in the dir
        tmp_scale_factors = [float(fname[6:-5]) for fname in fnames] #pull out scale factors

        #if the user passed in stuff, have to check a bunch of things
        if 'filenames' not in kwargs:
            kwargs['filenames'] = fnames
        elif 'scale_factors' in kwargs:#don't know why this case would ever be true
            assert len(kwargs['filenames'] ) == len(kwargs['scale_factors'])
            for kw_fname in kwargs['filenames']:
                assert kw_fname in fnames
            #do nothing, we're good.
        else:
            kwargs['scale_factors'] = []
            for kw_fname in kwargs['filenames']:
                assert kw_fname in fnames
                kwargs['scale_factors'].append(tmp_scale_factors[fnames.index(kw_fname)]) #get teh matching scale factor

        if 'scale_factors' not in kwargs:
            kwargs['scale_factors'] = tmp_scale_factors
        else: #Don't have to do the both case, covered above
            kwargs['filenames'] = []
            for a in kwargs['scale_factors']:
                assert a in tmp_scale_factors
                kwargs['filenames'].append(fnames[tmp_scale_factors.index(a)])  # get teh matching scale factor

        kwargs['pmass']*=((kwargs['Lbox']/125.0)**3)*((1024.0/kwargs['npart'])**3) #correct factor for right pmass

        super(Chinchilla, self).__init__(**kwargs)

        self.cache_locs = ['/u/ki/swmclau2/des/halocats/hlist_%.2f.list.%s_%s.hdf5' % (a, self.simname, self.version_name)
                           for a in self.scale_factors] #make sure we don't have redunancies.


cat_dict = {'multidark':Multidark,'emu': Emu, 'fox': Fox, 'multidark_highres': MDHR, 'chinchilla': Chinchilla}
