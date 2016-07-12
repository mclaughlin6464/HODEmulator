#!bin/bash
# @Author: Sean McLaughlin
# This file contains objects that are little more than holders for data I need to use over and over
# Introducting a little heirarchy to reduce copy-pasting seems like a good idea to me, if a little overkill
# important object at the end is the cat_dict, which links simnames to the objects here.

from astropy import cosmology
from socket import gethostname

__all__ = ['Bolshoi', 'Multidark', 'Emu', 'Fox', 'MDHR', 'Chinchilla', 'Aardvark', 'Guppy', 'Chinchilla1050', 'Emu200',
           'Fox1050', 'cat_dict']

hostname = gethostname()
KILS = hostname[:-2] == 'ki-ls'
KILS = True  # TODO fixme

# TODO each cat should carry a default output script, to which specific information is added.

# set filepaths depending on which cluster we're on.
if KILS:
    DEFAULT_LOCS = {'emu': '/u/ki/swmclau2/des/emu/Box000/',
                    'emu200': '/u/ki/swmclau2/des/emu200/',
                    'fox': '/nfs/slac/g/ki/ki23/des/BCCSims/Fox/Lb400/halos/rockstar/output/hlists/',
                    'fox1050': '/u/ki/swmclau2/des/fox/',
                    'multidark_highres': '/nfs/slac/g/ki/ki20/cosmo/behroozi/MultiDark/hlists/',
                    'chinchilla': '/nfs/slac/g/ki/ki21/cosmo/yymao/sham_test/resolution-test/',
                    'aardvark': '/nfs/slac/g/ki/ki18/des/mbusha/simulations/Aardvark-2PS/Lb400/rockstar/hlists/',
                    'guppy': '/u/ki/swmclau2/des/guppy/',
                    'chinchilla1050': '/u/ki/swmclau2/des/rockstar_outputs/'}

    CACHE_LOCS = {'cat': '/u/ki/swmclau2/des/halocats/hlist_%.2f.list.%s.hdf5',
                  'chinchilla': '/u/ki/swmclau2/des/halocats/hlist_%.2f.list.%s_%s.hdf5'}

else:  # On Sherlock

    DEFAULT_LOCS = {'emu': '/scratch/users/swmclau2/hlists/emu/Box000/',
                    'fox': '/scratch/users/swmclau2/hlists/Fox/',
                    'multidark_highres': '/scratch/users/swmclau2/hlists/MDHR/',
                    'chinchilla': '/scratch/users/swmclau2/hlists/Chinchilla/',
                    'emu200': '/scratch/PI/kipac/yymao/highres_emu/Box000/hlists/'}

    CACHE_LOCS = {'cat': '/scratch/users/swmclau2/halocats/hlist_%.2f.list.%s.hdf5',
                  'chinchilla': '/scratch/users/swmclau2/halocats/hlist_%.2f.list.%s_%s.hdf5'}

HLIST_COLS = {'halo_id': (1, 'i8'), 'halo_upid': (6, 'i8'),
              'halo_x': (17, 'f4'), 'halo_y': (18, 'f4'), 'halo_z': (19, 'f4'),
              'halo_vx': (20, 'f4'), 'halo_vy': (21, 'f4'), 'halo_vz': (22, 'f4'),
              'halo_mvir': (10, 'f4'), 'halo_rvir': (11, 'f4'), 'halo_rs': (12, 'f4')}

OUTLIST_COLS = {'halo_id': (0, 'i8'), 'halo_upid': (36, 'i8'),
                'halo_x': (8, 'f4'), 'halo_y': (9, 'f4'), 'halo_z': (10, 'f4'),
                'halo_vx': (11, 'f4'), 'halo_vy': (12, 'f4'), 'halo_vz': (13, 'f4'),
                'halo_mvir': (2, 'f4'), 'halo_rvir': (5, 'f4'), 'halo_rs': (6, 'f4')}


class Cat(object):
    def __init__(self, simname='Cat',
                 loc='', columns_to_keep={},
                 halo_finder='rockstar', version_name='most_recent',
                 Lbox=1.0, pmass=1.0, scale_factors=[], cosmo=cosmology.WMAP5,
                 filenames=[], **kwargs):

        self.simname = simname
        self.loc = loc
        self.columns_to_keep = columns_to_keep
        self.columns_to_convert = set(["halo_rvir", "halo_rs"])
        self.columns_to_convert = list(self.columns_to_convert & set(self.columns_to_keep.keys()))
        self.halo_finder = halo_finder
        self.version_name = version_name
        self.Lbox = Lbox
        self.pmass = pmass

        self.scale_factors = scale_factors
        self.redshifts = [1.0 / a - 1 for a in self.scale_factors]  # TODO let user pass in redshift and get a

        self.cosmology = cosmo  # default cosmology
        self.h = self.cosmology.H(0).value / 100.0

        self.filenames = filenames
        for i, fname in enumerate(self.filenames):
            self.filenames[i] = self.loc + fname

        assert (len(self.filenames) == len(self.redshifts)) or len(self.filenames) == 0  # built ins have no filenames
        assert len(self.scale_factors) == len(self.redshifts)

        self.cache_locs = [CACHE_LOCS['cat'] % (a, self.simname)
                           for a in self.scale_factors]

    def __len__(self):
        return len(self.filenames)

    def __str__(self):
        output = []
        output.append(self.simname)
        output.append('-' * 25)
        output.append('Halo finder:\t %s' % self.halo_finder)
        output.append('Version name:\t%s' % self.version_name)
        output.append('Cosmology:\n%s' % self.cosmology)
        # output.append('Redshifts:\t%s'%str(self.redshifts))
        output.append('-' * 25)
        # output.append('Location:\t%s'%self.loc)
        output.append('Lbox:\t%.1f' % self.Lbox)
        output.append('Particle Mass:\t%f' % self.pmass)
        output.append('Columns to Keep:\n%s' % str(self.columns_to_keep))

        return '\n'.join(output)

    def update_lists(self, user_kwargs, tmp_fnames, tmp_scale_factors):
        '''If the user passes in a scale factor or filename, we have to do some cropping'''
        if 'filenames' not in user_kwargs:
            user_kwargs['filenames'] = tmp_fnames
        elif 'scale_factors' in user_kwargs:  # don't know why this case would ever be true
            assert len(user_kwargs['filenames']) == len(user_kwargs['scale_factors'])
            for kw_fname in user_kwargs['filenames']:
                assert kw_fname in tmp_fnames
                # do nothing, we're good.
        else:
            user_kwargs['scale_factors'] = []
            for kw_fname in user_kwargs['filenames']:
                assert kw_fname in tmp_fnames
                user_kwargs['scale_factors'].append(
                    tmp_scale_factors[tmp_fnames.index(kw_fname)])  # get teh matching scale factor

        if 'scale_factors' not in user_kwargs:
            user_kwargs['scale_factors'] = tmp_scale_factors
        else:  # both case covered above.
            user_kwargs['filenames'] = []
            for a in user_kwargs['scale_factors']:
                assert a in tmp_scale_factors
                user_kwargs['filenames'].append(tmp_fnames[tmp_scale_factors.index(a)])  # get teh matching scale factor


class Hlist(Cat):
    def __init__(self, **kwargs):

        if 'columns_to_keep' not in kwargs or kwargs['columns_to_keep'] is None:
            kwargs['columns_to_keep'] = {'halo_id': (1, 'i8'), 'halo_upid': (6, 'i8'),
                                         'halo_x': (17, 'f4'), 'halo_y': (18, 'f4'), 'halo_z': (19, 'f4'),
                                         'halo_vx': (20, 'f4'), 'halo_vy': (21, 'f4'), 'halo_vz': (22, 'f4'),
                                         'halo_mvir': (10, 'f4'), 'halo_rvir': (11, 'f4'), 'halo_rs': (12, 'f4')}

        if 'simname' not in kwargs or kwargs['simname'] is None:
            kwargs['simname'] = 'Hlist'

        super(Hlist, self).__init__(**kwargs)


class OutList(Cat):
    def __init__(self, **kwargs):
        # Default dict here?
        if 'columns_to_keep' not in kwargs or kwargs['columns_to_keep'] is None:
            kwargs['columns_to_keep'] = {'halo_id': (0, 'i8'), 'halo_upid': (36, 'i8'),
                                         'halo_x': (8, 'f4'), 'halo_y': (9, 'f4'), 'halo_z': (10, 'f4'),
                                         'halo_vx': (11, 'f4'), 'halo_vy': (12, 'f4'), 'halo_vz': (13, 'f4'),
                                         'halo_mvir': (2, 'f4'), 'halo_rvir': (5, 'f4'), 'halo_rs': (6, 'f4')}

        if 'simname' not in kwargs or kwargs['simname'] is None:
            kwargs['simname'] = 'Outlist'

        super(OutList, self).__init__(**kwargs)


class Multidark(Cat):
    'Builtin, so needs very little. Will never be cached!'

    def __init__(self, **kwargs):
        if 'simname' not in kwargs or kwargs['simname'] is None:
            kwargs['simname'] = 'multidark'

        if 'version_name' not in kwargs or kwargs['version_name'] is None:
            kwargs['version_name'] = 'halotools_alpha_version2'

        if 'scale_factors' not in kwargs or kwargs['scale_factors'] is None:
            kwargs['scale_factors'] = [2.0 / 3, 1.0]

        super(Multidark, self).__init__(**kwargs)


class Bolshoi(Cat):
    'Builtin, so needs very little. Will never be cached!'

    def __init__(self, **kwargs):
        if 'simname' not in kwargs or kwargs['simname'] is None:
            kwargs['simname'] = 'bolshoi'

        if 'version_name' not in kwargs or kwargs['version_name'] is None:
            kwargs['version_name'] = 'halotools_alpha_version2'

        if 'scale_factors' not in kwargs or kwargs['scale_factors'] is None:
            kwargs['scale_factors'] = [1.0]

        super(Bolshoi, self).__init__(**kwargs)


class Emu(OutList):
    # TODO define as Box000
    # Actually could subclass boxes. Or with Chichilla, handle that as version info
    def __init__(self, **kwargs):

        defaults = {'simname': 'emu', 'loc': DEFAULT_LOCS['emu'],
                    'Lbox': 1050.0, 'pmass': 3.9876e10,
                    'filenames': ['out_%d.list' % i for i in xrange(10)],
                    'cosmo': cosmology.core.wCDM(H0=63.36569, Om0=0.340573, Ode0=0.659427, w0=-0.816597),
                    'scale_factors': [0.25, 0.333, 0.5, 0.540541, 0.588235, 0.645161, 0.714286, 0.8, 0.909091, 1.0]}

        tmp_scale_factors = defaults['scale_factors']
        tmp_fnames = defaults['filenames']

        self.update_lists(kwargs, tmp_fnames, tmp_scale_factors)

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        super(Emu, self).__init__(**kwargs)


class Chinchilla1050(OutList):
    # pretty different from chinchilla, so for the time being writing as a separate object.
    def __init__(self, **kwargs):

        defaults = {'simname': 'chinchilla1050', 'loc': DEFAULT_LOCS['chinchilla1050'],
                    'Lbox': 1050.0, 'pmass': 3.34881e+10,
                    'filenames': ['out_%d.list' % i for i in xrange(10)],
                    'cosmo': cosmology.core.LambdaCDM(H0=100 * 0.7, Om0=0.286, Ode0=0.714),
                    'scale_factors': [0.1429, 0.1667, 0.2, 0.25, 0.3333, 0.4, 0.5, 0.6667, 0.8, 1.0]}

        tmp_scale_factors = defaults['scale_factors']
        tmp_fnames = defaults['filenames']

        self.update_lists(kwargs, tmp_fnames, tmp_scale_factors)

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        super(Chinchilla1050, self).__init__(**kwargs)


class Guppy(OutList):
    def __init__(self, **kwargs):

        defaults = {'simname': 'guppy', 'loc': DEFAULT_LOCS['guppy'],
                    'Lbox': 1050.0, 'pmass': 3.45420e+10,
                    'filenames': ['out_%d.list' % i for i in xrange(10)],
                    'cosmo': cosmology.core.LambdaCDM(H0=100 * 0.6881, Om0=0.295, Ode0=0.705),
                    'scale_factors': [0.1429, 0.1667, 0.2, 0.25, 0.3333, 0.4, 0.5, 0.6667, 0.8, 1.0]}

        tmp_scale_factors = defaults['scale_factors']
        tmp_fnames = defaults['filenames']

        self.update_lists(kwargs, tmp_fnames, tmp_scale_factors)

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        super(Guppy, self).__init__(**kwargs)


class Fox(Hlist):
    def __init__(self, **kwargs):
        defaults = {'simname': 'fox', 'loc': DEFAULT_LOCS['fox'],
                    'Lbox': 400.0, 'pmass': 6.58298e8,
                    'filenames': ['hlist_%d' % n for n in [46, 57, 73, 76, 79, 82, 86, 90, 95, 99]],
                    'cosmo': cosmology.core.LambdaCDM(H0=100 * 0.6704346, Om0=0.318340, Ode0=0.681660),
                    'scale_factors': [0.25, 0.333, 0.5, 0.540541, 0.588235, 0.645161, 0.714286, 0.8, 0.909091, 1.0]}

        tmp_scale_factors = defaults['scale_factors']
        tmp_fnames = defaults['filenames']

        self.update_lists(kwargs, tmp_fnames, tmp_scale_factors)

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        super(Fox, self).__init__(**kwargs)


class MDHR(Hlist):
    def __init__(self, **kwargs):
        defaults = {'simname': 'multidark_highres', 'loc': DEFAULT_LOCS['multidark_highres'],
                    'Lbox': 1e3, 'pmass': 8.721e9,
                    'scale_factors': [0.25690, 0.34800, 0.49990, 0.53030, 0.65180, 0.71250, 0.80370, 0.91000, 1.00110]}

        # if 'filenames' not in kwargs or kwargs['filenames'] is None:
        #    kwargs['filenames'] = ['hlist_%.5f.list' % a for a in kwargs['scale_factors']]
        tmp_scale_factors = defaults['scale_factors']
        tmp_fnames = ['hlist_%.5f.list' % a for a in tmp_scale_factors]

        self.update_lists(kwargs, tmp_fnames, tmp_scale_factors)

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        super(MDHR, self).__init__(**kwargs)


# NOTE not sure how to subclass this. It's like a less complex Chinchilla.
# could roll both of them under a BCC class or just leave them separate.
# Will start writing for now and see what happens as I proceed.

class Aardvark(Hlist):
    # Lbox technically required, but I don't even have access to anything besides 400. Ignore for now.
    def __init__(self, **kwargs):

        defaults = {'simname': 'aardvark', 'loc': DEFAULT_LOCS['aardvark'],
                    'cosmo': cosmology.core.LambdaCDM(H0=100 * 0.73, Om0=0.23, Ode0=0.77, Ob0=0.047),
                    'pmass': 4.75619e+08, 'Lbox': 400.0}

        from glob import glob

        tmp_fnames = glob(kwargs['loc'] + 'hlist_*.list')  # snag all the hlists
        tmp_fnames = [fname[len(kwargs['loc']):] for fname in tmp_fnames]  # just want the names in the dir
        tmp_scale_factors = [float(fname[6:-5]) for fname in tmp_fnames]  # pull out scale factors

        # Looked into a way to put this in the global init.
        # However, the amount of copy-pasting that would save would be minimal, it turns out.
        self.update_lists(kwargs, tmp_fnames, tmp_scale_factors)

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        super(Aardvark, self).__init__(**kwargs)


class Emu200(Hlist):
    # TODO redesign the heiharchy here so this will be within Emu

    def __init__(self, **kwargs):

        defaults = {'simname': 'emu200', 'loc': DEFAULT_LOCS['emu200'],
                    'cosmo': cosmology.core.wCDM(H0=100 * 0.6616172, Om0=0.309483642394, Ode0=0.690516357606,
                                                 w0=-0.8588491),
                    'pmass': 6.3994e8, 'Lbox': 200.0}

        from glob import glob
        tmp_fnames = glob(kwargs['loc'] + 'hlist_*.list')  # snag all the hlists
        tmp_fnames = [fname[len(kwargs['loc']):] for fname in tmp_fnames]  # just want the names in the dir
        tmp_scale_factors = [float(fname[6:-5]) for fname in tmp_fnames]  # pull out scale factors

        # Looked into a way to put this in the global init.
        # However, the amount of copy-pasting that would save would be minimal, it turns out.
        self.update_lists(kwargs, tmp_fnames, tmp_scale_factors)

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        super(Emu200, self).__init__(**kwargs)


class Fox1050(OutList):
    def __init__(self, **kwargs):
        defaults = {'simname': 'fox1050', 'loc': DEFAULT_LOCS['fox1050'],
                    'Lbox': 1050.0, 'pmass': 3.72749e+10,
                    'filenames': ['out_%d.list' % i for i in xrange(10)],
                    'cosmo': cosmology.core.LambdaCDM(H0=100 * 0.6704346, Om0=0.318340, Ode0=0.681660),
                    'scale_factors': [0.2, 0.222222, 0.25, 0.287514, 0.333333, 0.4, 0.5, 0.666667, 0.8, 1.0]}

        tmp_scale_factors = defaults['scale_factors']
        tmp_fnames = defaults['filenames']

        self.update_lists(kwargs, tmp_fnames, tmp_scale_factors)

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        super(Fox1050, self).__init__(**kwargs)


class Chinchilla(Hlist):
    # Lbox and npart are required!
    def __init__(self, **kwargs):

        assert 'Lbox' in kwargs and 'npart' in kwargs

        from glob import glob
        # NOTE not sure if loc should be in default, or pmass for that matter
        defaults = {'simname': 'chinchilla', 'loc': DEFAULT_LOCS['chinchilla'],
                    'cosmo': cosmology.core.LambdaCDM(H0=100 * 0.7, Om0=0.286, Ode0=0.714),
                    'pmass': 1.44390e+08}  # mass for 125-1024}
        # TODO make it possible to do cuts on scale factor like "use only cats for z < 1".

        for key, value in defaults.iteritems():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = value

        valid_version_names = {'Lb125-1024', 'Lb125-2048', 'Lb250-1024', 'Lb250-128',
                               'Lb250-196', 'Lb250-2048', 'Lb250-2560', 'Lb250-320',
                               'Lb250-512', 'Lb250-768', 'Lb250-85', 'Lb400-1024',
                               'Lb400-136', 'Lb400-2048', 'Lb400-210', 'Lb400-315',
                               'Lb400-512', 'Lb400-768'}

        if 'version_name' not in kwargs:
            kwargs['version_name'] = 'Lb%d-%d' % (int(kwargs['Lbox']), kwargs['npart'])
        else:  # check version name is legit
            split_vname = kwargs['version_name'].split('-')
            assert int(split_vname[0][2:]) == kwargs['Lbox']
            assert int(split_vname[1]) == kwargs['npart']

        assert kwargs['version_name'] in valid_version_names
        # raise ValueError('%s is not a valid version of %s'%(kwargs['version_name'], kwargs['simname']))

        kwargs['loc'] += 'c%d-%d/' % (int(kwargs['Lbox']), kwargs['npart'])

        if KILS:  # differences in how the files are stored.
            kwargs['loc'] += '/rockstar/hlists/'

        tmp_fnames = glob(kwargs['loc'] + 'hlist_*.list')  # snag all the hlists
        tmp_fnames = [fname[len(kwargs['loc']):] for fname in tmp_fnames]  # just want the names in the dir
        tmp_scale_factors = [float(fname[6:-5]) for fname in tmp_fnames]  # pull out scale factors

        # if the user passed in stuff, have to check a bunch of things
        self.update_lists(kwargs, tmp_fnames, tmp_scale_factors)

        kwargs['pmass'] *= ((kwargs['Lbox'] / 125.0) ** 3) * (
        (1024.0 / kwargs['npart']) ** 3)  # correct factor for right pmass

        super(Chinchilla, self).__init__(**kwargs)

        self.cache_locs = [CACHE_LOCS['chinchilla'] % (a, self.simname, self.version_name)
                           for a in self.scale_factors]  # make sure we don't have redunancies.


cat_dict = {'bolshoi': Bolshoi, 'multidark': Multidark, 'emu': Emu, 'fox': Fox, 'multidark_highres': MDHR,
            'chinchilla': Chinchilla,
            'aardvark': Aardvark, 'guppy': Guppy, 'chinchilla1050': Chinchilla1050, 'emu200': Emu200,
            'fox1050': Fox1050}
