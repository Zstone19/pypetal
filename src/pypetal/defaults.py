import numpy as np
from astropy.table import Table




def set_general(input_args, fnames):
    
    default_kwargs = {
        'verbose': False,
        'plot': False,
        'time_unit': 'd',
        'lc_unit': 'Arbitrary Units',
        'file_fmt': 'csv',
        'lag_bounds': 'baseline',
        'threads': 1
    }
    
    params = { **default_kwargs, **input_args }

    verbose = params['verbose']
    plot = params['plot']
    time_unit = params['time_unit']
    lc_unit = params['lc_unit']
    file_fmt = params['file_fmt']
    lag_bounds = params['lag_bounds']
    threads = params['threads']
    
    if isinstance(lc_unit, str):
        lc_unit = list( np.full( len(fnames), lc_unit ) )
        
     
     
     
    len_lag_bounds = len(lag_bounds)
    if len_lag_bounds == 2:
        if not ( isinstance(lag_bounds[0], list) or ( lag_bounds[0] == 'baseline' ) ):
            len_lag_bounds = 1
        else:
            len_lag_bounds = 2


    if len_lag_bounds != len(fnames)-1:
        
        lag_bounds_og = lag_bounds
        lag_bounds = []

        for i in range(len(fnames)-1):
            lag_bounds.append(lag_bounds_og)
            
    if (len(fnames) == 2) and ( not isinstance(lag_bounds[0], list) ):
        lag_bounds = [lag_bounds]
            
            
            
    xvals_tot = []
    baselines = []
    for i in range(len(fnames)):
        
        try:
            dat = Table.read(fnames[i], format=file_fmt)
        except:
            dat = Table.read(fnames[i], format='ascii')
            
        colnames = dat.colnames       
        x = np.array( dat[colnames[0]] )

        sort_ind = np.argsort(x)
        x = x[sort_ind]

        xvals_tot.append(x)
        
    for i in range(len(fnames)-1):
        baseline = np.max([ xvals_tot[0].max(), xvals_tot[i+1].max() ]) - np.min([ xvals_tot[0].min(), xvals_tot[i+1].min() ])
        baselines.append(baseline)    
    
    
    for i in range(len(fnames)-1):
        
        if (lag_bounds[i] is None) | (lag_bounds[i] == 'baseline'):
            lag_bounds[i] = [-baselines[i], baselines[i]]

            
    #Return dict so it can be passed to other functions
    output = {
        'verbose': verbose,
        'plot': plot,
        'time_unit': time_unit,
        'lc_unit': lc_unit,
        'file_fmt': file_fmt,
        'lag_bounds': lag_bounds,
        'threads': threads
    }
            
    return output
        
        
    


def set_drw_rej(input_args, fnames):
    
    default_kwargs = {
        'jitter' : True,
        'nsig' : 3,
        'nwalker': 100,
        'nburn': 300,
        'nchain': 1000,
        'clip': np.full( len(fnames), True),
        'reject_data': np.hstack([ [True], np.full( len(fnames)-1, False) ]),
        'use_for_javelin': False
    }
    
    params = { **default_kwargs, **input_args }
    
    jitter = params['jitter']
    nsig = params['nsig']
    nwalker = params['nwalker']
    nburn = params['nburn']
    nchain = params['nchain']
    clip = params['clip']
    reject_data = params['reject_data']
    use_for_javelin = params['use_for_javelin']
    
    if isinstance(clip, bool):
        clip = np.full( len(fnames), clip)
        
    if isinstance(reject_data, bool):
        reject_data = np.full( len(fnames), reject_data)
        
    if (use_for_javelin) & (not reject_data[0]):
        raise Exception('Cannot use continuum for Javelin without fitting it to a DRW first')
    
    
    return jitter, nsig, nwalker, nburn, nchain, clip, \
        reject_data, use_for_javelin


def set_detrend(input_args):

    default_kwargs = {
        'K': 2,
        'nchain': 4,
        'miniter': 5000,
        'maxiter': 10000
    }

    params = { **default_kwargs, **input_args }

    K = params['K']
    nchain = params['nchain']
    miniter = params['miniter']
    maxiter = params['maxiter']

    return K, nchain, miniter, maxiter


def set_pyccf(input_args, fnames):
    
    default_kwargs = {
        'interp': 2 + 1e-10,
        'nsim': 3000,
        'mcmode': 0,
        'sigmode': 0.2,
        'thres': 0.8,
        'nbin': 50        
    }
    
    params = { **default_kwargs, **input_args }
    
    interp = params['interp']
    nsim = params['nsim']
    mcmode = params['mcmode']
    sigmode = params['sigmode']
    thres = params['thres']
    nbin = params['nbin']
    
    return interp, nsim, mcmode, sigmode, thres, nbin



def set_pyzdcf(input_args, fnames):
    
    default_kwargs = {
        'nsim': 1000,
        'minpts': 0,
        'uniform_sampling': False,
        'omit_zero_lags': True,
        'sparse': 'auto',
        'prefix': 'zdcf',
        'run_plike': False,
        'plike_dir': None,
    }
    
    params = { **default_kwargs, **input_args }
    
    nsim = params['nsim']
    minpts = params['minpts']
    uniform_sampling = params['uniform_sampling']
    omit_zero_lags = params['omit_zero_lags']
    sparse = params['sparse']
    prefix = params['prefix']
    run_plike = params['run_plike']
    plike_dir = params['plike_dir']
    
    
    return nsim, minpts, uniform_sampling, omit_zero_lags, \
        sparse, prefix, run_plike, plike_dir




def set_javelin(input_args, fnames):
    
    default_kwargs = {
        'lagtobaseline': 0.3,
        'fixed': None,
        'p_fix': None,
        'subtract_mean': True,
        'nwalker': 100,
        'nburn': 100,
        'nchain': 100,
        'threads': 1,
        'output_chains': True,
        'output_burn': True,
        'output_logp': True,
        'nbin': 50,
        'metric': 'med',
        'together': False,
        'rm_type': 'spec'
    }
    
    params = { **default_kwargs, **input_args }
    
    lagtobaseline = params['lagtobaseline']
    fixed = params['fixed']
    p_fix = params['p_fix']
    subtract_mean = params['subtract_mean']
    nwalkers = params['nwalker']
    nburn = params['nburn']
    nchain = params['nchain']
    threads = params['threads']
    output_chains = params['output_chains']
    output_burn = params['output_burn']
    output_logp = params['output_logp']
    nbin = params['nbin']
    metric = params['metric']
    together = params['together']
    rm_type = params['rm_type']
    
    
    
    if (rm_type == 'phot') & (together):
        print('ERROR: JAVELIN cannot do phtotometric RM with more than two lines.')
        print('Setting together=False')
        together = False    
                
                
        
    if not together:
        
        if fixed is not None:
            if len(fixed) < len(fnames)-1:
            
                fixed_og = fixed
                p_fix_og = p_fix
                
                fixed = []
                p_fix = []
                for i in range(len(fnames)-1):
                    fixed.append(fixed_og)
                    p_fix.append(p_fix_og)
        
        else:
            fixed = np.full( len(fnames)-1, None )
            p_fix = np.full( len(fnames)-1, None )    
            
        assert len(fixed) == len(fnames)-1
        
    else:
        if fixed is not None:
            assert len(fixed) == 2 + 3*( len(fnames) - 1 )
                
                
    return lagtobaseline, fixed, p_fix, subtract_mean, \
        nwalkers, nburn, nchain, threads, output_chains, \
            output_burn, output_logp, nbin, metric, together, rm_type
            
            
def set_weighting(input_args):
    
    default_kwargs = {
        'gap_size': 30,
        'k': 2,
        'width': 15,
        'zoom': True
    }
    
    params = { **default_kwargs, **input_args }
    
    gap_size = params['gap_size']
    k = params['k']
    width = params['width']
    zoom = params['zoom']
    
    return gap_size, k, width, zoom
