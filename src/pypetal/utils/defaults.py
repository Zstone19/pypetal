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

    if isinstance(lag_bounds, list):
        if len(lag_bounds) == 1:
            assert isinstance(lag_bounds[0], list) or isinstance(lag_bounds[0], str) or ( lag_bounds[0] is None )

        elif len(lag_bounds) == 2:
            pass

        else:
            for i in range(len(lag_bounds)):
                if isinstance(lag_bounds[i], list):
                    assert ( len(lag_bounds[i]) == 2 ) | ( lag_bounds[i] is None ) | ( lag_bounds[i] == 'baseline' )
                else:
                    assert ( lag_bounds[i] is None ) | ( lag_bounds[i] == 'baseline' )


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
        if not isinstance(lag_bounds[0], str):
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


def set_pyccf(input_args):

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


def set_pyroa(input_args, nlc):

    default_kwargs = {
        'nchain': 20000,
        'nburn': 15000,
        'init_tau': None,
        'subtract_mean': True,
        'div_mean': False,
        'add_var': True,
        'delay_dist': True,
        'psi_types': 'Gaussian',
        'together': True,
        'objname': None,
        'prior_func': None
    }

    params = { **default_kwargs, **input_args }

    nchain = params['nchain']
    nburn = params['nburn']
    init_tau = params['init_tau']
    subtract_mean = params['subtract_mean']
    div_mean = params['div_mean']
    add_var = params['add_var']
    delay_dist = params['delay_dist']
    psi_types = params['psi_types']
    together = params['together']
    objname = params['objname']
    prior_func = params['prior_func']

    if init_tau is None:
        init_tau = [10.] * (nlc-1)

    if objname is None:
        objname = 'pyroa'

    if isinstance(psi_types, str):
        psi_types = [psi_types] * (nlc-1)

    if not together:
        if isinstance(delay_dist, bool):
            delay_dist = [delay_dist] * (nlc-1)
        if isinstance(add_var, bool):
            add_var = [add_var] * (nlc-1)


    return nchain, nburn, init_tau, subtract_mean, div_mean, \
            add_var, delay_dist, psi_types, together, objname, prior_func


def set_mica2(input_args):
    
    default_kwargs = {
        'type_tf': 'gaussian',
        'max_num_saves': 2000,
        'flag_uniform_var_params': False,
        'flag_uniform_tranfuns': False,
        'flag_trend': 0,
        'flag_lag_positivity': False,
        'flag_negative_resp': False,
        'number_components': [1,1],
        'width_limit': None,
        'flag_con_sys_err': False,
        'flag_line_sys_err': False,
        'type_lag_prior': 0,
        'lag_prior': None,
        'num_particles': 1,
        'thread_steps_factor': 1,
        'new_level_interval_factor': 1,
        'save_interval_factor': 1,
        'lam': 10,
        'beta': 100,
        'ptol': 0.1,
        'max_num_levels': 0,
        'together': False
    }
    
    params = { **default_kwargs, **input_args }

    type_tf = params['type_tf']
    max_num_saves = params['max_num_saves']
    flag_uniform_var_params = params['flag_uniform_var_params']
    flag_uniform_tranfuns = params['flag_uniform_tranfuns']
    flag_trend = params['flag_trend']
    flag_lag_positivity = params['flag_lag_positivity']
    flag_negative_resp = params['flag_negative_resp']
    number_components = params['number_components']
    width_limit = params['width_limit']
    flag_con_sys_err = params['flag_con_sys_err']
    flag_line_sys_err = params['flag_line_sys_err']
    type_lag_prior = params['type_lag_prior']
    lag_prior = params['lag_prior']
    num_particles = params['num_particles']
    thread_steps_factor = params['thread_steps_factor']
    new_level_interval_factor = params['new_level_interval_factor']
    save_interval_factor = params['save_interval_factor']
    lam = params['lam']
    beta = params['beta']
    ptol = params['ptol']
    max_num_levels = params['max_num_levels']
    together = params['together']
    
    return type_tf, max_num_saves, flag_uniform_var_params, flag_uniform_tranfuns, \
        flag_trend, flag_lag_positivity, flag_negative_resp, number_components, \
        width_limit, flag_con_sys_err, flag_line_sys_err, type_lag_prior, lag_prior, \
        num_particles, thread_steps_factor, new_level_interval_factor, save_interval_factor, \
        lam, beta, ptol, max_num_levels, together
        
    


def set_weighting(input_args):

    default_kwargs = {
        'gap_size': 30,
        'k': 2,
        'width': 15,
        'rel_height': .99,
        'zoom': True
    }

    params = { **default_kwargs, **input_args }

    gap_size = params['gap_size']
    k = params['k']
    width = params['width']
    rel_height = params['rel_height']
    zoom = params['zoom']

    return gap_size, k, width, rel_height, zoom
