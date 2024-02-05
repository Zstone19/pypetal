import os
from functools import partial

import numpy as np

import pypetal.mica2.utils as utils
from pypetal.utils import defaults
from pypetal.utils.petalio import print_subheader, print_error

def mica2_tot(cont_fname, line_fnames, line_names, output_dir, general_kwargs, mica2_params, comm, rank):
    
    if line_fnames is str:
        line_fnames = [line_fnames]
        line_names = [line_names]
    
    
    #--------------------------------------------------
    #Read general kwargs

    verbose = general_kwargs['verbose']
    plot = general_kwargs['plot']
    time_unit = general_kwargs['time_unit']
    lc_unit = general_kwargs['lc_unit']
    lag_bounds = general_kwargs['lag_bounds']
    threads = general_kwargs['threads']

    
    #--------------------------------------------------
    #Read kwargs

    type_tf, max_num_saves, flag_uniform_var_params, flag_uniform_tranfuns, \
        flag_trend, flag_lag_posivity, flag_negative_resp, number_component, \
        width_limit, flag_con_sys_err, flag_line_sys_err, type_lag_prior, lag_prior, \
        num_particles, thread_steps_factor, new_level_interval_factor, save_interval_factor, \
        lam, beta, ptol, max_num_levels, together, no_order = defaults.set_mica2(mica2_params)
    
    #-------------------------------------------

    if verbose:
        if rank == 0:
            print_subheader('Running MICA2', 35, mica2_params)


    mica2_func = partial(utils.run_mica2, 
                         type_tf=type_tf, max_num_saves=max_num_saves, 
                         flag_uniform_var_params=flag_uniform_var_params,
                         flag_uniform_tranfuns=flag_uniform_tranfuns,
                         flag_trend=flag_trend, flag_lag_posivity=flag_lag_posivity,
                         flag_negative_resp=flag_negative_resp, number_component=number_component,
                         width_limit=width_limit, flag_con_sys_err=flag_con_sys_err,
                         flag_line_sys_err=flag_line_sys_err, type_lag_prior=type_lag_prior,
                         lag_prior=lag_prior, num_particles=num_particles,
                         thread_steps_factor=thread_steps_factor, new_level_interval_factor=new_level_interval_factor,
                         save_interval_factor=save_interval_factor, lam=lam, beta=beta, ptol=ptol,
                         max_num_levels=max_num_levels, 
                         comm=comm, rank=rank,
                         together=together, show=plot)

    cwd = os.getcwd()
    cont_lc = np.loadtxt( cont_fname, delimiter=',', usecols=[0,1,2] )

    if together:
        line_lcs = []
        for i in range(len(line_fnames)):
            line_lc = np.loadtxt( line_fnames[i], delimiter=',', usecols=[0,1,2] )
            line_lcs.append(line_lc)
        
        os.chdir(output_dir + 'mica2/')
        
        res_tot = mica2_func([os.getcwd()+'/'], line_names, 
                             cont_lc, line_lcs, lag_limit=lag_bounds[0])
        os.chdir(cwd)  

    else:
        
        if no_order:
            res_tot = []

            for i in range(len(line_fnames)):
                os.chdir(output_dir + line_names[i+1] + '/mica2/')
                
                line_lc = np.loadtxt( line_fnames[i], delimiter=',', usecols=[0,1,2] )
                res = mica2_func([os.getcwd()+'/'], [line_names[0], line_names[i+1]], 
                                cont_lc, [line_lc], lag_limit=lag_bounds[i])

                res_tot.append(res)
                os.chdir(cwd)
                
        else:
            lagmin = np.min([ x[0] for x in lag_bounds ])
            lagmax = np.max([ x[1] for x in lag_bounds ])
            
            line_lcs = []
            for i in range(len(line_fnames)):
                line_lc = np.loadtxt( line_fnames[i], delimiter=',', usecols=[0,1,2] )
                line_lcs.append(line_lc)
            
            
            outdirs = []
            for i in range(len(line_fnames)):
                outdirs.append(output_dir + line_names[i+1] + '/mica2/')

            os.chdir(output_dir + 'mica2/')
            res_tot = mica2_func(outdirs, line_names, 
                                cont_lc, line_lcs, lag_limit=[lagmin, lagmax])
            os.chdir(cwd)


    return res_tot