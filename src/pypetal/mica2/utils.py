from mpi4py import MPI
import numpy as np
import pymica


def run_mica2(cont_lc, line_lcs, 
              #MICA2
              type_tf='gaussian', max_num_saves=2000, 
              flag_uniform_var_params=False, flag_uniform_transfuns=False,
              flag_trend=0, flag_lag_positivity=False,
              flag_negative_resp=False,
              lag_limit=[0,100], number_components=[1,1], width_limit=None,
              flag_con_sys_err=False, flag_line_sys_err=False,
              type_lag_prior=0, lag_prior=None, 
              #CDNEST
              num_particles=1, thread_steps_factor=1,
              new_level_interval_factor=1, save_interval_factor=1,
              lam=10, beta=100, ptol=0.1, max_num_levels=0,
              together=False, show=False):


    # initiate MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        
        data_input = {}
        
        if together:
            data_input["set1"] = [cont_lc]
            for i in range(len(line_lcs)):
                data_input["set1"].append(line_lcs[i])
            
        else:
            for i in range(len(line_lcs)):
                data_input["set{}".format(i+1)] = [cont_lc, line_lcs[i]]

    else:
        data_input = None


    data_input = comm.bcast(data_input, root=0)


    #Run MICA2
    model = pymica.gmodel()
    model.setup(data=data_input, type_tf=type_tf, max_num_saves=max_num_saves,
                flag_uniform_var_params=flag_uniform_var_params, flag_uniform_transfuns=flag_uniform_transfuns,
                flag_trend=flag_trend, flag_lag_positivity=flag_lag_positivity,
                flag_negative_resp=flag_negative_resp,
                lag_limit=lag_limit, number_components=number_components, width_limit=width_limit,
                flag_con_sys_err=flag_con_sys_err, flag_line_sys_err=flag_line_sys_err,
                type_lag_prior=type_lag_prior, lag_prior=lag_prior,
                num_particles=num_particles, thread_steps_factor=thread_steps_factor,
                new_level_interval_factor=new_level_interval_factor, save_interval_factor=save_interval_factor,
                lam=lam, beta=beta, ptol=ptol, max_num_levels=max_num_levels)
    model.run()
    

    #Make plots
    model.plot_results(doshow=show)
    model.post_process(doshow=show)

    return model

