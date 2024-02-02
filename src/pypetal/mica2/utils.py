from mpi4py import MPI
import numpy as np
import pymica

import os





def run_mica2(outdir, line_names, cont_lc, line_lcs, 
              #MICA2
              type_tf='gaussian', max_num_saves=2000, 
              flag_uniform_var_params=False, flag_uniform_tranfuns=False,
              flag_trend=0, flag_lag_posivity=False,
              flag_negative_resp=False,
              lag_limit=[0,100], number_component=[1,1], width_limit=None,
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
                flag_uniform_var_params=flag_uniform_var_params, flag_uniform_tranfuns=flag_uniform_tranfuns,
                flag_trend=flag_trend, flag_lag_posivity=flag_lag_posivity,
                flag_negative_resp=flag_negative_resp,
                lag_limit=lag_limit, number_component=number_component, width_limit=width_limit,
                flag_con_sys_err=flag_con_sys_err, flag_line_sys_err=flag_line_sys_err,
                type_lag_prior=type_lag_prior, lag_prior=lag_prior,
                num_particles=num_particles, thread_steps_factor=thread_steps_factor,
                new_level_interval_factor=new_level_interval_factor, save_interval_factor=save_interval_factor,
                lam=lam, beta=beta, ptol=ptol, max_num_levels=max_num_levels)
    model.run()
    

    #Make plots
    model.plot_results(doshow=show)
    model.post_process()
    
    #Save data
    cwd = os.getcwd() 
    
    if type_tf == 'gaussian':
        typetf = 0
    else:
        typetf = 1
    
    get_mica2_data(outdir, line_names, cwd+'/data/', 'data_input.txt',
                   lag_limit[0], lag_limit[1], typetf,
                   np.max(number_component), flag_trend, flag_uniform_var_params,
                   flag_uniform_tranfuns, flag_negative_resp)

    return model









def get_mica2_data(outdir, line_names, fdir, fname, tau_upp, tau_low, typetf,
                   ngau, flagtrend, flagvar, flagtran, flagnegresp):

    sample = np.atleast_2d(np.loadtxt(fdir+"/data/posterior_sample1d.txt_%d"%ngau))
    data = np.loadtxt(fdir+fname)
    sall = np.loadtxt(fdir+"/data/pall.txt_%d"%ngau)

    if flagtrend > 0:
       trend = np.loadtxt(fdir+"/data/trend.txt_%d"%ngau)
       
       
    fp = open(fdir+fname)
    # read numbe of datasets
    line = fp.readline()
    nd = int(line[1:])
    if flagvar == 1:
        num_params_var = 3
    else:
        num_params_var = 3*nd
       
       
    # number of parameters for long-term trend
    nq = flagtrend + 1
    
    
    # read number of data points in each dataset
    nl = []
    for i in range(nd):
        line = fp.readline()
        ls = line[1:].split(":")
        ns = np.array([int(i) for i in ls])
        nl.append(ns)
    fp.close()
    
    
    # read number of points of reconstructions
    fp = open(fdir+"/data/pall.txt_%d"%ngau, "r")
    line = fp.readline()
    nl_rec = []
    for i in range(nd):
        line = fp.readline()
        ls = line[1:].split(":")
        ns = np.array([int(i) for i in ls])
        nl_rec.append(ns)
    fp.close()
    
    
    # assign index of cont data
    indx_con_data = []
    indx_con_rec = []
    indx_con_data.append(0)
    indx_con_rec.append(0)
    for i in range(1, nd):
        ns = nl[i-1]
        ns_rec = nl_rec[i-1]
        indx_con_data.append(np.sum(ns) + indx_con_data[i-1])
        indx_con_rec.append(np.sum(ns_rec) + indx_con_rec[i-1])
    
    
    # assign index of the parmaeter for the first line of each dataset 
    indx_line = []
    indx_line.append(num_params_var)
    for i in range(1, nd):
        if flagtran == 1:
            indx_line.append(num_params_var)
        else:
            indx_line.append(indx_line[i-1] + (len(nl[i-1])-1)*(1+ngau*3))
       
       
    # # print time lags, median, and 68.3% confidence limits
    # if flagnegresp == False:
    #     sample_lag = np.zeros(sample.shape[0])
    #     weight_lag = np.zeros(sample.shape[0])
        
    #     #loop over datasets
    #     for m in range(nd):

    #         #loop over lines
    #         ns = nl[m]
    #         for j in range(1, len(ns)):
    #             sample_lag[:] = 0.0
    #             weight_lag[:] = 0.0

    #             for k in range(ngau):
    #                 if flagnegresp == 0:  # no negative response
    #                     sample_lag[:] +=  sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] * np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
    #                     weight_lag[:] +=  np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
    #                 else:
    #                     sample_lag[:] +=  sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
    #                     weight_lag[:] +=  1.0


    #             lag, err1, err2 = np.quantile(sample_lag/weight_lag, q=(0.5, (1.0-0.683)/2.0, 1.0-(1.0-0.683)/2.0))
    #             err1 = lag-err1
    #             err2 = err2 - lag


    # dtau = tau_upp - tau_low 
    ntau = 500
    tran = np.zeros((sample.shape[0], ntau))
    # shift = 0.0
    

    idx_q = 0 # index for long-term trend parameters


    for m in range(nd):
        ns = nl[m]
        ns_rec = nl_rec[m]
        
        
        # Get reconstructed continuum
        sall_con0 = sall[indx_con_rec[m]:(indx_con_rec[m]+ns_rec[0]), :] 
        np.savetxt(outdir + 'cont_recon.dat', sall_con0, delimiter=',',
                   header=['time', 'flux', 'errlo', 'errhi'])
        
        # Get trend
        if flagtrend > 0:
            x = np.linspace(sall_con0[0, 0], sall_con0[-1, 0], 100)
            y = np.zeros(100)
            for j in range(nq):
                y += trend[idx_q + j, 0] * x**(j)
        
            np.savetxt( outdir + 'trend.dat', [x, y], delimiter=',',
                        header='time, trend')


        # set time lag range for Gaussian centers and centriods
        tau1 = 1.0e10
        tau2 = -1.0e10
        for j in range(1, len(ns)):      
            for k in range(ngau):
                tau1 = np.min((tau1, np.min(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1])))
                tau2 = np.max((tau2, np.max(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1])))



        # set time lag range for transfer function 
        tau1_tf = 1.0e10
        tau2_tf = -1.0e10
        for j in range(1, len(ns)):   
            for k in range(ngau):
                tau1_tf = np.min((tau1_tf, np.min(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                                    -np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]))))
                tau2_tf = np.max((tau2_tf, np.max(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                                    +np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]))))
        
        tau1_tf = np.min((tau_low, tau1_tf))
        tau2_tf = np.max((tau_upp, tau2_tf))



        for j in range(1, len(ns)):
            
            #Get reconstructed line
            sall_hb = sall[(indx_con_rec[m] + np.sum(ns_rec[:j])):(indx_con_rec[m] + np.sum(ns_rec[:j+1])), :]
            np.savetxt(sall_hb, outdir + "{}_recon.dat".format(line_names[j]), delimiter=',',
                       header=['time', 'flux', 'errlo', 'errhi'])

            #Get centers
            for k in range(ngau):
                cen = sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                np.savetxt( outdir + "{}_centers_{}.dat".format(line_names[j], k), cen, delimiter=',' )


            #Get centroids
            if flagnegresp == False:
                
                cent = np.zeros(sample.shape[0])
                norm = np.zeros(sample.shape[0])
                for k in range(ngau):

                    if flagnegresp == 0: # no negative responses
                        norm += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
                        cent += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]) \
                                * sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                    else:
                        norm += 1.0
                        cent += sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
            
                    np.savetxt(outdir + "{}_centroids_{}.dat".format(line_names[j], k), cent, delimiter=',' )


        #Get transfer function
        tau = np.linspace(tau1_tf, tau2_tf, ntau)
        tran[:, :] = 0.0
        
        if typetf == 0: # gaussian
            for i in range(sample.shape[0]):
                # loop over gaussians
                for k in range(ngau):

                    if flagnegresp == 0:
                        amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
                    else:
                        amp =      sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]

                        cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                        sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
                        tran[i, :] += amp/sig * np.exp(-0.5*(tau - cen)**2/sig**2)

        else: #tophats
            
            for i in range(sample.shape[0]):
            # loop over tophats
                for k in range(ngau):

                    if flagnegresp == 0:
                        amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
                    else:
                        amp =      sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]


                    cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                    sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
                    
                    tran[i, :] += amp/sig/2.0 *(np.heaviside(sig-np.abs(tau-cen), 1.0))


        tran_best = np.percentile(tran, 50.0, axis=0)
        tran1 = np.percentile(tran, (100.0-68.3)/2.0, axis=0)
        tran2 = np.percentile(tran, 100.0-(100.0-68.3)/2.0, axis=0)

        errlo = tran_best - tran1
        errhi = tran2 - tran_best
        np.savetxt(outdir + "{}_transfunc.dat".format(line_names[j]), [tau, tran_best, errlo, errhi], delimiter=',',
                   header='tau, tf, errlo, errhi')

    return