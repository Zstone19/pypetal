import pypetal.utils as utils
import pypetal.plotting as plotting
from pypetal.petalio import write_data
import pypetal.defaults as defaults

import os
import shutil
import multiprocessing as mp
from functools import partial
import itertools

import numpy as np
from astropy.table import Table
import astropy.units as u
import matplotlib.pyplot as plt

#For individual functions, see utils.py
#For plotting tools, see plotting.py
#For the pipeline, see pipeline.py


#For multiprocessing
def mp_map(func, args, threads):
    n_inputs = len(args)
    
    if (threads > 1) & (n_inputs > 1):
        pool = mp.Pool(threads)
        res = pool.starmap(func, args)
        
        pool.close()
        pool.join()
        
    else:
        res = list(itertools.starmap(func, args))

    return res





#################################################################
######################## DRW REJECTION ##########################
#################################################################

def str2unit(x):
    
    if (x == 'Arbitrary Units'):
        unit = u.dimensionless_unscaled
    else:
        unit = u.Unit(x)
    
    return unit

def drw_rej_tot(cont_fname, line_fnames, line_names, output_dir, 
                general_kwargs, kwargs):


    #--------------------------------------------------
    #Read general kwargs
    
    verbose = general_kwargs['verbose']
    plot = general_kwargs['plot']
    time_unit = general_kwargs['time_unit']
    lc_unit = general_kwargs['lc_unit']
    threads = general_kwargs['threads']
        
    #--------------------------------------------------
    #Read kwargs
    
    jitter, nsig, nwalker, nburn, nchain, clip, \
        reject_data, use_for_javelin = defaults.set_drw_rej(kwargs, 
                                                            np.hstack([ [cont_fname], line_fnames ])
                                                            )
            
    #--------------------------------------------------
    #Get units
    
    time_unit = str2unit(time_unit)
    
    for i in range(len(lc_unit)):
        lc_unit[i] = str2unit(lc_unit[i]) 
    
    #--------------------------------------------------
    if verbose:
        
        if len(clip) > 2:
            clip_str = 'array'
        else:
            clip_str = clip
        
        txt_str = """
Performing DRW rejection
------------------------
jitter: {}
nsig: {}
nwalker: {}
nburn: {}
nchain: {}
clip: {}
reject_data: {}
use_for_javelin: {}
------------------------
        """.format( jitter, nsig, nwalker, nburn, 
                    nchain, clip_str, reject_data, use_for_javelin)
        
        print(txt_str)
    
    
    sigmas = []
    taus = []
    
    if jitter:
        jitters = []
    
    rejecting = np.any(reject_data)
    fnames = np.hstack([ [cont_fname], line_fnames ])
 
 
    drw_rej_func = partial( utils.drw_flag, nwalkers=nwalker, nburn=nburn, nsamp=nchain,
                            nsig=nsig, jitter=jitter, plot=plot)
       
       
       
    arg1 = []
    arg2 = []
    arg3 = []
    arg4 = line_names
    arg5 = []
    
    for i in range(len(fnames)):
        if reject_data[i]:
            x, y, yerr = np.loadtxt( fnames[i], delimiter=',', unpack=True, usecols=[0,1,2] )
    
            arg1.append( x*time_unit )
            arg2.append( y*lc_unit[i] )
            arg3.append( yerr*lc_unit[i] )
            arg5.append( output_dir + line_names[i] + '/drw_rej/' + line_names[i] + '_drw_fit.pdf' )
       
       
    args = list(zip(arg1, arg2, arg3, arg4, arg5))
    res_tot = mp_map(drw_rej_func, args, threads)
                        
    masks_tot = []
    n = 0
    for i in range(len(fnames)):
        
        if reject_data[i]:            
            res = res_tot[n]
            
            sigmas.append( np.median(res['sigma']) )
            taus.append( np.median(res['tau']) )
            
            if jitter:
                jitters.append( np.median(res['jitter']) )
            
            mask = res['mask']
            
            #Write mask
            dat_fname = output_dir + line_names[i] + '/drw_rej/' + line_names[i] + '_mask.dat'
            write_data(mask, dat_fname)
                    
            #Write LC fit
            dat_fname = output_dir + line_names[i] + '/drw_rej/' + line_names[i] + '_drw_fit.dat'
            header = '#time,light curve,error'
            write_data( [ res['fit_x'], res['fit_y'], res['fit_err'] ], dat_fname, header )        
            
            #Write the MCMC chains
            dat_fname = output_dir + line_names[i] + '/drw_rej/' + line_names[i] + '_chain.dat'

            if jitter:
                header = '#sigma,tau,jitter'
                write_data( [ res['sigma'], res['tau'], res['jitter'] ], dat_fname, header )
            else:
                header = '#sigma,tau'
                write_data( [ res['sigma'], res['tau'] ], dat_fname, header )        
        

            x_new = x[~mask]
            y_new = y[~mask]
            yerr_new = yerr[~mask]
            
            dat_fname = output_dir + 'processed_lcs/' + line_names[i] + '_data.dat'
            write_data( [x_new, y_new, yerr_new], dat_fname )

            
            n += 1 
        
        else:
            x, y, yerr = np.loadtxt( fnames[i], delimiter=',', unpack=True, usecols=[0,1,2] )
            
            mask = np.full( len(x), False )
            taus.append(None)
            sigmas.append(None)
            
            if jitter:
                jitters.append(None)
                
            if rejecting:
                dat_fname = output_dir + 'processed_lcs/' + line_names[i] + '_data.dat'
                write_data([x, y, yerr], dat_fname)       
                
            
        masks_tot.append(mask)
    
    output = {
        'masks': masks_tot,
        'reject_data': reject_data,
        'taus': taus,
        'sigmas': sigmas
    }
    
    if jitter:
        output['jitters'] = jitters
    
    return output


#################################################################
############################# pyCCF #############################
#################################################################

def pyccf_tot(cont_fname, line_fnames, line_names, output_dir, 
              general_kwargs, kwargs):
    
    if line_fnames is str:
        line_fnames = [line_fnames]  
        
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
    interp, nsim, mcmode, sigmode, thres, nbin = defaults.set_pyccf(kwargs, 
                                                                    np.hstack([ [cont_fname], line_fnames ])
                                                                    )

    
    #-------------------------------------------
    #Run pyCCF for each line
    
    if verbose:
            
        if len(lag_bounds) > 2:
            lag_bounds_str = 'array'
        else:
            lag_bounds_str = lag_bounds
        
        txt_str = """
Running pyCCF
-----------------
lag_bounds: {}
interp: {}
nsim: {}
mcmode: {}
sigmode: {}
thres: {}
nbin: {}
-----------------
        """.format( lag_bounds_str, interp, nsim,
                    mcmode, sigmode, thres, nbin)
        
        print(txt_str)
     
     
    res_tot = []
    for i in range(len(line_fnames)):
        res = utils.get_pyccf_lags( cont_fname, line_fnames[i], lag_bounds[i],
                                   interp=interp, nsim=nsim, 
                                   mcmode=mcmode, sigmode=sigmode, thres=thres,
                                   threads=threads, verbose=verbose)

        #Write CCF to file
        dat_fname = output_dir + line_names[i+1] + r'/pyccf/' + line_names[i+1] + '_ccf.dat'
        header = '#Lag,CCF'
        write_data( [ res['CCF_lags'], res['CCF'] ], dat_fname, header )        
        
        #Write CCCD, CCPD to file
        dat_fname = output_dir + line_names[i+1] + r'/pyccf/' + line_names[i+1] + '_ccf_dists.dat'
        header = '#CCCD,CCPD'
        write_data( [ res['CCCD_lags'], res['CCPD_lags'] ], dat_fname, header )        
        
        
        res['name'] = line_names[i]
    
        x1, y1, yerr1 = np.loadtxt( cont_fname, delimiter=',', unpack=True, usecols=[0,1,2] )
        x2, y2, yerr2 = np.loadtxt( line_fnames[i], delimiter=',', unpack=True, usecols=[0,1,2] )
        
        plotting.plot_pyccf_results(x1, y1, yerr1, x2, y2, yerr2,
                                    res['CCF_lags'], res['CCF'],
                                    res['CCCD_lags'], res['CCPD_lags'],
                                    nbin=nbin, time_unit=time_unit, lc_unit=[lc_unit[0], lc_unit[i+1]],
                                    lc_names=[line_names[0], line_names[i+1]],
                                    fname = output_dir + line_names[i+1] + r'/pyccf/' + line_names[i+1] + '_ccf.pdf', 
                                    show=plot)
    
    
        res_tot.append(res)
    
    
    return res_tot 


#################################################################
############################ pyZDCF #############################
#################################################################

def pyzdcf_tot(cont_fname, line_fnames, line_names, output_dir, 
               general_kwargs, kwargs):
    
    
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
            
    nsim, minpts, uniform_sampling, omit_zero_lags, \
        sparse, prefix, run_plike, plike_dir = defaults.set_pyzdcf(kwargs, 
                                                                   np.hstack([ [cont_fname], line_fnames ])
                                                                   )
        
    #-------------------------------------------
    if (run_plike) & (plike_dir is None):
        print('Error: plike_dir must be specified if run_plike=True')
        print('Skipping PLIKE')
        run_plike = False
    
    input_dir = os.path.dirname( os.path.realpath(cont_fname) )
    
    for i in range(len(line_fnames)):
        if input_dir != os.path.dirname( os.path.realpath(line_fnames[i]) ):
            print('ERROR: All light curve files must be in the same directory')
            return {}
            
    
    if verbose:
        
        txt_str = """
Running pyZDCF
----------------------
nsim: {}
minpts: {}
uniform_sampling: {}
omit_zero_lags: {}
sparse: {}
prefix: {}
run_plike: {}
plike_dir: {}
----------------------
        """.format( nsim, minpts, uniform_sampling, omit_zero_lags,
                    sparse, prefix, run_plike, plike_dir)
        
        print(txt_str)
    
    
    
    input_dir += r'/'
    cont_fname_short = os.path.basename(cont_fname)
    line_fnames_short = [ os.path.basename(x) for x in line_fnames ]
        
    
    pyzdcf_func = partial( utils.get_zdcf, num_MC=nsim, minpts=minpts, 
                           uniform_sampling=uniform_sampling, omit_zero_lags=omit_zero_lags,
                           sparse=sparse, sep=',', verbose=False)
    
    arg1 = np.full( len(line_fnames), input_dir )
    arg2 = np.full( len(line_fnames), cont_fname_short )
    arg3 = line_fnames_short
    arg4 = [ output_dir + x + r'/pyzdcf/' for x in line_names[1:] ]
    arg5 = [ x + '_' + prefix for x in line_names[1:] ]
    
    args = list(zip(arg1, arg2, arg3, arg4, arg5))
    res_tot = mp_map(pyzdcf_func, args, threads)
    
    plike_tot = []
    for i in range(len(line_fnames)):        
        
        x1, y1, yerr1 = np.loadtxt( cont_fname, delimiter=',', unpack=True, usecols=[0,1,2] )
        x2, y2, yerr2 = np.loadtxt( line_fnames[i], delimiter=',', unpack=True, usecols=[0,1,2] )
        
        plike_dict = None
        if run_plike:
                                                        
            utils.run_plike( arg4[i] + arg5[i] + '.dcf', lag_bounds[i], plike_dir,
                            verbose=verbose)   
                    
            plike_fname = plike_dir + 'plike.out'     
            output_fname = output_dir + line_names[i+1] + r'/pyzdcf/' + line_names[i+1] + '_plike.out'
            assert os.path.exists( os.path.dirname(output_fname) )
            
            shutil.move( plike_fname, output_fname )
            plike_dat = Table.read( output_fname, 
                                    format='ascii',
                                    names=['num', 'lag', 'r', '-dr', '+dr', 'likelihood'])
            
            
            file = open(output_fname, 'r')
            output_str = list(file)[-3:]
            
            ml_lag = float( output_str[1].split()[7] )
            ml_lag_err_hi = np.abs( float( output_str[1].split()[8] )  )
            ml_lag_err_lo = np.abs( float( output_str[1].split()[9] )  )
        
            plike_dict = {
                'output': plike_dat,
                'ML_lag': ml_lag,
                'ML_lag_err_hi': ml_lag_err_hi,
                'ML_lag_err_lo': ml_lag_err_lo
            }
  
            plike_tot.append( plike_dict )  
            
        plotting.plot_pyzdcf_results(x1, y1, yerr1, x2, y2, yerr2,
                                     res_tot[i], plike_dict,
                                     time_unit=time_unit, lc_unit=[lc_unit[0], lc_unit[i+1]],
                                     lc_names=[line_names[0], line_names[i+1]],
                                     fname=output_dir + line_names[i+1] + r'/pyzdcf/' + line_names[i+1] + '_zdcf.pdf', 
                                     show=plot)    
    
    return res_tot, plike_tot


#################################################################
############################# JAVELIN ###########################
#################################################################

def javelin_tot(cont_fname, line_fnames, line_names, output_dir, general_kwargs, kwargs):
    
    if line_fnames is str:
        line_fnames = [line_fnames]    
        
    #--------------------------------------------------
    #Read general kwargs
    
    verbose = general_kwargs['verbose']
    plot = general_kwargs['plot']
    time_unit = general_kwargs['time_unit']
    lc_unit = general_kwargs['lc_unit']   
    laglimit = general_kwargs['lag_bounds']
    threads = general_kwargs['threads']
                
    #--------------------------------------------------
    #Read kwargs
    
    lagtobaseline, fixed, p_fix, subtract_mean, \
        nwalkers, nburn, nchain, threads, output_chains, \
            output_burn, output_logp, nbin, metric, \
                together, rm_type = defaults.set_javelin( kwargs,
                                                         np.hstack([ [cont_fname], line_fnames ])
                                                         )
    
    #--------------------------------------------------
    #Account for parameters if javelin['together'] = False
    if not together:
          
        if ( type(laglimit) is str ) | ( laglimit is None ):
            laglimit = np.full( len(line_fnames), laglimit )                    

    #--------------------------------------------------
    
    if ( not isinstance(laglimit, str) ) & ( laglimit is not None ):
        if len(laglimit) > 2:
            laglimit_str = 'array'
        else:
            laglimit_str = laglimit
    else:
        laglimit_str = laglimit
    
    
    if verbose:
        txt_str = """
Running JAVELIN
--------------------
rm_type: {}
lagtobaseline: {}
laglimit: {}
fixed: {}
p_fix: {}
subtract_mean: {}
nwalker: {}
nburn: {}
nchain: {}
output_chains: {}
output_burn: {}
output_logp: {}
nbin: {}
metric: {}
together: {}
--------------------
        """.format( rm_type, lagtobaseline, laglimit_str, not (fixed is None), not (fixed is None),
                    subtract_mean, nwalkers, nburn, nchain, output_chains,
                    output_burn, output_logp, nbin, metric, together )
        
        print(txt_str)
    
    if together:
        res = utils.run_javelin(cont_fname, line_fnames, line_names, 
                                rm_type=rm_type,
                                output_dir=output_dir + r'javelin/', 
                                lagtobaseline=lagtobaseline, laglimit=laglimit,
                                fixed=fixed, p_fix=p_fix, subtract_mean=subtract_mean,
                                nwalkers=nwalkers, nburn=nburn, nchain=nchain, threads=threads,
                                output_chains=output_chains, output_burn=output_burn, output_logp=output_logp,
                                nbin=nbin, verbose=verbose, plot=plot)


        #Plot histograms
        fig, ax = plotting.plot_javelin_hist( res, fixed=fixed, nbin=nbin,
                                             time_unit=time_unit,
                                             remove_fixed=False,
                                             fname= output_dir + 'javelin/javelin_histogram.pdf' )
        
        if plot:
            plt.show()
            
        plt.cla()
        plt.clf()
        plt.close()
        
        
        if (fixed is None) | (fixed is np.ones( 2 + 3*len(line_fnames) )):
                    
            #Corner plot
            fig, ax = plotting.javelin_corner(res,
                                            fname= output_dir + 'javelin/javelin_corner.pdf' )
                
                
            if plot:
                plt.show()
                
            plt.cla()
            plt.clf()
            plt.close()



        
        bestfit_model = utils.javelin_pred_lc( res['rmap_model'], 
                                            res['tot_dat'].jlist[0], res['tot_dat'].jlist[1:],
                                            metric=metric)
        
        res['bestfit_model'] = bestfit_model
        
        
        #Plot model fits
        fig, ax = plotting.plot_javelin_bestfit(res, bestfit_model, time_unit=time_unit, lc_unit=lc_unit,
                                                fname= output_dir + 'javelin/javelin_bestfit.pdf' )

        if plot:
            plt.show()
        
        plt.cla()
        plt.clf()
        plt.close()




        #Write fits to light curves
        for i in range(len(line_names)):
            dat_fname = output_dir + r'javelin/' + line_names[i] + '_lc_fits.dat'
            dat = [ bestfit_model.jlist[i], 
                    bestfit_model.mlist[i] + bestfit_model.blist[i],
                    bestfit_model.elist[i] ]
            write_data( dat, dat_fname )            


        return res
    
    
    else:
        res_tot = []
        
        for i in range(len(line_fnames)):
            names_i = [line_names[0], line_names[i+1]]
            
            if isinstance(laglimit[i], str):
                input_laglimit = laglimit[i]
            else:
                input_laglimit = [laglimit[i]]
            
            res = utils.run_javelin(cont_fname, line_fnames[i], names_i, 
                                    output_dir=output_dir + names_i[1] + r'/javelin/', 
                                    lagtobaseline=lagtobaseline, laglimit=input_laglimit,
                                    fixed=fixed[i], p_fix=p_fix[i], subtract_mean=subtract_mean,
                                    nwalkers=nwalkers, nburn=nburn, nchain=nchain, threads=threads,
                                    output_chains=output_chains, output_burn=output_burn, output_logp=output_logp,
                                    verbose=verbose, plot=plot)
            
            res_tot.append(res)
            
            #Plot histograms
            fig, ax = plotting.plot_javelin_hist( res, fixed=fixed[i], nbin=nbin,
                                                  time_unit=time_unit,
                                                  remove_fixed=False,
                                                  fname= output_dir + line_names[i+1] + r'/javelin/javelin_histogram.pdf' )
            
            if plot:
                plt.show()
                
            plt.cla()
            plt.clf()
            plt.close()


            
            #Corner plot
            if (fixed[i] is None) | (fixed[i] is np.ones( 2 + 3*len(line_fnames) )):
                fig, ax = plotting.javelin_corner(res,
                                                fname= output_dir + line_names[i+1] + '/javelin/javelin_corner.pdf' )
                    
                if plot:
                    plt.show()
                    
                plt.cla()
                plt.clf()
                plt.close()
            
            
            bestfit_model = utils.javelin_pred_lc( res['rmap_model'], 
                                                res['tot_dat'].jlist[0], res['tot_dat'].jlist[1:],
                                                metric=metric)
            
            
            #Plot model fits
            fig, ax = plotting.plot_javelin_bestfit(res, bestfit_model, time_unit=time_unit, 
                                                    lc_unit=[lc_unit[0], lc_unit[i+1]],
                                                    fname= output_dir + line_names[i+1] + '/javelin/javelin_bestfit.pdf' )

            if plot:
                plt.show()
            
            plt.cla()
            plt.clf()
            plt.close()
            
            
            
            #Write fits to light curves        
            for j in range( bestfit_model.nlc ):
            
                if j == 0:
                    name = line_names[0]
                else:
                    name = line_names[i+1]
            
                dat_fname = output_dir + line_names[i+1] + r'/javelin/' + name + '_lc_fits.dat'
                dat = [ bestfit_model.jlist[j], 
                        bestfit_model.mlist[j] + bestfit_model.blist[j],
                        bestfit_model.elist[j] ]
                write_data( dat, dat_fname )            
            
        return res_tot