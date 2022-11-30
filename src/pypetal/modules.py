import pypetal.utils as utils
import pypetal.plotting as plotting
from pypetal.formatting import write_data

import os

import numpy as np
from astropy.table import Table
import astropy.units as u
import matplotlib.pyplot as plt

#For individual functions, see utils.py
#For plotting tools, see plotting.py
#For the pipeline, see pipeline.py

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
    time_unit = general_kwargs['time_unit']
    lc_unit = general_kwargs['lc_unit']
        
    #--------------------------------------------------
    #Read kwargs
    
    default_kwargs = {
        'jitter' : True,
        'nsig' : 1,
        'nwalkers': 100,
        'nburn': 300,
        'nchain': 1000,
        'clip': np.full( len(line_fnames) + 1, True),
        'reject_data': np.full( len(line_fnames) + 1, True),
        'use_for_javelin': False
    }
    
    params = { **default_kwargs, **kwargs }
    
    jitter = params['jitter']
    nsig = params['nsig']
    nwalkers = params['nwalkers']
    nburn = params['nburn']
    nchain = params['nchain']
    clip = params['clip']
    reject_data = params['reject_data']
    use_for_javelin = params['use_for_javelin']
            
    if type(clip) is bool:
        clip = np.full( len(line_fnames) + 1, clip)
        
    if type(reject_data) is bool:
        reject_data = np.full( len(line_fnames) + 1, True)
        
    #--------------------------------------------------
    #Get units
    
    time_unit = str2unit(time_unit)
    lc_unit = str2unit(lc_unit) 
    
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
        """.format( jitter, nsig, nwalkers, nburn, 
                    nchain, clip_str, reject_data, use_for_javelin)
        
        print(txt_str)
    
    
    sigmas = []
    taus = []
    
    if jitter:
        jitters = []
    
    
    #Do continuum
    x, y, yerr = np.loadtxt( cont_fname, delimiter=',', unpack=True, usecols=[0,1,2] )

    res = utils.drw_flag( x*time_unit, y*lc_unit, yerr*lc_unit, 
                                nwalkers=nwalkers, nburn=nburn, nsamp=nchain,
                                nsig=nsig, jitter=jitter, clip=clip[0], 
                                target=line_names[0], 
                                fname=output_dir + line_names[0] + '/drw_rej/' + line_names[0] + '_drw_fit.pdf', plot=verbose)
    
    cont_mask = res['mask']
    
    sigmas.append( np.median(res['sigma']) )
    taus.append( np.median(res['tau']) )
    
    if jitter:
        jitters.append( np.median(res['jitter']) )
        
        
        
    #Write mask
    dat_fname = output_dir + line_names[0] + '/drw_rej/' + line_names[0] + '_mask.dat'
    write_data( cont_mask, dat_fname )
        
    #Write LC fit
    dat_fname = output_dir + line_names[0] + '/drw_rej/' + line_names[0] + '_drw_fit.dat'
    header = '#time,light curve,error'
    write_data( [ res['fit_x'], res['fit_y'], res['fit_err'] ], dat_fname, header )    
    
    #Write the MCMC chains
    dat_fname = output_dir + line_names[0] + '/drw_rej/' + line_names[0] + '_chain.dat'
    
    if jitter:
        header = '#sigma,tau,jitter'
        write_data( [ res['sigma'], res['tau'], res['jitter'] ], dat_fname, header )
    else:
        header = '#sigma,tau'
        write_data( [ res['sigma'], res['tau'] ], dat_fname, header )
    
    
    
    if reject_data[0]:
        x_new = x[~cont_mask]
        y_new = y[~cont_mask]
        yerr_new = yerr[~cont_mask]
        
        dat_fname = output_dir + line_names[0] + '_data.dat'
        write_data([x_new, y_new, yerr_new], dat_fname)
        
            
    line_masks = []
    
    for i in range(len(line_fnames)):
        x, y, yerr = np.loadtxt( line_fnames[i], delimiter=',', unpack=True, usecols=[0,1,2] )
        
        res = utils.drw_flag( x*time_unit, y*lc_unit, yerr*lc_unit, 
                                    nwalkers=nwalkers, nburn=nburn, nsamp=nchain,
                                    nsig=nsig, jitter=jitter, clip=clip[i+1],
                                    target=line_names[i+1],
                                    fname=output_dir + line_names[i+1] + '/drw_rej/' + line_names[i+1] + '_drw_fit.pdf', plot=verbose)
        
        sigmas.append( np.median(res['sigma']) )
        taus.append( np.median(res['tau']) )
        
        if jitter:
            jitters.append( np.median(res['jitter']) )
        
        line_mask = res['mask']
        line_masks.append(line_mask)
        
        #Write mask
        dat_fname = output_dir + line_names[i+1] + '/drw_rej/' + line_names[i+1] + '_mask.dat'
        write_data(line_mask, dat_fname)
                
        #Write LC fit
        dat_fname = output_dir + line_names[i+1] + '/drw_rej/' + line_names[i+1] + '_drw_fit.dat'
        header = '#time,light curve,error'
        write_data( [ res['fit_x'], res['fit_y'], res['fit_err'] ], dat_fname, header )        
        
        #Write the MCMC chains
        dat_fname = output_dir + line_names[i+1] + '/drw_rej/' + line_names[i+1] + '_chain.dat'

        if jitter:
            header = '#sigma,tau,jitter'
            write_data( [ res['sigma'], res['tau'], res['jitter'] ], dat_fname, header )
        else:
            header = '#sigma,tau'
            write_data( [ res['sigma'], res['tau'] ], dat_fname, header )        
        
        if reject_data[i+1]:
            x_new = x[~line_mask]
            y_new = y[~line_mask]
            yerr_new = yerr[~line_mask]
            
            dat_fname = output_dir + line_names[i+1] + '_data.dat'
            write_data( [x_new, y_new, yerr_new], dat_fname )
    
    tot_masks = []
    tot_masks.append(cont_mask)
    
    for i in range(len(line_masks)):
        tot_masks.append( line_masks[i] )
    
    output = {
        'masks': tot_masks,
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
    time_unit = general_kwargs['time_unit']
    lc_unit = general_kwargs['lc_unit']
    
    #--------------------------------------------------
    #Read kwargs
    
    default_kwargs = {
        'lag_bounds': np.full( len(line_fnames), None),
        'interp': None,
        'nsim': 1000,
        'mcmode': 0,
        'sigmode': .2,
        'thres': .8,
        'nbin': 50        
    }
        
    params = { **default_kwargs, **kwargs }
    
    lag_bounds = params['lag_bounds']
    interp = params['interp']
    nsim = params['nsim']
    mcmode = params['mcmode']
    sigmode = params['sigmode']
    thres = params['thres']
    nbin = params['nbin']

    
    if lag_bounds is None:
        lag_bounds = np.full( len(line_fnames), None)
        
    if len(lag_bounds) != len(line_fnames):    
        lag_bounds_og = lag_bounds
        lag_bounds = []
        
        for i in range(len(line_fnames)):
            lag_bounds.append( lag_bounds_og )

    
    #-------------------------------------------
    #Run pyCCF for each line
    
    if verbose:
        
        if interp is None:
            interp_str = 'mean/2'
        else:
            interp_str = interp
            
            
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
        """.format( lag_bounds_str, interp_str, nsim,
                    mcmode, sigmode, thres, nbin)
        
        print(txt_str)
        
    
    res_tot = []
    for i in range(len(line_fnames)):
        res = utils.get_pyccf_lags( cont_fname, line_fnames[i], 
                                   lag_bounds=lag_bounds[i], 
                                   interp=interp, nsim=nsim, mcmode=mcmode, 
                                   sigmode=sigmode, thres=thres)


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
                                    nbin=nbin, time_unit=time_unit, lc_unit=lc_unit,
                                    lc_names=[line_names[0], line_names[i+1]],
                                    plot_weights=False,
                                    fname = output_dir + line_names[i+1] + r'/pyccf/' + line_names[i+1] + '_ccf.pdf', show=verbose)

    
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
    time_unit = general_kwargs['time_unit']
    lc_unit = general_kwargs['lc_unit']
        
    #--------------------------------------------------
    #Read kwargs
    
    default_kwargs = {
        'nsim': 500,
        'minpts': 0,
        'uniform_sampling': False,
        'omit_zero_lags': True,
        'sparse': 'auto',
        'prefix': 'zdcf',
        'run_plike': False,
        'plike_dir': None,
        'lag_bounds': None
    }
        
    params = { **default_kwargs, **kwargs }
    
    nsim = params['nsim']
    minpts = params['minpts']
    uniform_sampling = params['uniform_sampling']
    omit_zero_lags = params['omit_zero_lags']
    sparse = params['sparse']
    prefix = params['prefix']
    run_plike = params['run_plike']
    plike_dir = params['plike_dir']
    lag_bounds = params['lag_bounds']
        
    if (lag_bounds == 'baseline') or (lag_bounds is None):
        lag_bounds = np.full( len(line_fnames), 'baseline' )
        

    len_lag_bounds = len(lag_bounds)
    if len_lag_bounds == 2:
        if not ( ( type(lag_bounds[0]) is type([]) ) or ( lag_bounds[0] == 'baseline' ) ):
            len_lag_bounds = 1
        else:
            len_lag_bounds = 2

    if len_lag_bounds != len(line_fnames):
        
        lag_bounds_og = lag_bounds
        lag_bounds = []

        for i in range(len(line_fnames)):
            lag_bounds.append(lag_bounds_og)
        
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
        
        
    if run_plike:
        
        if len(lag_bounds) > 2:
            lag_bounds_str = 'array'
        else:
            lag_bounds_str = lag_bounds
    
    input_dir += r'/'
    cont_fname_short = os.path.basename(cont_fname)
    line_fnames_short = [ os.path.basename(x) for x in line_fnames ]
        
    
    res_tot = []
    plike_tot = []
    
    for i in range(len(line_fnames)):
        dcf_df = utils.get_zdcf( input_dir, cont_fname_short, line_fnames_short[i], 
                                output_dir + line_names[i+1] + r'/pyzdcf/',
                                num_MC=nsim, minpts=minpts, 
                                uniform_sampling=uniform_sampling, omit_zero_lags=omit_zero_lags,
                                sparse=sparse, 
                                sep=',',
                                prefix=line_names[i+1] + '_' + prefix, 
                                verbose=verbose)
        
        res_tot.append(dcf_df)
        
        
        x1, y1, yerr1 = np.loadtxt( cont_fname, delimiter=',', unpack=True, usecols=[0,1,2] )
        x2, y2, yerr2 = np.loadtxt( line_fnames[i], delimiter=',', unpack=True, usecols=[0,1,2] )
        
        plike_dict = None
        if run_plike:
            
            if lag_bounds[i] == 'baseline':    
                lag_range = [ np.min([x1.min(), x2.min()]) ,
                               np.max([x1.max(), x2.max()]) ]

            else:
                lag_range = lag_bounds[i]
                
                        
            utils.run_plike( output_dir + line_names[i+1] + r'/pyzdcf/' + line_names[i+1] + '_' + prefix + '.dcf', lag_range, plike_dir,
                            verbose=verbose)   
                    
            plike_fname = plike_dir + 'plike.out'        
            os.rename( plike_fname, 
                      output_dir + line_names[i+1] + r'/pyzdcf/' + line_names[i+1] + '_plike.out' )
            
            plike_dat = Table.read( output_dir + line_names[i+1] + r'/pyzdcf/' + line_names[i+1] + '_plike.out', 
                                    format='ascii',
                                    names=['num', 'lag', 'r', '-dr', '+dr', 'likelihood'])
            
            
            file = open(output_dir + line_names[i+1] + r'/pyzdcf/' + line_names[i+1] + '_plike.out', 'r')
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
                                     dcf_df, plike_dict,
                                     time_unit=time_unit, lc_unit=lc_unit,
                                     lc_names=[line_names[0], line_names[i+1]],
                                     fname=output_dir + line_names[i+1] + r'/pyzdcf/' + line_names[i+1] + '_zdcf.pdf', 
                                     show=verbose)    
    
    return res_tot, plike_tot


#################################################################
############################# JAVELIN ###########################
#################################################################

def javelin_tot(cont_fname, line_fnames, line_names, output_dir, general_kwargs, kwargs):
    
    if line_fnames is str:
        line_fnames = [line_fnames]    
        line_names = [line_names]
        
    #--------------------------------------------------
    #Read general kwargs
    
    verbose = general_kwargs['verbose']
    time_unit = general_kwargs['time_unit']
    lc_unit = general_kwargs['lc_unit']    
                
    #--------------------------------------------------
    #Read kwargs
        
    default_kwargs = {
        'lagtobaseline': .3,
        'laglimit': 'baseline',
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
        'together': True,
        'rm_type': 'spec',
        'use_weights': False
    }
    
    
    params = { **default_kwargs, **kwargs }
    
    lagtobaseline = params['lagtobaseline']
    laglimit = params['laglimit']
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
    use_weights = params['use_weights']
    
    #--------------------------------------------------
    #Account for parameters if javelin['together'] = False
    if not together:
          
        if ( type(laglimit) is str ) | ( laglimit is None ):
            laglimit = np.full( len(line_fnames), laglimit )
            
        if ( len(line_fnames) == 1 ) & ( laglimit is not None ):
            laglimit = [laglimit]
    
        if fixed is None:
            fixed = np.full( len(line_fnames), fixed )
            p_fix = np.full( len(line_fnames), p_fix )
    
        if ( len(fixed) < len(line_fnames) ) or ( len(fixed) == 5 ):
            fixed_og = fixed
            p_fix_og = p_fix
            
            fixed = []
            p_fix = []
            for i in range(len(line_fnames)):
                fixed.append( fixed_og )
                p_fix.append( p_fix_og )

    #--------------------------------------------------
    
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
threads: {}
output_chains: {}
output_burn: {}
output_logp: {}
nbin: {}
metric: {}
together: {}
use_weights: {}
--------------------
        """.format( rm_type, lagtobaseline, laglimit, not (fixed is None), not (fixed is None),
                    subtract_mean, nwalkers, nburn, nchain, threads, output_chains,
                    output_burn, output_logp, nbin, metric, together, use_weights )
        
        print(txt_str)
    
    if together:
        res = utils.run_javelin(cont_fname, line_fnames, line_names, 
                                rm_type=rm_type,
                                output_dir=output_dir + r'javelin/', 
                                lagtobaseline=lagtobaseline, laglimit=laglimit,
                                fixed=fixed, p_fix=p_fix, subtract_mean=subtract_mean,
                                nwalkers=nwalkers, nburn=nburn, nchain=nchain, threads=threads,
                                output_chains=output_chains, output_burn=output_burn, output_logp=output_logp,
                                verbose=verbose)


        #Plot histograms
        fig, ax = plotting.plot_javelin_hist( res, fixed=fixed, nbin=nbin,
                                             time_unit=time_unit, lc_unit=lc_unit,
                                             plot_weights=use_weights, remove_fixed=False,
                                             fname= output_dir + 'javelin/javelin_histogram.pdf' )
        
        if verbose:
            plt.show()
            
        plt.cla()
        plt.clf()
        plt.close()
        
        
        if (fixed is None) | (fixed is np.ones( 2 + 3*len(line_fnames) )):
                    
            #Corner plot
            fig, ax = plotting.javelin_corner(res, plot_weights=use_weights, 
                                            fname= output_dir + 'javelin/javelin_corner.pdf' )
                
            plt.cla()
            plt.clf()
            plt.close()



        
        bestfit_model = utils.javelin_pred_lc( res['rmap_model'], 
                                            res['tot_dat'].jlist[0], res['tot_dat'].jlist[1:],
                                            nbin=nbin, prob_weight=use_weights,
                                            metric=metric)
        
        res['bestfit_model'] = bestfit_model
        
        
        #Plot model fits
        fig, ax = plotting.plot_javelin_bestfit(res, bestfit_model, time_unit=time_unit, lc_unit=lc_unit,
                                                fname= output_dir + 'javelin/javelin_bestfit.pdf' )

        if verbose:
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
            
            res = utils.run_javelin(cont_fname, line_fnames[i], names_i, 
                                    output_dir=output_dir + names_i[1] + r'/javelin/', 
                                    lagtobaseline=lagtobaseline, laglimit=laglimit[i],
                                    fixed=fixed[i], p_fix=p_fix[i], subtract_mean=subtract_mean,
                                    nwalkers=nwalkers, nburn=nburn, nchain=nchain, threads=threads,
                                    output_chains=output_chains, output_burn=output_burn, output_logp=output_logp,
                                    verbose=verbose)
            
            res_tot.append(res)
            
            #Plot histograms
            fig, ax = plotting.plot_javelin_hist( res, fixed=fixed[i], nbin=nbin,
                                                  time_unit=time_unit, lc_unit=lc_unit,
                                                  plot_weights=use_weights, remove_fixed=False,
                                                  fname= output_dir + line_names[i+1] + r'/javelin/javelin_histogram.pdf' )
            
            if verbose:
                plt.show()
                
            plt.cla()
            plt.clf()
            plt.close()


            
            #Corner plot
            fig, ax = plotting.javelin_corner(res, plot_weights=use_weights, 
                                              fname= output_dir + 'javelin/javelin_corner.pdf' )
                
            plt.cla()
            plt.clf()
            plt.close()
            
            
            bestfit_model = utils.javelin_pred_lc( res['rmap_model'], 
                                                res['tot_dat'].jlist[0], res['tot_dat'].jlist[1:],
                                                nbin=nbin, prob_weight=use_weights,
                                                metric=metric)
            
            
            #Plot model fits
            fig, ax = plotting.plot_javelin_bestfit(res, bestfit_model,
                                                    fname= output_dir + line_names[i+1] + '/javelin/javelin_bestfit.pdf' )

            if verbose:
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