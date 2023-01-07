import os
import subprocess
import time
import glob

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.units as u

from pyzdcf import pyzdcf

from javelin.zylc import get_data
from javelin.lcio import writelc
from javelin.lcmodel import Cont_Model, Rmap_Model, Pmap_Model, SPmap_Model

import pypetal.drw_funcs as drw
import pypetal.pyccf as pyccf


##############################################################
####################### DRW REJECTION ########################
##############################################################

def drw_flag(times, data, error,
             target=None, fname=None, 
             nwalkers=32, nburn=300, nsamp=1000,
             nsig=1, jitter=True, clip=True,
             plot=True):
    
    
    """Fits the light curve to a DRW model using celerite with MCMC sampling from emcee. 
    celerite will then predict the light curve with the best fit parameters. Points lying above
    or below the mean prediction by ``nsig`` standard deviations will be flagged as outliers in
    an output mask.


    Parameters
    ----------
    
    times : array_like
        Times for the light curve.
        
    data : array_like
        Data for the light curve.
    
    error : array_like
        Uncertainty in the light curve.
        
    nwalkers : int, optional
        Number of walkers for the MCMC sampler. Default is 32.
        
    nburn : int, optional
        Number of burn-in steps for the MCMC sampler. Default is 300.
        
    nsamp : int, optional
        Number of samples for the MCMC sampler. Default is 1000.
        
    nsig : int, optional
        Number of standard deviations above or below the mean prediction to flag as an outlier.
        Default is 1.
        
    jitter : bool, optional
        If True, will fit for a noise (jitter) term in the data. Default is ``True``.
        
    clip : bool, optional
        If true, will clip data points which are too close in time to another data point. Specifically, 
        will remove data points that are closer than 1e-8 days apart. Default is ``True``.
        
    target : str, optional
        The name of the target. Default is ``None``.
        
    fname : str, optional
        The name of the file to output the resulting plot. If ``None``, the plot will
        not be saved. Default is ``None``.
        
    plot : bool, optional
        If true, the plot will be shown. Default is ``True``.
        


    Returns
    --------
        
    res : dict 
        The output of the DRW fitting and rejection. The keys are:
        
        * mask : array_like
            The DRW rejection mask. Values marked as True are outliers.
            
        * tau : array_like
            The DRW tau chain from the MCMC sampler.
            
        * sigma : array_like
            The DRW sigma chain from the MCMC sampler.
            
        * fit_x : array_like
            The times for the predicted light curve.
            
        * fit_y : array_like
            The predicted light curve.
            
        * fit_yerr : array_like
            The uncertainty in the predicted light curve.
            
        * jitter : array_like
            The jitter chain from the MCMC sampler. Only present if ``jitter`` is True. 
    

    .. note:: The input times, data, and error must be astropy.units.Quantity objects.
    
    """
    
    assert times.shape == data.shape == error.shape
    assert len(times.shape) == 1
    
    if type(times[0]) == u.Quantity:
        time_unit = times[0].unit
    if type(data[0]) == u.Quantity:
        data_unit = data[0].unit
    
    flag_mask = np.zeros( times.shape, dtype=bool )
    
    tau_err = np.zeros( 2 )
    sigma_err = np.zeros( 2 )
    jitter_err = np.zeros( 2 )


    #Assume 'times' and 'fluxes' have the same shape
    
    #Fit to the DRW model   
    samples, gp, statuses = drw.MCMC_fit(times, data, error,
                                         nwalkers=nwalkers, nburn=nburn, nsamp=nsamp,
                                         jitter=jitter, clip=clip)


    fig, ax = drw.plot_outcome(times, data, error, samples, gp, data_unit, 
                                nsig=nsig, target=target, show_mean=True, 
                                filename=fname, jitter=jitter, show=plot)

    plt.cla()
    plt.clf()
    plt.close()

    tau_vals = 1/np.exp(samples[:, 1])
    sig_vals = np.sqrt( np.exp(samples[:, 0])/2 )
    
    tau_med = np.median(tau_vals)
    sigma_med = np.median(sig_vals)
    
    if jitter:
        jitter_vals = np.exp(samples[:, 2])
    
    
    #Get fit to light curve data
    baseline = times[-1] - times[0]
    extra_t = int(baseline.value//10)
    
    t = np.linspace( times[0].value - extra_t, times[-1].value + extra_t, 1000 )
    mu, var = gp.predict(data.value, t, return_var=True)
    std = np.sqrt(var)
    
    mu_flag, var_flag = gp.predict(data.value, times.value, return_var=True)
    std_flag = np.sqrt(var_flag)
        
        
    #Reject if data is beyond nsig*sig of fit mean
    flag_mask = np.abs(data.value - mu_flag) > nsig*std_flag    
        
    res = {
        'mask': flag_mask,
        'tau' : tau_vals,
        'sigma': sig_vals,
        'fit_x': t,
        'fit_y': mu,
        'fit_err': std
    }
    
    if jitter:
        res['jitter'] = jitter_vals
        
    return res


##############################################################
######################### ZDCF ###############################
##############################################################

def run_plike(dcf_fname, lag_bounds, plike_dir, verbose=False):
    
    """
    Runs the PLIKE algorithm to compute the maximum likelihood peak of the ZDCF.
    Must have PLIKE (https://ui.adsabs.harvard.edu/abs/2013arXiv1302.1508A/abstract).
    Will output a file containing the PLIKE results ('plike.out') in ``plike_dir``.
    
    
    Parameters
    ----------
    
    plike_dir : str
            Path to the directory with the PLIKE executable.
            
    lag_bounds : (2,) array_like
            Lower and upper bounds of lags to search for ML peak.
            
    dcf_name :str
            Path to the ZDCF file, usually ouptput by pyZDCF.
            
    verbose : bool, optional
            If True, will read output of PLIKE. Default is False.
        
        
    Returns
    -------
    
    """
    
    cwd = os.getcwd()    
    os.chdir(plike_dir)
    
    #Make sure dcf file exists
    assert dcf_fname in glob.glob( plike_dir + '*' )
    
    if verbose:
        print('Executing PLIKE')
    
    exec_str =  r"./plike <<< $'" + dcf_fname + r"\n" + str(lag_bounds[0]) + r"\n" +  str(lag_bounds[1]) + r"'"
    res = subprocess.check_output(exec_str, shell=True)
    
    #Make sure plike.out exists
    while plike_dir + 'plike.out' not in glob.glob( plike_dir + '*' ):
        time.sleep(1)
    
    os.chdir(cwd)
        
    return


def get_zdcf(input_dir, fname1, fname2, out_dir, prefix='zdcf', num_MC=500, minpts=0, 
             uniform_sampling=False, autocf=False, omit_zero_lags=True, 
             sparse='auto', sep=',', verbose=False):
    
    """
    Runs the pyZDCF algorithm to compute the Z-Transformed Discrete Correlation Function (ZDCF).
    For more information on the algorithm, see pyZDCF (https://pyzdcf.readthedocs.io/).
    The algorithm will take in two light curves and output a file containing the ZDCF ('zdcf.dcf')
    in the specified output directory. 
    
    
    Parameters
    ----------
    
    input_dir : str
            Path to the directory containing the light curves
   
    fname1 : str
            Name of the first light curve file
   
    fname2 : str
            Name of the second light curve file
   
    out_dir : str
            Path to the directory to place the output ZDCF file    
   
    num_MC : float, optional
            The number of Monte Carlo simulations to run. Default is 500.
   
    minpts : int, optional
            The minimum number of points to use in each bin when computing the ZDCF. 
            Must be larger than 11. If set to 0, it will be set to 11. Default is 0
            
    uniform_sampling: bool, optional
            If True, the light curves will be assumed to be uniformly sampled.
            Default is ``False``.
            
    autocf : bool, optional
            If True, the auto-correlation function for the first light curve will be computed.
            If False, the ZDCF will be computed between the light curves.
            Default is ``False``.
   
    omit_zero_lags : bool, optional
            If True, will omit the points with zero lags when computing the ZDCF.
            Default is ``True``.
            
    sparse : (bool, str), optional
            Determines whether to use a sparse matrix implementation for reduced RAM usage.
            This feature is suitable for longer light curves (> 3000 data points). If True, will
            use sparse matrix implementation. If set to 'auto', will use sparse matrix implementation
            if there are more than 3000 data points per light curve. Default is 'auto'
   
    sep : str, optional
            The delimiter used in the light curve files. Default is ',' for CSV files.            

    verbose : bool, optional
            If ``True``, will output progress of pyZDCF. Default is ``False``.



    Returns
    -------

    dcf_df : pandas.DataFrame
        A pandas DataFrame object containing the ZDCF and its errors. The 
        columns within the DataFrame match the order of columns within the 
        output ZDCF file 'zdcf.dcf'.
    
    """
    
    
    #Files need to be csv
    params = dict(        
        autocf = autocf,
        prefix = prefix,
        uniform_sampling = uniform_sampling,
        omit_zero_lags = True,
        minpts = 0,
        num_MC = num_MC,
        lc1_name = fname1,
        lc2_name = fname2
    )
    
    dcf_df = pyzdcf(input_dir=input_dir, output_dir=out_dir, intr=False,
                    verbose=verbose, parameters=params, sep=sep, sparse=sparse )

    return dcf_df
    
    
    
def get_zdcf_ml_lag(fname1, fname2, temp_dir, plike_dir, lag_bounds, sep=','):
    
    input_dir = os.path.dirname( os.path.realpath(fname1) )
    assert os.path.dirname( os.path.realpath(fname2) ) == input_dir
    
    input_dir += r'/'
    
    lc_name1 = os.path.basename(fname1)
    lc_name2 = os.path.basename(fname2)
    
    dcf_df = get_zdcf(input_dir, lc_name1, lc_name2, temp_dir, sep=sep)    
    dcf_fname = os.path.join(temp_dir, 'zdcf.dcf')

    run_plike(dcf_fname, lag_bounds, plike_dir)
    plike_fname = os.path.join(plike_dir, 'plike.out')
    
    plike_dat =  Table.read( plike_fname, format='ascii',
                             names=['num', 'lag', 'r', '-dr', '+dr', 'likelihood'])
    
    
    #Get peak likelihood lag
    file = open(plike_fname, 'r')
    output_str = list(file)[-3:]
    
    ml_lag = float( output_str[1].split()[7] )
    ml_lag_err_hi = np.abs( float( output_str[1].split()[8] )  )
    ml_lag_err_lo = np.abs( float( output_str[1].split()[9] )  )
    
    
    return {
        'output': plike_dat,
        'ML_lag': ml_lag,
        'ML_lag_err_hi': ml_lag_err_hi,
        'ML_lag_err_lo': ml_lag_err_lo
    }
    
    
##############################################################
######################### pyCCF ##############################
##############################################################
    
    
def get_pyccf_lags(fname1, fname2, 
                   lag_bounds=None, interp=None,
                   nsim=1000, mcmode=0, sigmode=.2, thres=.8,
                   threads=1, verbose=False):
    
    """Obtain time lags from the correlation between two light curves using the pyCCF algorithm: https://ui.adsabs.harvard.edu/abs/2018ascl.soft05032S/abstract.
    In short, pyCCF creates the Interpolated Cross-Correlation Function (ICCF) between two light curves, 
    and then samples the peak and centroid of the ICCF using Monte Carlo simulations with 
    Flux Randomization (FR) and/or Random Subset Resampling (RSS). These Monte Carlo simulations are then
    used to create the cross-correlation centroid distribution (CCCD) and peak distribution (CCPD). 



    Parameters
    ----------    

    fname1 : str
        Path to the first light curve file 
            
    fname2 : str
        Path to the second light curve file 
            
    file_fmt : str, optional
        Format of the file. Default is 'csv'.
            
    lag_bounds : (2,) array_like, optional
        The bounds of times to search for the lag. The first element is the minimum lag and the second is the maximum.
        If set to ``None", the lag bounds will be set to (-baseline, baseline). The default is ``None``.
    
    interp : float, optional
        The interval with which pyCCF will interpolate the ligh curves to form the ICCF. This value must be 
        shorter than the average cadence of the ligh curves. Setting this value too low can introduce noise.
        If set to ``None``, will be set to half of the average cadence of the light curves. The default is ``None``.
    
    nsim : int, optional
        The number of Monte Carlo simulations to run. The default is 1000.
    
    mcmode : int, optional
        The type of resampling to do for the Monte Carlo Simulations. 0 performs both FR and RSS, 1 performs FR, and 2 performs RSS.
        The default is 0.
    
    sigmode : float, optional
        The threshold for considering a measurement in the ICCF significant when computing peaks and centroids. Must be within the 
        interval (0,1). All peaks and centroids with correlation coefficient r_max <= sigmode will be considered as "failed". 
        If set to 0, will exclude all peaks based on a p-value significance test (see pyCCF documentation). The default is 0.2. 
    
    thres : float, optional
        The lower limit of correlation coefficient used when calculating the centroid of the ICCF. Must be within 
        the interval (0,1). The default is 0.8.
            

    .. note:: Both light curve files must be in CSV format with the following columns in order: time, value, uncertainty



    Returns
    -------
    
    result: dict 
        Dict of output results, containing:
        
        * 'CCF' : (N,) array_like
            The ICCF
            
        * 'CCF_lags' : (N,) array_like
            The lags corresponding to the ICCF

        * 'centroid' : float
            The median of the CCCD

        * 'centroid_err_hi' : float
            The upper error of the centroid

        * 'centroid_err_lo' : float
            The lower error of the centroid

        * 'peak' : float
            The median of the CCPD

        * 'peak_err_hi' : float
            The upper error of the peak

        * 'peak_err_lo' : float
            The lower error of the peak

        * 'CCCD_lags' : (nsim,) array_like
            The CCCD

        * 'CCPD_lags' : (nsim,) array_like
            The CCPD
                
    """
    
    #Files need to be csv
    x1, y1, yerr1 = np.loadtxt( fname1, delimiter=',', unpack=True, usecols=[0,1,2] )
    x2, y2, yerr2 = np.loadtxt( fname2, delimiter=',', unpack=True, usecols=[0,1,2] )
    
    #Sort by time
    sort_ind = np.argsort(x1)
    x1 = x1[sort_ind]
    y1 = y1[sort_ind]
    yerr1 = yerr1[sort_ind]
    
    sort_ind = np.argsort(x2)
    x2 = x2[sort_ind]
    y2 = y2[sort_ind]
    yerr2 = yerr2[sort_ind]
    
    #Get mean cadence and baseline
    diff1 = np.diff(x1)
    diff2 = np.diff(x2)
    mean_diff = np.mean( [ np.mean(diff1), np.mean(diff2) ] )
    baseline = np.max([ x1[-1], x2[-1] ]) - np.min([ x1[0], x2[0] ])

    #Define pyCCF parameters
    if lag_bounds is None:
        lag_bounds = [-baseline, baseline]
      
    if interp is None:   
       interp = mean_diff/2


    if np.any( np.abs(lag_bounds) > baseline ):
        print('Lag bounds are larger than baseline')

    
    #Run algorithm
    _, status_peak, _, status_centroid, ccf_pack, \
        _, status_rval, _ = pyccf.peakcent(x1, y1, x2, y2, lag_bounds[0],
                                           lag_bounds[1], interp)
        
    # assert status_peak == 1
    # assert status_centroid == 1
    
    tlags_peak, tlags_centroid, nsuccess_peak, nfail_peak, \
        nsuccess_centroid, nfail_centroid, max_rvals, nfail_rvals, \
            pvals = pyccf.xcor_mc(x1, y1, np.abs(yerr1), x2, y2, 
                                np.abs(yerr2), lag_bounds[0], lag_bounds[1], interp, 
                                nsim=nsim, mcmode=mcmode, sigmode=sigmode,
                                thres=thres, threads=threads, verbose=verbose)

    ccf = ccf_pack[0]
    ccf_lags = ccf_pack[1]
    
    #Get peak and centroid values/uncertainties        
    cent_med = np.median( tlags_centroid )
    cent_hi = np.percentile( tlags_centroid, 84 )
    cent_lo = np.percentile( tlags_centroid, 16 )

    peak_med = np.median( tlags_peak )
    peak_hi = np.percentile( tlags_peak, 84 )
    peak_lo = np.percentile( tlags_peak, 16 )   
    
    
    return {
        'CCF': ccf,
        'CCF_lags': ccf_lags,
        'centroid': cent_med,
        'centroid_err_hi': cent_hi - cent_med,
        'centroid_err_lo': cent_med - cent_lo,
        'peak': peak_med,
        'peak_err_hi': peak_hi - peak_med,
        'peak_err_lo': peak_med - peak_lo,
        'CCCD_lags': tlags_centroid,
        'CCPD_lags': tlags_peak
    }                                                                


############################################################################################################
############################################## JAVELIN #####################################################
############################################################################################################


def get_javelin_filenames(output_chain, output_burn, output_logp, prefix, output_dir=None):
    tot_fnames = np.full(3, None)
    
    if output_chain:
        tot_fnames[0] = 'chain_' + prefix + '.txt'
    
    if output_burn:
        tot_fnames[1] = 'burn_' + prefix + '.txt'
    
    if output_logp:
        tot_fnames[2] = 'logp_' + prefix + '.txt'
        
    if output_dir is not None:
        if output_dir[-1] != r'/':
            output_dir += r'/'
        
        for i in range(len(tot_fnames)):
            
            if tot_fnames[i] is not None:
                tot_fnames[i] = output_dir + tot_fnames[i]
        
    return tot_fnames



def run_javelin(cont_fname, line_fnames, line_names, 
                rm_type='spec',
                lagtobaseline=0.3, laglimit='baseline',
                fixed=None, p_fix=None, subtract_mean=True,
                nwalkers=100, nburn=100, nchain=100, threads=1, output_chains=False,
                output_burn=False, output_logp=False, output_dir=None, 
                nbin=50, verbose=False, plot=False):
    
    """Run JAVELIN on a set of light curves.
    
    
    
    Parameters
    ----------
    
    cont_fname : str
        The filename of the continuum light curve.
    
    line_fnames: str, list of str
        The filename(s) of the line light curve(s).
        
    rm_type : str, optional
        The type of analysis (and JAVELIN model) to use. May either be 'spec for 
        spectroscopic RM or 'phot' for photometric RM. Default is 'spec'.
        
    lagtobaseline : float, optional
        JAVELIN will use a log prior on the lag and penalize lags larger x*baseline,
        where x is the input value and baseline is the baseline of the light curves.
        Default is 0.3.
        
    laglimit : (2,) list of floats, str, optional
         The range of lags to search for the lag between light curves in the following form:
         [lower, upper]. If set to 'baseline', the range will be set to [-baseline, baseline].
         Default is 'baseline'.
         
    fixed : None or list of floats, optional
        An array defining which parameters are fixed and which are allowed to vary. The length of 
        the array must match the number of parameters ( 2 + 3*len(line_fnames) ). In the array, 
        1 corresponds to a variable parameter and 0 corresponds to a fixed parameter. If set to ``None``, 
        all parameters will be allowed to vary. Default is ``None``.
        
    p_fix : None or list of floats, optional
        An array defining the values of the fixed parameters. The length of the array must match the
        number of parameters and of the ``fixed" array. If set to ``None``, all parameters will be allowed
        to vary. Must be defined if ``fixed`` is. Default is ``None``.
        
    subtract_mean : bool, optional
        If True, will subtract the mean of the light curves before analysis. Default is ``True``.
        
    nwalkers : int, optional
        The number of walkers to use in the MCMC. Default is 100.
    
    nburn : int, optional
        The number of burn-in steps to use in the MCMC. Default is 100.
        
    nchain : int, optional
        The number of steps to use in the MCMC. Default is 100.
        
    threads : int, optional
        The number of parallel threads to use in the MCMC. Default is 1.    
        
    output_chains : bool, optional
        If ``True``, will output the MCMC chains to a file. Default is ``False``.
        
    output_burn : bool, optional
        If ``True``, will output the MCMC burn-in chains to a file. Default is ``False``.
        
    output_logp : bool, optional
        If ``True``, will output the MCMC log-probability values to a file. Default is ``False``.
         
         
        
    Returns
    -------

    res : dict
        A dictionary containing the results of the JAVELIN analysis. Has the following keys:
        
        * cont_hpd : (3,2) array_like
            The highest posterior density (HPD) interval for the continuum light curve DRW fits.

        * tau : array_like
            The MC chain for the DRW tau parameter.
            
        * sigma : array_like
            The MC chain for the DRW sigma parameter.
            
        * tophat_params : array_like
            The MC chains for the tophat parameters (lag, width, scale), 3 for each light curve.
            
        * hpd : array_like
            The HPD intervals for the DRW and tophat parameters for each light curve.
        
        * cont_model : javelin.lcmodel.Cont_Model
            The continuum model object.
            
        * rmap_model : javelin.lcmodel.Rmap_Model or javelin.lcmodel.Pmap_Model
            The RM model object.
            
        * cont_dat : javelin.zylc.LightCurve
            The continuum light curve object.
            
        * tot_dat : javelin.zylc.LightCurve
            The object containing all light curves (continuum + lines).
            
    """
    
    
    total_fnames = np.hstack( [cont_fname, line_fnames] )
    total_names = line_names

    #Deal with csv files    
    xi, yi, erri = np.loadtxt(cont_fname, delimiter=',', unpack=True, usecols=[0,1,2])                
    writelc( [[xi, yi, erri]] , output_dir + 'cont_lcfile.dat' )
    
    tot_lc_arr = []
    for i in range(len(total_fnames)):
        xi, yi, erri = np.loadtxt(total_fnames[i], delimiter=',', unpack=True, usecols=[0,1,2])            
        tot_lc_arr.append( [xi, yi, erri] )

    writelc( tot_lc_arr, output_dir + 'tot_lcfile.dat' )

    cont_fname = output_dir + 'cont_lcfile.dat'
    total_fnames = output_dir + 'tot_lcfile.dat'    
    
    con_dat = get_data(cont_fname, names=[total_names[0]], set_subtractmean=subtract_mean)
    tot_dat = get_data(total_fnames, names=total_names, set_subtractmean=subtract_mean)
    
    if fixed is not None:
        if (fixed[0] == 0) & (fixed[1] == 0):
            skip_init = True
        else:
            skip_init = False
    else:
        skip_init = False    
            

    if not skip_init:
        
        #Run continuum fit to get priors on DRW values
        fnames = get_javelin_filenames(output_chains, output_burn, output_logp, 'cont', output_dir)
        
        if fixed is not None:
            fixed_cont = fixed[:2]
            p_fix_cont = p_fix[:2]
        else:
            fixed_cont = None
            p_fix_cont = None
        
        cmod = Cont_Model(con_dat)
        cmod.do_mcmc(fixed=fixed_cont, p_fix=p_fix_cont, 
                    nwalkers=nwalkers, nburn=nburn, nchain=nchain, threads=threads,
                    fchain=fnames[0], fburn=fnames[1], flogp=fnames[2], 
                    set_verbose=verbose)
        
        #Get HPD from continuum fit
        cmod.get_hpd(set_verbose=False)
        conthpd = cmod.hpd
        
        if fixed is not None:
            if fixed[0] == 0:
                conthpd[0,0] = conthpd[1,0] - 1e-6
                conthpd[2,0] = conthpd[1,0] + 1e-6
                
            if fixed[1] == 0:
                conthpd[0,1] = conthpd[1,1] - 1e-6
                conthpd[2,1] = conthpd[1,1] + 1e-6
                
        #Get histogram figure from continuum fit
        if plot:
            cmod.show_hist(bins=nbin)
    else:
        conthpd = None
        cmod = None
    
    
    
        
    
    
    #Run fit on continuum + line(s)
    fnames = get_javelin_filenames(output_chains, output_burn, output_logp, 'rmap', output_dir)
    
    if rm_type == 'spec':
        rmod = Rmap_Model(tot_dat)
    elif rm_type == 'phot':
        rmod = Pmap_Model(tot_dat)
        
    if len(total_fnames) == 2:
        laglimit = [laglimit]
    
    rmod.do_mcmc(conthpd=conthpd, fixed=fixed, p_fix=p_fix, lagtobaseline=lagtobaseline, laglimit=laglimit,
                 nwalkers=nwalkers, nburn=nburn, nchain=nchain, 
                 threads=threads, fchain=fnames[0], fburn=fnames[1], flogp=fnames[2], 
                 set_verbose=verbose)
    
    #Get HPD from continuum + line(s) fit
    rmod.get_hpd(set_verbose=False)
    rmap_hpd = rmod.hpd
        
    #plot predicted light curves from DRW fits
    med_params = rmap_hpd[1,:]
    rmap_bestfit = rmod.do_pred(med_params)
    
        
    tau = np.exp( rmod.flatchain[:,1] )
    sigma = np.exp( rmod.flatchain[:,0] )
    tophat_params = rmod.flatchain[:, 2:].T    



    return {
        'cont_hpd': conthpd,
        'tau': tau,
        'sigma': sigma,
        'tophat_params': tophat_params,
        'hpd': rmap_hpd,
        'cont_model': cmod,
        'rmap_model': rmod,
        'cont_dat': con_dat,
        'tot_dat': tot_dat
    }        
    
    
from scipy.stats import binned_statistic   
def javelin_pred_lc(rmod, t_cont, t_lines, nbin=None, metric='med'):
    
    """Predict light curve(s) from a JAVELIN RM fit object.


    Parameters
    ----------

    rmod : javelin.lcmodel.Rmap_Model or javelin.lcmodel.Pmap_Model
        The output RM model from ``run_javelin``, after fitting the data.
        
    t_cont : array_like
        The time array for the continuum light curve.
        
    t_lines : array_like
        The time array for the line light curve(s).

    metric : str, optional
        The metric to use to get the bestfit parameter from the parameter distributions.
        Can be either 'mean' or 'med'. Default is 'med'.


    Returns
    -------

    rmap_bestfit : javelin.lcmodel.Rmap_Model or javelin.lcmodel.Pmap_Model
        The RM model using the best fit parameters to predict the light curve(s).
    
    """
    
    tau = rmod.flatchain[:,1]
    sigma = rmod.flatchain[:,0]
    tophat_params = rmod.flatchain[:, 2:].T
    
    Nlc = len(tophat_params)//3
    
    if metric == 'med':
        func = np.median
    if metric == 'mean':
        func = np.mean
            
    
    
    tau_best = func(tau)
    sigma_best = func(sigma)
    
    tophat_best = np.zeros( Nlc*3 )
    
    for i in range(Nlc):
        for j in range(3):
            
            
            if j == 0:
                mc_vals = tophat_params[3*i + j, :]               
                tophat_best[ 3*i + j ] = func(mc_vals)
                
            else:
                tophat_best[ 3*i + j ] = func(tophat_params[ 3*i + j, : ])


    bestfit_vals = np.concatenate([ [sigma_best], [tau_best], tophat_best ])
    rmap_bestfit = rmod.do_pred(bestfit_vals)

    return rmap_bestfit
