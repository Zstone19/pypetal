import os
import subprocess
import time
import threading
import sys
import signal

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from pyzdcf import pyzdcf

import pypetal.drw_funcs as drw
import pypetal.pyccf as pyccf

import pypetal.pyroa_funcs as pyroa
import PyROA

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

    times : list of astropy.units.Quantity
        Times for the light curve.

    data : list of astropy.units.Quantity
        Data for the light curve.

    error : list of astropy.units.Quantity
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

    t = np.linspace( times[0].value - extra_t, times[-1].value + extra_t, 1000 ).tolist()
    for i in range(len(times)):
        t.append(times[i].value)

    sort_ind = np.argsort(t)
    t = np.array(t)[sort_ind]

    mu, var = gp.predict(data.value, t, return_var=True)
    std = np.sqrt(var)

    mu_flag = []
    for i in range(len(times)):
        ind = np.argwhere(t == times[i].value).T[0][0]
        mu_flag.append(mu[ind])

    #Reject if data is beyond nsig*sig of fit mean
    flag_mask = np.abs(data.value - mu_flag) > nsig*error.value

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

    #Make sure plike dir exists
    assert os.path.exists( plike_dir )
    plike_dir = os.path.abspath(plike_dir) + r'/'

    cwd = os.getcwd()
    os.chdir(plike_dir)

    #Delete old plike.out if it exists
    if os.path.exists( plike_dir + 'plike.out' ):
        os.remove( plike_dir + 'plike.out' )

    #Make sure dcf file exists
    assert os.path.exists(dcf_fname)
    dcf_fname = os.path.abspath(dcf_fname)


    #Make sure there are two lag bounds
    if len(lag_bounds) != 2:
        print('Lag bounds: ', lag_bounds)
        raise ValueError('Must provide two lag bounds.')


    #Make file with the arguments to pass to plike
    with open('args.txt', 'w') as f:
        f.write(dcf_fname + '\n')
        f.write(str(lag_bounds[0]) + '\n')
        f.write(str(lag_bounds[1]))


    if verbose:
        print('Executing PLIKE')

    exec_str = './plike < args.txt'
    res = subprocess.Popen(exec_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output, error = res.communicate()
    if res.returncode != 0:
        raise Exception("File handling failed %d %s %s" % (res.returncode, output, error))


    #Delete args.txt
    os.remove('args.txt')


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
            Path to the directory containing the light curves.

    fname1 : str
            Name of the first light curve file.

    fname2 : str
            Name of the second light curve file.

    out_dir : str
            Path to the directory to place the output ZDCF file.

    num_MC : float, optional
            The number of Monte Carlo simulations to run. Default is 500.

    minpts : int, optional
            The minimum number of points to use in each bin when computing the ZDCF.
            Must be larger than 11. If set to 0, it will be set to 11. Default is 0.

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
            if there are more than 3000 data points per light curve. Default is 'auto'.

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
        Path to the first light curve file .

    fname2 : str
        Path to the second light curve file .

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


    .. note:: Both light curve files must be in CSV format with the following columns in order: time, value, uncertainty.



    Returns
    -------

    result: dict
        Dict of output results, containing:

        * 'CCF' : list of floats
            The ICCF.

        * 'CCF_lags' : list of floats
            The lags corresponding to the ICCF.

        * 'centroid' : float
            The median of the CCCD.

        * 'centroid_err_hi' : float
            The upper error of the centroid.

        * 'centroid_err_lo' : float
            The lower error of the centroid.

        * 'peak' : float
            The median of the CCPD.

        * 'peak_err_hi' : float
            The upper error of the peak.

        * 'peak_err_lo' : float
            The lower error of the peak.

        * 'CCCD_lags' : list of floats
            The CCCD.

        * 'CCPD_lags' : list of floats
            The CCPD.

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


##############################################################
######################### PyROA ##############################
##############################################################

def close_pool():
    global pool
    pool.close()
    pool.terminate()
    pool.join()

def term(*args,**kwargs):
    sys.stderr.write('\nStopping...')
    # httpd.shutdown()
    stoppool=threading.Thread(target=close_pool)
    stoppool.daemon=True
    stoppool.start()


def run_pyroa(fnames, lc_dir, line_dir, line_names,
              nburn=10000, nchain=15000, lag_bounds=None, 
              init_tau=None, init_delta=10, sig_level=100,
              together=True, subtract_mean=True, div_mean=False,
              add_var=False, delay_dist=False, psi_types='Gaussian',
              objname=None):
    
    
    """Run PyROA for a number of input light curves.
    NOTE: This will assume a log-Gaussian for the delay distribution.
    
    
    Parameters
    ----------
    
    fnames : list of str
        A list of paths to the light curves to be used in PyROA. The first light curve will be assumed to be the continuum. Must be in CSV format.
        
    lc_dir : str
        The output directory to put the light curves for PyROA to use.
        
    line_dir : str, list of str
        The output directory to store the data products output from PyROA. 
        If together=False, this should be a list of paths for each line. If together=True, this should be a single path.
        
    line_names : list of str
        A list of the line names corresponding to the light curves.
        
    nburn : int, optional
        The number of burn-in samples to discard. Default is 10000.
        
    nchain : int, optional
        The number of samples to get per walker. Default is 15000.
        
    lag_bounds : list of str, optional
        The range of lags to consider for PyROA, one set of bounds for each line (excluding the continuum). If ``None``, the baseline of
        each light curve will be used as the lag bounds. Default is ``None``.
        
    init_tau : list of float
        A list of inital values for the time delay (tau) for each light curve. If ``None``, will set the inital value to 10.0 for each line.
        Default is ``None``.

    sig_level : float
        The number to use for sigma clipping in PyROA. Default is 100.
        
    together : bool, optional
        Whether to fit all light curves together or not. If ``together=False``, the ``line_dirs`` argument must be set. Default is ``True``.
        
    subtract_mean : bool, optional
        Whether to subtract the mean from all light curves before using PyROA or not. Will occur after ``div_mean`` if set to ``True``. Default is ``True``.
        
    div_mean : bool, optional
        Whether to divide each light curve by its mean before using PyROA. Will occur before ``subtract_mean`` if set to ``True``. Default is ``False``.
        
    add_var : bool or list of bool, optional
        Whether or not to add additional uncertainty in the data, same as the PyROA argument. If ``together=False``, multiple values may be given for each line. 
        If only one value is given, it will be assumed for all lines. Default is ``True``.
        
    delay_dist : bool or list of bool, optional
        Same as the ``delay_dist`` argument for PyROA. If ``together=False``, multiple values may be given for each line. 
        If only one value is given, it will be assumed for all lines. Default is ``True``.
        
    objname : str, optional
        The name of the object, will be used for plot and for the saved PyROA light curve data. If ``None``, will be set to "pyroa".
        
    line_names : list of str
        A list of directories to place the output PyROA data for each of the lines (excluding the continuum). Must be set if ``together=False``, and will only be used in such a case.
        Default is ``None``.
        
        
    Returns
    -------
    
    fit : PyROA.Fit or list of pyROA.Fit
        The PyROA.Fit object output from PyROA. If ``together=False``, this will be an array of the ``Fit`` objects.

    
    """

    
    if objname is None:
        objname = 'pyroa'
    
    if init_tau is None:
        init_tau = np.full( len(fnames), 10. )
        
    if isinstance(psi_types, str):
        psi_types = [psi_types] * ( len(fnames)-1 )
        
    #If no lag bounds given, use the baseline of the light curves
    if lag_bounds is None:
        lag_bounds = []
        
        x_cont, _, _ = np.loadtxt( fnames[0], unpack=True, usecols=[0,1,2] )
        for i in range(1, len(fnames)):
            x, _, _ = np.loadtxt( fnames[i], unpack=True, usecols=[0,1,2] )
            
            max_x = np.max([ np.max(x), np.max(x_cont) ])
            min_x = np.min([ np.min(x), np.min(x_cont) ])
            bl = max_x - min_x
            
            lag_bounds.append([-bl, bl])
            
            
    if (not together) & (len(fnames)==2) & isinstance(line_dir, str):
        line_dir = [line_dir]
        
            
    if isinstance(line_dir, str):
        os.makedirs(line_dir, exist_ok=True)
    elif isinstance(line_dir, list):
        for ldir in line_dir:
            os.makedirs(ldir, exist_ok=True)
    else:
        raise ValueError('line_dir must be a string or list of strings')
    
    
    os.makedirs(lc_dir, exist_ok=True)
    
    prior_arr = pyroa.get_priors(fnames, lag_bounds, subtract_mean=subtract_mean, div_mean=div_mean, together=together, delimiter=',')
    _ = pyroa.save_lines(fnames, line_names, lc_dir, objname=objname, subtract_mean=subtract_mean, div_mean=div_mean, delimiter=',')

    cwd = os.getcwd()

    if not together:
        
        assert isinstance(line_dir, list), 'Must provide multiple line_dir if together=False'
        
        if isinstance(add_var, bool):
            add_var = np.full( len(fnames)-1, add_var )
            
        if isinstance(delay_dist, bool):
            delay_dist = np.full( len(fnames)-1, delay_dist )
        
        
        fit_arr = []
        
        for i in range(len(fnames)-1):
            
            filters = [line_names[0], line_names[i+1]]
            fit = PyROA.Fit(lc_dir, objname, filters, prior_arr[i,:,:], add_var=add_var[i],
                            init_tau=[init_tau[i]], init_delta=init_delta, sig_level=sig_level,
                            delay_dist=delay_dist[i], psi_types=[psi_types[i]],
                            Nsamples=nchain, Nburnin=nburn)

            pyroa.move_output_files(cwd, line_dir[i])

            fit_arr.append(fit)
            
            signal.signal(signal.SIGTERM, term)
            signal.signal(signal.SIGINT, term)
            signal.signal(signal.SIGQUIT, term)
 
        return fit_arr
        
    else:    
        
        assert isinstance(line_dir, str), 'Must provide one line_dir if together=True'
                
        fit = PyROA.Fit(lc_dir, objname, line_names, prior_arr, add_var=add_var,
                    init_tau=init_tau, init_delta=init_delta, sig_level=sig_level,
                    delay_dist=delay_dist, psi_types=psi_types,
                    Nsamples=nchain, Nburnin=nburn)
        
        pyroa.move_output_files(cwd, line_dir)
    
        return fit
