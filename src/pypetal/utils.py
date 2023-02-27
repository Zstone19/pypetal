import glob
import os
import shlex
import subprocess
import time

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from javelin.lcio import writelc
from javelin.lcmodel import Cont_Model, Pmap_Model, Rmap_Model
from javelin.zylc import get_data
from pyzdcf import pyzdcf

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
