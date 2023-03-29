import warnings

import celerite
import emcee
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stat
from astropy.timeseries import LombScargle
from astropy.visualization import quantity_support
from celerite import terms
from scipy.optimize import curve_fit, differential_evolution, minimize

from pypetal.drw_rej.utils import MCMC_fit
from pypetal.drw_rej.plotting import plot_outcome

quantity_support()

##############################################################
###################### ASSIST FUNCTIONS ######################
##############################################################

def celerite_fit(x, y, yerr, kernel, nwalkers, nburn, nsamp,
                 solver='minimize', suppress_warnings=True, jitter=True):

    """Fit time-series data to a given Gaussian process kernel using celerite.


    Parameters
    ----------

    x : list of float
        Time values of the data.

    y : list of float
        Values of the data.

    yerr : list of float
        Uncertainties in the data.

    kernel : celerite.terms.Term
        Kernel to fit to the data.

    nwalkers : int
        Number of walkers to use in the MCMC.

    nburn : int
        Number of burn-in steps to use in the MCMC.

    nsamp : int
        Number of samples to use in the MCMC.

    solver : str, optional
        Solver to use for the GP fit. Options are "minimize" for ``scipy.optimize.minimize`` and
        "diff_evo" for ``scipy.optimize.differential_evolution``. Default is "minimize".

    suppress_warnings : bool, optional
        Suppress warnings from the GP fit. Default is ``True``.

    jitter : bool, optional
        If true, fit for a noise (jitter) term in the GP. Default is ``True``.



    Returns
    -------

    samples : array_like
        Samples of the kernel parameters from the MCMC fit.

    gp : celerite.GP
        GP object used for the fit.

    statuses : list of bool
        Status of the fit to the light curve. There are three statuses given:

        - baseline_good
            If ``True``, tau > baseline/10.

        - cadence_good
            If ``True``, tau > mean cadence of the light curve.

        - stn_good
            If ``True``, the DRW sigma parameter is greater than the noise in the light curve.

    """


    if suppress_warnings:
        warnings.filterwarnings("ignore")

    #Set up gp for celerite
    gp = celerite.GP(kernel, np.mean(y.value), fit_mean=False)
    gp.compute(x.value, yerr.value)


    # Define a cost function
    def neg_ll(params, y, gp):
        gp.set_parameter_vector(params)
        return -gp.log_likelihood(y)

    def grad_neg_ll(params, y, gp):
        gp.set_parameter_vector(params)
        return -gp.grad_log_likelihood(y)[1]

    # Define the log probablity
    def log_probability(params):
        gp.set_parameter_vector(params)
        lp = gp.log_prior()
        if not np.isfinite(lp):
            return -np.inf
        return gp.log_likelihood(y) + lp


    #Bounds and values for input parameters
    init_params = gp.get_parameter_vector()
    bounds = gp.get_parameter_bounds()

    #Get soln from extrema of maximum likelihood
    if solver == 'minimize':
        soln = minimize(neg_ll, init_params, jac=grad_neg_ll,
                        method="L-BFGS-B", bounds=bounds, args=(y.value, gp))
        initial = np.array(soln.x)
    elif solver == 'diff_evo':
        soln = differential_evolution(neg_ll, bounds=bounds, args=(y.value, gp))
        initial = np.array(soln.x)


    #Set parameter vector to fit params
    gp.set_parameter_vector(initial)

    #Set up MCMC
    ndim = len(initial)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability)

    #Run burn-in on the MCMC
    p0 = initial + 1e-8 * np.random.randn(nwalkers, ndim)
    p0, lp, _ = sampler.run_mcmc(p0, nburn)

    #Run the MCMC
    sampler.reset()
    sampler.run_mcmc(p0, nsamp)


    #Get samples from the posterior dist
    samples = sampler.flatchain

    #Get median values for the samples from the posterior
    s = np.median(samples, axis=0)
    gp.set_parameter_vector(s)

    tau = 1/np.exp(s[1])
    sig_drw = np.sqrt( np.exp(s[0])/2 )

    if jitter is True:
        sig_n = np.exp(s[2])
    else:
        sig_n = 0

    baseline = x[-1] - x[0]
    if tau < .1*baseline.value:
        baseline_good = True
    else:
        baseline_good = False


    cadence = np.mean(np.diff(x.value))
    if tau > cadence:
        cadence_good = True
    else:
        cadence_good = False

    sig_drw = np.exp( s[1] )
    dy = np.mean(yerr.value)

    if sig_drw**2 > sig_n**2 + dy**2:
        stn_good = True
    else:
        stn_good = False

    return samples, gp, [baseline_good, cadence_good, stn_good]



def MCMC_fit(x, y, yerr, nwalkers=32, nburn=300, nsamp=1000,
             solver='minimize', suppress_warnings=True, jitter=True,
             clip=True):

    """Fit time-series data to a DRW using celerite and emcee.


    Parameters
    ----------

    x : list of astropy.units.Quantity
        Time values of the data.

    y : list of astropy.units.Quantity
        Values of the data.

    yerr : list of astropy.units.Quantity
        Uncertainties on the data.

    nwalkers : int, optional
        Number of walkers to use in the MCMC. Default is 32.

    nburn : int, optional
        Number of burn-in steps to use in the MCMC. Default is 300.

    nsamp : int, optional
        Number of samples to use in the MCMC. Default is 1000.

    solver : str, optional
        Solver to use for the GP fit. Options are "minimize" for ``scipy.optimize.minimize`` and
        "diff_evo" for ``scipy.optimize.differential_evolution``. Default is "minimize".

    suppress_warnings : bool, optional
        Suppress warnings from the GP fit. Default is ``True``.

    jitter : bool, optional
        If true, fit for a noise (jitter) term in the GP. Default is ``True``.

    clip : bool, optional
        If ``True``, clip data points too close together in time (< 1e-8 days). Default is ``True``.



    Returns
    -------

    samples : list of float
        Samples of the kernel parameters from the MCMC fit

    gp : celerite.GP
        GP object used for the fit

    statuses : list of bool
        Status of the fit to the light curve. There are three statuses given:

        - baseline_good
            If True, tau > baseline/10

        - cadence_good
            If true, tau > mean cadence of the light curve

        - stn_good
            If true, the DRW sigma parameter is greater than the noise in the light curve

    """


    #Set up GP in celerite
    baseline = x[-1] - x[0]


    #Bounds and value for a
    min_prec = np.min(yerr.value)
    amp = np.max( y.value + yerr.value ) - np.min(y.value - yerr.value)

    amin = np.log(.001*min_prec)
    amax = np.log(10*amp)

    aval = np.mean([amin, amax])


    #Bounds and value for c
    if clip is True:
        min_cadence = np.clip(np.min( np.diff(x.value) ), 1e-8, None)
    else:
        min_cadence = np.min( np.diff(x.value) )

    cmin = np.log( 1/10/baseline.value )
    cmax = np.log(1/min_cadence)

    cval = np.mean([cmin, cmax])

    bounds = dict(log_a=(amin, amax), log_c=(cmin, cmax))
    kernel = terms.RealTerm(log_a=aval, log_c=cval, bounds=bounds )


    if jitter is True:
        #Bounds and value for s (jitter)
        smin = -10
        smax = np.log(amp)
        sval = np.mean([smin, smax])

        bounds = dict(log_sigma=(smin, smax))
        kernel += terms.JitterTerm(log_sigma=sval, bounds=bounds)

    samples, gp, statuses = celerite_fit(x, y, yerr, kernel, nwalkers,
                                         nburn, nsamp, solver,
                                         suppress_warnings, jitter)

    return samples, gp, statuses



def psd_from_gp(fLS, powerLS, gp, samples, baseline):

    """Calculate the PSD of a light curve using the DRW fit from celerite.
    This assumes that a Lomb-Scargle periodogram has been made from the
    light curve, giving a list of frequencies and powers.


    Parameters
    ----------

    fLS : list of astropy.units.Quantity
        Frequencies from the Lomb-Scargle periodogram.

    powerLS : list of astropy.units.Quantity
        Powers from the Lomb-Scargle periodogram.

    gp : celerite.GP
        GP object used for the fit.

    samples : array_like
        Samples of the kernel parameters from the MCMC fit.

    baseline : astropy.units.Quantity
        Baseline of the light curve.



    Returns
    -------

    f_eval : list of astropy.units.Quantity
        Frequencies at which the PSD is evaluated

    psd_credint : list of astropy.units.Quantity
        An (n,3) array of the PSD, with the first columns being
        the 16th percentile, the second column being the median,
        and the third column being the 84th percentile.

    """


    #Get freqs
    f_eval = np.logspace( np.log10(fLS[0].value), np.log10(fLS[-1].value), 200 )

    psd_samples = np.empty( (len(f_eval), len(samples)) )
    kernel = gp.kernel

    #Get the PSD for each parameter-set in the sample
    for i, s in enumerate(samples):
        gp.set_parameter_vector(s)
        psd_samples[:, i] = kernel.get_psd(2*np.pi*f_eval)/(2*np.pi)


    #Get credibility intervals (16, 50, 84)
    psd_credint = np.empty((len(f_eval), 3))
    psd_credint[:, 0] = np.percentile(psd_samples, 16, axis=1)
    psd_credint[:, 2] = np.percentile(psd_samples, 84, axis=1)
    psd_credint[:, 1] = np.median(psd_samples, axis=1)


    #Normalize the PSD
    f_norm = np.max(powerLS.value[fLS.value>1/(2*np.pi*0.2*baseline.value)])/psd_credint[0, 1]
    #Max[ power[ f > 1/2pi/baseline/5 ] ]

    psd_credint[:, 0] = psd_credint[:, 0]*f_norm
    psd_credint[:, 2] = psd_credint[:, 2]*f_norm
    psd_credint[:, 1] = psd_credint[:, 1]*f_norm

    return f_eval, psd_credint





def binLS(fLS, powerLS_samp, num_bins):

    """Bin the Lomb-Scargle periodogram by frequency.

    Parameters
    ----------
    fLS : list of astropy.units.Quantity
        Frequencies from the Lomb-Scargle periodogram.

    powerLS_samp : list of astropy.units.Quantity
        Samples of the Lomb-Scargle periodogram.

    num_bins : int
        Number of bins to use for the binned periodogram.



    Returns
    -------

    binCenters : list of astropy.units.Quantity
        Bin centers for the binned periodogram.

    bin_vals : list of astropy.units.Quantity
        The values of the binned periodogram.

    bin_errs : list of astropy.units.Quantity
        The upper error on the binned periodogram.

    lower_err : list of astropy.units.Quantity
        The lower error on the binned periodogram.

    """


    #Get the credibility intervals for each point on the PSD
    bins_credint = np.empty((len(fLS), 3))
    bins_credint[:, 0] = np.percentile(powerLS_samp, 16, axis=0)
    bins_credint[:, 2] = np.percentile(powerLS_samp, 84, axis=0)
    bins_credint[:, 1] = np.median(powerLS_samp, axis=0)

    binEdges = np.logspace( np.log10(fLS[0].value), np.log10(fLS[-1].value), num_bins+1)

    #Get data for bins (16th, 50th, 84th percentile for each bin mean)
    draws_binned = np.empty((len(binEdges)-1, 3))
    draws_binned[:,0], _, _ = stat.binned_statistic(fLS, bins_credint[:, 0],
                                                                     statistic='mean',
                                                                     bins=binEdges)
    draws_binned[:,2], _, _ = stat.binned_statistic(fLS, bins_credint[:, 2],
                                                                     statistic='mean',
                                                                     bins=binEdges)
    draws_binned[:,1], _, _ = stat.binned_statistic(fLS, bins_credint[:, 1],
                                                                     statistic='mean',
                                                                     bins=binEdges)

    #Get error (upper and lower) for each bin value
    bin_err = np.array([draws_binned[:, 1] - draws_binned[:, 0],
                        draws_binned[:, 2] - draws_binned[:, 1]])

    binEdges[0] *= .99


    binCenters = []

    for i in range(len(binEdges) - 1):

        #Get center of the bin
        cent = (np.log10(binEdges[i]) + np.log10(binEdges[i+1]) )/2

        binCenters.append( 10**cent )

    bin_vals = draws_binned[:, 1]
    lower_err = bin_err[1]
    bin_err = bin_err[0]

    binCenters = np.concatenate( ([fLS[0].value], binCenters ) )
    bin_vals = np.concatenate( ([bin_vals[0]], bin_vals ) )
    bin_err = np.concatenate( ([0.], bin_err ) )
    lower_err = np.concatenate( ([0.], lower_err ) )

    binCenters = np.append( binCenters, fLS[-1].value )
    bin_vals = np.append( bin_vals, bin_vals[-1] )
    bin_err = np.append(bin_err, 0.)
    lower_err = np.append(lower_err, 0.)


    #Mask out infinite bins
    mask_finite = np.isfinite(bin_vals)
    binCenters = np.array(binCenters)[mask_finite]
    bin_vals = np.array(bin_vals)[mask_finite]
    bin_err = np.array(bin_err)[mask_finite]
    lower_err = np.array(lower_err)[mask_finite]

    return binCenters, bin_vals, bin_err, lower_err




def smoothly_broken_power_law(f, A=1, f_br=1e-3, alpha=0, beta=2):

    """A smoothly broken power law:

    .. math:: P(f) = \frac{A}{(f / f_{br})^\alpha + (f/f_{br})^\beta}



    Parameters
    ----------

    f : list of float
        Frequencies at which to evaluate the power law.

    A : float
        The amplitude of the power law.

    f_br : float
        The break frequency of the power law.

    alpha : float
        The slope of the power law below the break frequency.

    beta : float
        The slope of the power law above the break frequency.



    Returns
    -------

    power : list of float
        The power law.

    """

    return A/((f/f_br)**alpha + (f/f_br)**beta)



def psd_sbpl(f, psd, err, p0, bounds):

    """Generate a smoothly broken power law fit to the PSD.


    Parameters
    ----------

    f : list of float
        Frequencies for the input PSD.

    psd : list of float
        The input PSD.

    err : list of float
        The error on the input PSD.

    p0 : list of float
        The initial guess for the parameters of the smoothly broken power law.

    bounds : list of float
        The bounds on the parameters of the smoothly broken power law.




    Returns
    -------

    fit_vals : list of float
        The best fit parameters for the smoothly broken power law.

    fit_errs : list of float
        The errors on the best fit parameters for the smoothly broken power law.

    """


    try:
        soln = curve_fit(smoothly_broken_power_law, f, psd,
                         sigma=err, p0=p0, bounds=bounds, maxfev=10000)
        fit_vals = soln[0]
        fit_err = np.sqrt(np.diag(soln[1]))
    except:
        soln = curve_fit(smoothly_broken_power_law, f, psd,
                         sigma=None, p0=None, bounds=bounds, maxfev=10000)
        fit_vals = soln[0]
        fit_err = np.sqrt(np.diag(soln[1]))

    return fit_vals, fit_err


def psd_data(x, y, yerr, samples, gp, nsamp=20):


    """Generate all PSD data for the summary plot.


    Parameters
    ----------

    x : list of astropy.units.Quantity
        The light curve times.

    y : list of astropy.units.Quantity
        The light curve values.

    yerr : list of astropy.units.Quantity
        The uncertainty in the light curve.

    samples : list of float
        The samples from the MCMC fit.

    gp : celerite.GP
        The GP object used to fit the light curve.

    nsamp : int
        The number of samples of the Lomb-Scargle periodogram to use.



    Returns
    -------

    fLS : list of astropy.unit.Quantity
        The frequencies for the Lomb-Scargle periodogram.

    powerLS : list of astropy.unit.Quantity
        The power for the Lomb-Scargle periodogram.

    f_eval : list of astropy.unit.Quantity
        The frequencies for the PSD from the celerite fit.

    psd_credint : list of astropy.unit.Quantity
        The (16th, 50th, 84th) percentiles of the PSD from the celerite fit.

    bin_vals : list of astropy.unit.Quantity
        The values of the Lomb-Scargle periodogram.

    bin_err : list of astropy.unit.Quantity
        The upper error on the binned Lomb-Scargle periodogram.

    lower_err : list of astropy.unit.Quantity
        The lower error on the binned Lomb-Scargle periodogram.

    binCenters : list of astropy.unit.Quantity
        The centers of the bins for the binned Lomb-Scargle periodogram.

    fit_vals : list of astropy.unit.Quantity
        The best fit parameters for the smoothly broken power law fit to
        the Lomb-Scargle periodogram.

    fit_errs : list of astropy.unit.Quantity
        The errors on the best fit parameters for the smoothly broken
        power law fit to the Lomb-Scargle periodogram.

    """

    baseline = x[-1] - x[0]
    powerLS_samps = []


    for _ in range(nsamp):
        y_samp = y.value + np.random.normal(0, yerr.value)
        y_samp *= y.unit

        #Get Lomb-Scargle
        fLS, powerLS = LombScargle(x, y_samp, yerr).autopower(normalization='psd')

        #Put power into lightkurve units (flux variance / frequency unit)
        fs = (1./(np.min(np.diff(x.value)[np.diff(x)>0])))
        powerLS *= 2/(len(x)*fs)


        ind = np.argsort(fLS.value)
        fLS = fLS[ind]
        powerLS = powerLS[ind]

        powerLS_samps.append(powerLS)

    f_eval, psd_credint = psd_from_gp(fLS, powerLS, gp, samples, baseline)

    binCenters, bin_vals, bin_err, lower_err = binLS(fLS, powerLS_samps, 15)


    s = np.median(samples, axis=0)
    tau = 1/np.exp(s[1])
    f = 1/2/np.pi/tau

    bounds = np.array([(.1*np.nanmax(bin_vals), 10*np.nanmax(bin_vals)),
                       (fLS[0].value, fLS[-1].value),
                       (-3., 0.), (.1, 5.)]).T
    p0 = [np.nanmax(bin_vals), f, 0., 2.]

    fit_vals, fit_err = psd_sbpl(binCenters[1:-4], bin_vals[1:-4],
                                 None, p0, bounds)

    return fLS, powerLS, f_eval, psd_credint, \
           bin_vals, bin_err, binCenters, lower_err, \
           fit_vals, fit_err


##############################################################
##############################################################
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


    if type(data[0]) == u.Quantity:
        data_unit = data[0].unit

    flag_mask = np.zeros( times.shape, dtype=bool )


    #Assume 'times' and 'fluxes' have the same shape

    #Fit to the DRW model
    samples, gp, statuses = MCMC_fit(times, data, error,
                                         nwalkers=nwalkers, nburn=nburn, nsamp=nsamp,
                                         jitter=jitter, clip=clip)


    fig, ax = plot_outcome(times, data, error, samples, gp, data_unit,
                                nsig=nsig, target=target, show_mean=True,
                                filename=fname, jitter=jitter, show=plot)

    plt.cla()
    plt.clf()
    plt.close()

    tau_vals = 1/np.exp(samples[:, 1])
    sig_vals = np.sqrt( np.exp(samples[:, 0])/2 )

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
