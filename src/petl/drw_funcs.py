import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.gridspec as gridspec
from astropy import units as u

import corner
import celerite
import celerite.terms as terms
import emcee

from scipy.optimize import minimize, curve_fit, differential_evolution
import scipy.stats as stat
from astropy.timeseries import LombScargle

import os

from astropy.visualization import quantity_support
quantity_support()

import warnings

def str2float_col(array, dtype):
    
    if type(array) == type(str):
        array = array[1:-1].split()  
        array = np.array(array, dtype=dtype)
        return array
    else:
        return array

def hampel_filter(x, y, window_size, n_sigmas=3):
    """
    Perform outlier rejection using a Hampel filter
    
    x: time (list or np array)
    y: value (list or np array)
    window_size: window size to use for Hampel filter
    n_sigmas: number of sigmas to reject outliers past
    
    returns: x, y, mask [lists of cleaned data and outlier mask]
        
    Adapted from Eryk Lewinson
    https://towardsdatascience.com/outlier-detection-with-hampel-filter-85ddf523c73d
    """
    
    # Ensure data are sorted
    if np.all(np.diff(x) > 0):
        ValueError('Data are not sorted!')
        
    x0 = x[0]
    
    n = len(x)
    outlier_mask = np.zeros(n)
    k = 1.4826 # MAD scale factor for Gaussian distribution
    
    # Loop over data points
    for i in range(n):
        # Window mask
        mask = (x > x[i] - window_size) & (x < x[i] + window_size)
        if len(mask) == 0:
            idx.append(i)
            continue
        # Compute median and MAD in window
        y0 = np.median(y[mask])
        S0 = k*np.median(np.abs(y[mask] - y0))
        # MAD rejection
        if (np.abs(y[i] - y0) > n_sigmas*S0):
            outlier_mask[i] = 1
            
    outlier_mask = outlier_mask.astype(np.bool)
    
    return np.array(x)[~outlier_mask], np.array(y)[~outlier_mask], outlier_mask


def celerite_fit(x, y, yerr, kernel, nwalkers, nburn, nsamp, solver='minimize', suppress_warnings=True, jitter=True):
    
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
    
    if jitter == True:
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



def MCMC_fit(x, y, yerr, nwalkers=32, nburn=300, nsamp=1000, solver='minimize', suppress_warnings=True, jitter=True, clip=True):

    #Set up GP in celerite
    baseline = x[-1] - x[0]


    #Bounds and value for a
    min_prec = np.min(yerr.value)
    amp = np.max( y.value + yerr.value ) - np.min(y.value - yerr.value)
    
    amin = np.log(.001*min_prec)
    amax = np.log(10*amp)

    aval = np.mean([amin, amax])


    #Bounds and value for c
    if clip == True:
        min_cadence = np.clip(np.min( np.diff(x.value) ), 1e-8, None)
    else:
        min_cadence = np.min( np.diff(x.value) )

    cmin = np.log( 1/10/baseline.value )
    cmax = np.log(1/min_cadence)

    cval = np.mean([cmin, cmax])

    bounds = dict(log_a=(amin, amax), log_c=(cmin, cmax))
    kernel = terms.RealTerm(log_a=aval, log_c=cval, bounds=bounds )


    if jitter == True:
        #Bounds and value for s (jitter)
        smin = -10
        smax = np.log(amp)
        sval = np.mean([smin, smax])

        bounds = dict(log_sigma=(smin, smax))
        kernel += terms.JitterTerm(log_sigma=sval, bounds=bounds)

    samples, gp, statuses = celerite_fit(x, y, yerr, kernel, nwalkers, nburn, nsamp, solver, suppress_warnings, jitter)

    return samples, gp, statuses







def CARMA_fit(x, y, yerr, p=2, nburn=300, nsamp=1000, solver='minimize', suppress_warnings=True):

    #p = J, q = p-1

    #Assume bounds for kernel
    amin = -10; amax = 10; aval = 0
    bmin = -10; bmax = 10; bval = 0
    cmin = -10; cmax = 10; cval = 0
    dmin = -10; dmax = 10; dval = 0
    
    smin = -10; smax = 10; sval = 0

    bounds = dict( log_a=(amin, amax), log_b=(bmin, bmax), log_c=(cmin, cmax), log_d=(dmin, dmax) )
    kernel = terms.ComplexTerm(log_a=aval, log_b=bval, log_c=cval, log_d=dval, bounds=bounds)

    for j in range(2, p+1):
        kernel += terms.ComplexTerm(log_a=aval, log_b=bval, log_c=cval, log_d=dval, bounds=bounds)

    # Add jitter term
    bounds = dict( log_sigma=(smin, smax) )
    kernel += terms.JitterTerm(log_sigma=sval, bounds=bounds)



    samples, gp, statuses = celerite_fit(x, y, yerr, kernel, nburn, nsamp, solver, suppress_warnings)
    
    return samples, gp, statuses





def psd_from_gp(fLS, powerLS, gp, samples, baseline):

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
    
    #Get the credibility intervals for each point on the PSD
    bins_credint = np.empty((len(fLS), 3))
    bins_credint[:, 0] = np.percentile(powerLS_samp, 16, axis=0)
    bins_credint[:, 2] = np.percentile(powerLS_samp, 84, axis=0)
    bins_credint[:, 1] = np.median(powerLS_samp, axis=0)
    
    binEdges = np.logspace( np.log10(fLS[0].value), np.log10(fLS[-1].value), num_bins+1)
 
    #Get data for bins (16th, 50th, 84th percentile for each bin mean)
    draws_binned = np.empty((len(binEdges)-1, 3))
    draws_binned[:, 0], bin_edges, binnumber = stat.binned_statistic(fLS, bins_credint[:, 0], statistic='mean', bins=binEdges)
    draws_binned[:, 2], bin_edges, binnumber = stat.binned_statistic(fLS, bins_credint[:, 2], statistic='mean', bins=binEdges)
    draws_binned[:, 1], bin_edges, binnumber = stat.binned_statistic(fLS, bins_credint[:, 1], statistic='mean', bins=binEdges)    
    
    #Get error (upper and lower) for each bin value
    bin_err = np.array([draws_binned[:, 1] - draws_binned[:, 0], draws_binned[:, 2] - draws_binned[:, 1]])
    
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
    return A/((f/f_br)**alpha + (f/f_br)**beta)


def psd_sbpl(f, psd, err, p0, bounds):
    
    try:
        soln = curve_fit(smoothly_broken_power_law, f, psd, sigma=err, p0=p0, bounds=bounds, maxfev=10000)
        fit_vals = soln[0]
        fit_err = np.sqrt(np.diag(soln[1]))
    except:
        soln = curve_fit(smoothly_broken_power_law, f, psd, sigma=None, p0=None, bounds=bounds, maxfev=10000)
        fit_vals = soln[0]
        fit_err = np.sqrt(np.diag(soln[1]))
    
    return fit_vals, fit_err


def psd_data(x, y, yerr, samples, gp, nsamp=20, psd_fit='sbpl'):
    
    baseline = x[-1] - x[0]
    powerLS_samps = []
    
    
    for i in range(nsamp):
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
    
    fit_vals, fit_err = psd_sbpl(binCenters[1:-4], bin_vals[1:-4], None, p0, bounds)
    
    return fLS, powerLS, f_eval, psd_credint, bin_vals, bin_err, binCenters, lower_err, fit_vals, fit_err


def SF_from_data(x, y, yerr, group=False):
    dt = []
    dm = []
    dm_err = []
    
    if group == False:
        for i in range(len(x)-1):
            for j in range(len(x[i+1:])):
                dt.append(  np.abs(x[j]-x[i]) )
                dm.append( np.abs(y[j] - y[i]) )
                dm_err.append(  np.sqrt( yerr[j]**2 + yerr[i]**2 ) ) 
            
        return dt, dm, dm_err
    
    baseline = x[-1] - x[0]
    bin_num = 50
    
    bin_width = baseline/bin_num
    
    bin_vals = []
    bin_err = []
    for i in range(bin_num):
        bin_vals.append([])
        bin_err.append([])
    
    #Get bin edges and centers
    bin_edges = [ x[0] + n*x[-1]/bin_num for n in range(50) ]
    bin_centers = []
    for n in range(len(bin_edges)-1):
        bin_centers.append( (bin_edges[n+1] + bin_edges[n])/2)
    
    #For all dt bins, gather all dm**2 and dm**2_err within it 
    for i in range(len(x)-1):
        for j in range(len(x[i+1:])):
            
            for n in range(bin_num-1):
                left = bin_edges[n]
                right = bin_edges[n+1]
                
                if left < np.abs(x[i]-x[j]) <= right:
                    bin_vals[n].append(  (y[j]-y[i])**2 )
                    bin_err[n].append(  2*np.abs(y[j]-y[i])* np.sqrt(yerr[j]**2 + yerr[i]**2) )

    print(bin_vals)
    #Get rms_dm and rms_dm_err for each dt bin
    for i in range(len(bin_vals)):
        
        rms_mag_err = 0
        for j in range(len(bin_err[i])):
            rms_mag_err += 2*bin_vals[i][j]*bin_err[i][j]
        rms_mag_err = np.sqrt(rms_mag_err)
        bin_err[i] = rms_mag_err
        
        rms_mag = np.sqrt( np.sum(bin_vals[i])/len(bin_vals[i]) )
        bin_vals[i] = rms_mag
        
    return bin_centers, bin_vals, bin_err
        
    
    
    
def plot_binned_SF(centers, vals, err, group=False):
    
    fig, ax = plt.subplots()
    
    ax.errorbar(centers, vals, yerr=err, 
                marker=None, drawstyle='steps-mid', color='k', 
                linewidth=1.0, capsize=3)
    
    ax.set_ylabel('RMS mag')
    ax.set_xlabel('$\Delta t$')
    
    return fig, ax
    
    
    


def plot_outcome(x, y, yerr, samples, gp, unit, nsig=0, 
                 carma_samples=None, carma_gp=None, 
                 target=None, show_mean=True, psd_fit='sbpl', 
                 filename=None, jitter=True, show=False):
    #Plot probability dists of output params
    
    #tau = 1/c
    #sig = sqrt(a/2)

    baseline = x[-1] - x[0]
    extra_t = int(baseline.value//10)

    if unit == u.mag:
        unit_label = 'Magnitude'
    else:
        unit_label = 'Flux [' + str(unit) + ']'
        
    if target is None:
        target_label = ''
    else:
        target_label = 'Light Curve of ' + target

    tau_vals = 1/np.exp(samples[:, 1])
    sig_vals = np.sqrt( np.exp(samples[:, 0])/2 )
    
    if jitter == True:
        n=3
        jitter_vals = np.exp(samples[:, 2])
        sample_vals = np.vstack((np.log10(sig_vals), np.log10(tau_vals), np.log10(jitter_vals) )).T
        labels = [r'$\log_{10}\ (\sigma_{\rm DRW}$ /mag$)$', r'$\log_{10}\ (\tau_{\rm DRW, obs}$ /days$)$', r'$\log_{10}\ (\sigma_n$ /mag$)$']
    else:
        n=2
        sample_vals = np.vstack((np.log10(sig_vals), np.log10(tau_vals) )).T
        labels = [r'$\log_{10}\ (\sigma_{\rm DRW}$ /mag$)$', r'$\log_{10}\ (\tau_{\rm DRW, obs}$ /days$)$']
        

    #This is an array of the following: [[tau, sig, jit], [tau, sig, jit], ..., [tau, sig, jit]]

    fig, axs = plt.subplots(n,n, figsize=(5,5))
    
    titles = []
                                                                       
    sig = np.log10(np.median(sig_vals))
    sig_err_l = sig - np.log10( np.percentile(sig_vals,16) ) 
    sig_err_u = np.log10( np.percentile(sig_vals,84) ) - sig
    sig_title = r'$' + '{:.2f}'.format(sig) + '^{' + ' +{:.2f}'.format(sig_err_u) +  '}_{' + '-{:.2f}'.format(sig_err_l) + '}$'
    titles.append(sig_title)
               
    tau = np.log10(np.median(tau_vals))
    tau_err_l = tau - np.log10( np.percentile(tau_vals,16) ) 
    tau_err_u = np.log10( np.percentile(tau_vals,84) ) - tau
    tau_title = r'$' + '{:.2f}'.format(tau) + '^{' + ' +{:.2f}'.format(tau_err_u) +  '}_{' + '-{:.2f}'.format(tau_err_l) + '}$'
    titles.append(tau_title)
        
    if jitter == True:
        jit = np.log10(np.median(jitter_vals))
        jit_err_l = jit - np.log10( np.percentile(jitter_vals,16) ) 
        jit_err_u = np.log10( np.percentile(jitter_vals,84) ) - jit
        jit_title = r'$' + '{:.2f}'.format(jit) + '^{' + ' +{:.2f}'.format(jit_err_u) +  '}_{' + '-{:.2f}'.format(jit_err_l) + '}$'
        titles.append(jit_title)
                                                                       
    
    #Make corner plot
    fig = corner.corner(sample_vals, labels=labels,
                        quantiles=[.16, .84], use_math_text=True, titles=titles, fig=fig,
                        label_kwargs=dict(fontsize=12), title_kwargs=dict(fontsize=17))

    axs = np.array(fig.axes).reshape((n, n))

    for i in range(n):
        for j in range(n):
            ax = axs[i,j]
            ax.tick_params('both',labelsize=12)
            ax.tick_params(axis='both', which='both', direction='in')
            ax.tick_params(axis='both', which='major', length=6)
            ax.tick_params(axis='both', which='minor', length=3)
            
            
    axs[0,0].title.set_text(sig_title)
    axs[0,0].title.set_fontsize(17)
    
    axs[1,1].title.set_text(tau_title)
    axs[1,1].title.set_fontsize(17)
    if jitter == True:
        axs[2,2].title.set_text(jit_title)
        axs[2,2].title.set_fontsize(17)
        
    #Red out bad tau regions
    axs[1,1].axvspan(np.log10(0.2*baseline.value), np.log10(100*baseline.value), color= "red", zorder=-5, alpha=0.2)

    if jitter == True:
        axs[2,1].axvspan(np.log10(0.2*baseline.value), np.log10(100*baseline.value), color= "red", zorder=-5, alpha=0.2)

    axs[1,0].axhspan(np.log10(0.2*baseline.value), np.log10(100*baseline.value), color= "red", zorder=-5, alpha=0.2)

    gs = gridspec.GridSpec(ncols=4, nrows=4, figure=fig)


    #--------------------------------------------------------------------------
    #Plot data + posterior probability

    ax1 = fig.add_subplot(gs[0, :])

    box = ax1.get_position()
    box.x0 = box.x0 + 0.2
    box.x1 = box.x1 + 1.0
    box.y0 = box.y0 + 0.4
    box.y1 = box.y1 + 0.9
    ax1.set_position(box)


    t = np.linspace(x[0].value - extra_t, x[-1].value+extra_t, 1000)
    mu, var = gp.predict(y.value, t, return_var=True)
    std = np.sqrt(var)
    
    mu_flag, var_flag = gp.predict(y.value, x.value, return_var=True)
    std_flag = np.sqrt(var_flag)

    flag_mask = np.abs(y.value - mu_flag) > nsig*std_flag
    ax1.errorbar( (x.value - x[0].value)[flag_mask], 
                 y.value[flag_mask], yerr.value[flag_mask], 
                 fmt='.', color='DodgerBlue', capsize=1., alpha=.4, ms=7)
    ax1.errorbar( (x.value - x[0].value)[~flag_mask], 
                y.value[~flag_mask], yerr.value[~flag_mask], 
                fmt='.k', capsize=1., alpha=.7, ms=8)

    
    if show_mean == True:
        ax1.plot(t-x[0].value, mu, color='orange', zorder=-1)
    ax1.fill_between(t-x[0].value, mu+std, mu-std, color='orange', alpha=.3)

    ax1.set_xlabel('Time [days]', fontsize=20)
    ax1.set_ylabel(unit_label, fontsize=20)
    ax1.set_title(target_label, fontsize=21)
    
    ax1.tick_params(which='major', length=6)
    ax1.tick_params(which='minor', length=3)
    
    ax1.set_xlim(-extra_t/2, baseline.value + extra_t/2)

    if unit_label == 'Magnitude':
        ax1.invert_yaxis()

    ax1.tick_params('both', labelsize=17)

    #--------------------------------------------------------------------------
    #Plot LS PSD
    
    fLS, powerLS, f_eval, psd_credint, bin_vals, bin_err, binCenters, lower_err, fit_vals, fit_err = psd_data(x, y, yerr, samples, gp)
    
    if type(carma_samples) != type(None):
        _, _, _, _, _, _, _, _, carma_vals, carma_err = psd_data(x, y, yerr, carma_samples, carma_gp)
        
    ax2 = fig.add_subplot(gs[:, -1])
    fig.set_size_inches([6,6])

    # Move the subplot over
    box = ax2.get_position()
    box.x0 = box.x0 + 0.4
    box.x1 = box.x1 + 1.2
    ax2.set_position(box)

    #Plot
    ax2.loglog(fLS.value, powerLS.value, alpha=.3, label='Data-based PSD')

    #--------------------------------------------------------------------------
    #Get PSD from the celerite posterior dist

    #Fill the PSD b/w the 16th and 84th percentiles
    ax2.fill_between(f_eval, psd_credint[:, 2], psd_credint[:, 0], 
                     alpha=0.3, label='Model PSD', color='orange')

    #--------------------------------------------------------------------------
    #Bin Lomb Scargle

    ax2.errorbar(binCenters, bin_vals, yerr=(lower_err, bin_err), 
                 marker=None, drawstyle='steps-mid', color='k', 
                 linewidth=1.0, capsize=3, label='Binned PSD')
    ax2.legend(loc='upper right', fontsize=15)
    
    #--------------------------------------------------------------------------
    #Fit PSD to a smoothly broken power law    
    
    if psd_fit == 'sbpl':
        sbpl_dat = smoothly_broken_power_law(fLS.value, fit_vals[0], fit_vals[1], fit_vals[2], fit_vals[3])

        ax2.plot(fLS.value, sbpl_dat, color='red', label='SBPL Fit')

    
    #--------------------------------------------------------------------------
    
    #Plot values for f_br
    
    f_br = fit_vals[1]
    f_br_err = fit_err[1]
    
    f_br_plus = f_br + f_br_err
    f_br_minus = f_br - f_br_err
    
    if f_br_minus < 0.0 and np.isfinite(f_br_plus) == True:    
        ls_diff = np.log10(f_br_plus) - np.log10(f_br)
        f_br_minus = 10**( np.log10(f_br) - ls_diff )
    
    
    f_br_val = psd_credint[:, 2][0]   
    
    if np.isfinite(f_br_plus) == True:
        #Dat for f_br line
        t = np.logspace(np.log10(f_br_minus), np.log10(f_br_plus), 100)
        dat = np.full(len(t), f_br_val)
        ax2.plot(t, dat, color='red', linewidth=2)
    
        arrow_head = (f_br, 10**(np.log10(f_br_val)+.2) )
        arrow_base = (f_br, 10**(np.log10(f_br_val)+.6) )
        ax2.errorbar([f_br], [arrow_base[1]], yerr=[arrow_base[1] - arrow_head[1]], 
                     uplims=True, color='red', elinewidth=2, capthick=2)
        
        
        
    if type(carma_samples) != type(None):        
        f_br = carma_vals[1]
        f_br_err = carma_err[1]
    
        f_br_plus = f_br + f_br_err
        f_br_minus = f_br - f_br_err
    
        if f_br_minus < 0.0 and np.isfinite(f_br_plus) == True:    
            ls_diff = np.log10(f_br_plus) - np.log10(f_br)
            f_br_minus = 10**( np.log10(f_br) - ls_diff )
    
    
        f_br_val = 10**( np.log10(f_br_val)+.1 )   
    
        if np.isfinite(f_br_plus) == True:
            #Dat for f_br line
            t = np.logspace(np.log10(f_br_minus), np.log10(f_br_plus), 100)
            dat = np.full(len(t), f_br_val)
            ax2.plot(t, dat, color='orange', linewidth=2)
    
            arrow_head = (f_br, 10**(np.log10(f_br_val)+.2) )
            arrow_base = (f_br, 10**(np.log10(f_br_val)+.6) )
            ax2.errorbar([f_br], [arrow_base[1]], yerr=[arrow_base[1] - arrow_head[1]], 
                         uplims=True, color='orange', elinewidth=2, capthick=2)

    #--------------------------------------------------------------------------
    #Red out bad freqs

    ax2.axvspan( np.min(fLS.value) , 1/(2*np.pi*0.2*baseline.value), color='red', alpha=.3)

    fmax = (1./(2*np.pi*np.mean(np.diff(x)[np.diff(x)>0])))
    ax2.axvspan(fmax.value, np.max(fLS.value), color='red', alpha=0.2)


    ax2.set_xlabel(r'Frequency [days$^{-1}$]', fontsize=18)
    ax2.set_ylabel(r'Power [(mag)$^2$ (days)]', fontsize=18)

    ax2.set_ylim( np.min(powerLS.value) )
    ax2.set_xlim(np.min(fLS.value), np.max(fLS.value))

    ax2.set_xscale('log', nonpositive='clip')
    ax2.set_yscale('log', nonpositive='clip')
    
    ax2.tick_params('both', labelsize=15)
    ax2.tick_params('both', which='both', direction='in')
    ax2.tick_params('both', which='major', length=8)
    ax2.tick_params('both', which='minor', length=3)

    fig = plt.gcf()
    
    if filename != None:
        fig.savefig(f'{filename}', bbox_inches='tight')
    
    if show == True:
        plt.show()

    return fig, axs









#-----------------------------------------------------------------------------------------------------------------------------------------
#Plotting functions


def linmix_regression(x, y, xerr, yerr, nchains=10, parallelize=True, verbose=True):

    import linmix.linmix as linmix

    lm = linmix.LinMix(x,y,xerr,yerr, nchains=nchains, parallelize=parallelize)
    lm.run_mcmc(silent = ~verbose)

    output = lm.chain

    alpha = []
    beta = []
    sigsqr = []
    
    for samp in output:
        alpha.append( samp[0] )
        beta.append( samp[1] )
        sigsqr.append( samp[2] )
        
    sig = np.sqrt(sigsqr)
    
    av = np.median(alpha)
    ae_u = np.percentile(alpha, 84) - av
    ae_l = av - np.percentile(alpha, 16)
    
    bv = np.median(beta)
    be_u = np.percentile(beta, 84) - bv
    be_l = bv - np.percentile(beta, 16)

    sv = np.median(sig)
    se_u = np.percentile(sig, 84) - sv
    se_l = sv - np.percentile(sig, 16)
    
    
    print('a: {:.3f} +/- {:.3f}'.format( bv, (be_u+be_l)/2  ))
    print('b: {:.3f} +/- {:.3f}'.format( av, (ae_u+ae_l)/2  ))
    print('sig: {:.3f} +/- {:.3f}'.format( sv, (se_u+se_l)/2 ))
    
    return output

def linmix_ranges(output, xvals):
    sample_yvals = []
    
    for samp in output:
        sample_yvals.append( samp[1]*xvals + samp[0] )

    yvals_upper = np.percentile(sample_yvals, 84, axis=0)
    yvals_med = np.median(sample_yvals, axis=0)
    yvals_lower = np.percentile(sample_yvals, 16, axis=0)
    
    return yvals_upper, yvals_med, yvals_lower




def mbh_l_tau_plot(df, band, out_dir):

    fig, ax = plt.subplots(1,2, figsize=(15,6))

    markers, caps, bars = ax[0].errorbar(df['log_M_BH'], df['log_tau_drw_rest'],\
                                         xerr=df['log_M_BH_e'], \
                                         yerr=[ df['log_tau_drw_rest_err_l'], df['log_tau_drw_rest_err_u'] ], \
                                         fmt='.k')
    [bar.set_alpha(.1) for bar in bars]
    
    ax[0].set_xlabel('$\log_{10} (M_{BH}) $')
    ax[0].set_ylabel(r'$\log_{10} (\tau_{rest})$')

    markers, caps, bars = ax[1].errorbar(df['log_L5100'], df['log_tau_drw_rest'], \
                                         xerr=df['log_L5100_e'], \
                                         yerr=[ df['log_tau_drw_rest_err_l'], df['log_tau_drw_rest_err_u'] ], \
                                         fmt='.k')
    [bar.set_alpha(.1) for bar in bars]
    
    ax[1].set_xlabel('$\log_{10} (L_{5100}) $', fontsize=15)
    ax[1].set_ylabel(r'$\log_{10} (\tau_{rest})$', fontsize=15)
    
    ax[0].tick_params(axis='x', labelsize=12)
    ax[0].tick_params(axis='y', labelsize=12)

    ax[1].tick_params(axis='x', labelsize=12)
    ax[1].tick_params(axis='y', labelsize=12)
    
    print(out_dir + 'Mbh_v_tau_' + band + '.pdf')
    plt.savefig(out_dir + 'Mbh_v_tau_' + band + '.pdf')

    return fig, ax


def tau_z_plot(df, band, out_dir):
    fig, ax = plt.subplots()

    markers, caps, bars = ax.errorbar(df_r_copy['Z'], df_r_copy['log_tau_drw_obs'], \
                                      yerr=[ df_r_copy['log_tau_drw_obs_err_l'], df_r_copy['log_tau_drw_obs_err_u'] ],
                                      fmt='.k')
    
    [bar.set_alpha(.2) for bar in bars]
    
    baseline = df['MJD'][0][-1] - df['MJD'][0][0]
    
    ax.axhline(np.log10(baseline), color='r', ls='-', alpha=.5)
    
    y1, y2 = ax.get_ylim()
    
    ax.set_xlabel('z', fontsize=15)
    ax.set_ylabel(r'$\log_{10} (\tau_{obs})$', fontsize=15)
    
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)    
    
    ax.text(2.5, y2+.05, 
            r'$\log_{10} (\tau_{baseline}) = $ %.3f' % np.log10(baseline),
           fontsize=12)

    ax.tick_params('both', which='major', length=8)
    ax.tick_params('both', which='minor', length=3)
    
    plt.savefig( out_dir + 'redshift_v_tau_' + band + '.pdf')
   
    return fig, ax





def plot_tau_sf(tot_dict, out_dir, levels=5):
    import seaborn as sns

    for b in 'gri':

        sf_inf1 = np.log10( np.sqrt(2)*(10**tot_dict['SDSS'][b]['log_sigma_drw']) )
        tau_rest1 = tot_dict['SDSS'][b]['log_tau_drw_rest']

        sf_inf2 = np.log10( np.sqrt(2)*(10**tot_dict['SDSS+PS1+DES'][b]['log_sigma_drw']) )
        tau_rest2 = tot_dict['SDSS+PS1+DES'][b]['log_tau_drw_rest']

        sf_inf3 = np.log10( np.sqrt(2)*(10**tot_dict['SDSS+PS1'][b]['log_sigma_drw']) )
        tau_rest3 = tot_dict['SDSS+PS1'][b]['log_tau_drw_rest']


        fig, ax = plt.subplots(figsize=(7,5))

        sns.kdeplot(sf_inf1, tau_rest1, gridsize=1000,
                ax=ax, levels=levels,
                linewidths=1, label='SDSS',
                linestyles='--',
                color='red', fill=False)

        sns.kdeplot(sf_inf3, tau_rest3, gridsize=1000,
                ax=ax, levels=levels,
                linewidths=1.5, label='SDSS+PS1',
                linestyles='-.',
                color='black', fill=False)

        sns.kdeplot(sf_inf2, tau_rest2, gridsize=1000,
                ax=ax, levels=levels,
                linewidths=2, label='SDSS+PS1+DES',
                color='blue', fill=False)


#        if b == 'g':
#            sf_inf_c = np.log10( np.sqrt(2)*tot_dict['Colin']['Tot']['sigma_drw']  )
#            tau_rest_c = np.log10( tot_dict['Colin']['Tot']['tau_rest'] )
#
#            sns.kdeplot(sf_inf_c, tau_rest_c, gridsize=1000,
#                    ax=ax, levels=5,
#                    linewidths=2, label='Colin',
#                    color='orange', fill=False)



        ax.set_ylabel(r'$\log_{10}(\tau_{rest})$', fontsize=15)
        ax.set_xlabel(r'$\log_{10}(SF_{\infty})$', fontsize=15)

        ax.set_xlim(-1.75, 0.25)
        ax.set_ylim(.75, 4)

        ax.set_title(b + ' Band Taufit')
        plt.legend(loc='upper left')
        
        plt.savefig(out_dir + 'SF_v_tau_' + b + '.pdf')
        plt.show()
        
    return 






def lambda_plot(tot_dict, survey, colors, ls, out_dir, levels=[.3, .7]):

    import seaborn as sns
    from scipy.optimize import curve_fit
    
    def linear(x, a, b):
        return a*x + b

    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(7,6))
    
    
    dfs = [ tot_dict[survey]['g'], tot_dict[survey]['r'], tot_dict[survey]['i'] ]
    for i, df in enumerate(dfs):
    
        sf_vals = np.log10( np.sqrt(2) * 10**(df['log_sigma_drw']) )
        tau_vals = df['log_tau_drw_rest']
        l_vals =  np.log10(df['lambda_rest']/4000)
        
        sf_err = (df['log_sigma_drw_err_l'] + df['log_sigma_drw_err_u'])/2
        tau_err = (df['log_tau_drw_obs_err_l'] + df['log_tau_drw_obs_err_u'])/2/(1+df['Z'])
        
        sns.kdeplot(l_vals, tau_vals, gridsize=1000,
                    ax=ax[0], levels=levels,
                    linestyles=ls[i], color=colors[i],
                    linewidths=1.5)

        sns.kdeplot(l_vals, sf_vals, gridsize=1000,
                    ax=ax[1], levels=levels,
                    linestyles=ls[i], color=colors[i],
                    linewidths=1.5)
        
        if i == 0:
            tau_tot = tau_vals
            sf_tot = sf_vals
            l_tot = l_vals
            
            sf_err_tot = sf_err
            tau_err_tot = tau_err
        else:
            tau_tot = np.concatenate([tau_tot, tau_vals])
            sf_tot = np.concatenate([sf_tot, sf_vals])
            l_tot = np.concatenate([l_tot, l_vals])
            
            sf_err_tot = np.concatenate([sf_err_tot, sf_err])
            tau_err_tot = np.concatenate([tau_err_tot, tau_err])
    
    x1=-.7
    x2 = .2
    
    output = linmix_regression(l_tot , tau_tot, None, tau_err_tot)
    l_plt = np.linspace(x1, x2, 100)
    up, med, low = linmix_ranges(output, l_plt)

    ax[0].fill_between(l_plt, low, up, color='k', alpha=.25)
    ax[0].plot(l_plt, med, color='k', lw=1)
                            
        

    output = linmix_regression(l_tot , sf_tot, None, sf_err_tot)
    l_plt = np.linspace(x1, x2, 100)
    up, med, low = linmix_ranges(output, l_plt)

    ax[1].fill_between(l_plt, low, up, color='k', alpha=.25)
    ax[1].plot(l_plt, med, color='k', lw=1)
        
        
        
    ax[0].set_ylabel(r'$\log_{10}$( $\tau_{rest}$ [d] )', fontsize=12)
    ax[0].set_xlabel(r'$\log_{10}$( $\lambda_{rest} / 4000$ [$\AA$] )', fontsize=12)

    ax[1].set_ylabel(r'$\log_{10}$( $SF_{\infty}$ [mag] )', fontsize=12)
    ax[1].set_xlabel(r'$\log_{10}$( $\lambda_{rest} / 4000$ [$\AA$] )', fontsize=12)

    ax[1].tick_params(axis='x', labelsize=10)

    ax[0].tick_params(axis='y', labelsize=10)
    ax[1].tick_params(axis='y', labelsize=10)
    
    

    ax[0].set_ylim(1, 4)
    ax[1].set_ylim(-1.4, -.3)
    
    ax[1].set_xlim(x1, x2)

    fig.subplots_adjust(hspace=0)
    fig.suptitle(survey)

    plt.savefig(out_dir + 'lambda_plot.pdf')
    
    return fig, ax




def SF_from_data(x, y, xtot, bin_num=50, log=True):
    tot_baseline = xtot[-1] - xtot[0]
    
    good_ind = np.argsort(x)
    x = x[good_ind]
    y = y[good_ind]
 
    bin_vals = np.zeros(bin_num-1)
    Npairs = np.zeros(bin_num-1)

    if log == True:
        bin_edges = 10**np.linspace(-3, np.log10(tot_baseline), bin_num)
    else:
        bin_edges = np.linspace(0, tot_baseline, bin_num)
        
    bin_centers = []
    for n in range(len(bin_edges)-1):
        bin_centers.append( (bin_edges[n+1] + bin_edges[n])/2)
    
    
    #For all dt bins, gather all dm**2 
    for i in range(len(x)-1):
        for j in range( len(x[i+1:]) ):
            
            for n in range(bin_num-1):
                left = bin_edges[n]
                right = bin_edges[n+1]
                
                if left < np.abs(x[i]-x[j]) <= right:
                    bin_vals[n] +=  (y[j]-y[i])**2 
                    Npairs[n] += 1

    
    out_vals = np.zeros(len(bin_vals))
    #Get rms_dm for each dt bin
    for i in range(len(bin_vals)):
        if Npairs[i] == 0:
            out_vals[i] = np.nan
            continue    
        
        out_vals[i] = np.sqrt( bin_vals[i]/Npairs[i] ) 
    
    return np.array(bin_centers), np.array(out_vals)


def plot_binned_SF(centers, vals, out_dir):
    
    fig, ax = plt.subplots()
    
    ax.errorbar(centers, vals, 
                marker=None, drawstyle='steps-mid', color='k', 
                linewidth=1.0)
    
#    ax.errorbar(centers, vals, fmt='.k', linewidth=1.0)
    
    ax.set_ylabel('SF($\Delta t$)')
    ax.set_xlabel('$\Delta t_{rest}$')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    return fig, ax



def randomize_lc(y,yerr):    
    return np.random.normal(y, yerr)

def get_lc_sample(x, y, yerr, bin_num):

    for _ in range(1000):
        i = np.random.randint(0, len(y))
        j = np.random.randint(0, len(y))
        if i != j and i < j:
            break

    y = randomize_lc(y, yerr)    

    x_samp = x[i:j]
    y_samp = y[i:j]

    c, v = SF_from_data(x_samp, y_samp, x, bin_num=bin_num)
            
    return c, v
    
def bootstrap_SF(x, y, yerr, bin_num, n_samp=500, cores=40, percentages=[16, 84]):
    
    if n_samp < cores:
        cores = n_samp
    
    pool = mp.Pool(cores)
    
    arg1 = []
    arg2 = []
    arg3 = []
    arg4 = []
    
    for _ in range(n_samp):    
        arg1.append(x)
        arg2.append(y)
        arg3.append(yerr)
        arg4.append(bin_num)
        
    args = list(zip(arg1, arg2, arg3, arg4))
    
    output = pool.starmap(get_lc_sample, args)
    pool.close()
    
    samples = []
    for i in range(len(output)):
        samples.append( output[i][1] )
    
    c = output[0][0]
    
    lower = np.nanpercentile(samples, percentages[0], axis=0)
    upper = np.nanpercentile(samples, percentages[1], axis=0)
    med = np.nanmedian(samples, axis=0)
        
    return c, samples, lower, med, upper





def get_SF_data(tot_dict, survey, band, i, out_dir, bin_num=50, n_samp=500):

    x = tot_dict[survey][band]['MJD'][i]
    y = tot_dict[survey][band]['MAG'][i]
    yerr = tot_dict[survey][band]['MAG ERR'][i]

    z = tot_dict[survey][band]['Z'][i]
    dbid = tot_dict[survey][band]['DBID'][i]

    dt, samp, l, m, u = bootstrap_SF(x, y, yerr, bin_num, cores=48, n_samp=500)
    nan_ind = np.isnan(m)

    fig, ax = plt.subplots()

    ax.set_xscale('log')
    ax.set_yscale('log')
    markers, caps, bars = ax.errorbar(dt/(1+z), m, yerr=[m-l,u-m], fmt='.k', 
                                      capsize=2.5, elinewidth=.75, capthick=.75,
                                      ms=5)
#    ax.fill_between(dt[~nan_ind]/(1+z), l[~nan_ind], u[~nan_ind], color='r', alpha=.25)

    import matplotlib.ticker as ticker
    for axis in [ax.xaxis, ax.yaxis]:
        formatter = ticker.FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)

    ax.set_ylabel('SF [mag]', fontsize=15)
    ax.set_xlabel('$\Delta t_{rest}$ [days]', fontsize=13)

    ax.tick_params('both', which='both', labelsize=12)
    ax.tick_params('both', which='major', length=6)
    ax.tick_params('both', which='minor', length=2)

    ax.set_title('DBID: %s' % dbid)

    plt.savefig( out_dir + str(dbid) + '/' + str(survey) + '_SF_' + str(band) + '.pdf')
    plt.close(fig)
    
    return dt/(1+z), m, l, u 






def merge_plots(out_dir,email=None):

    """ Merge all the diagnostic pdf into one pdf """

    import glob
    
    f1 = glob.glob(out_dir + 'light_curve_g*')[0]
    f2 = glob.glob(out_dir + 'light_curve_r*')[0]
    f3 = glob.glob(out_dir + 'light_curve_i*')[0]
    
    f4 = glob.glob(out_dir + 'Total_taufit_result_g*')[0]
    f5 = glob.glob(out_dir + 'Total_taufit_result_r*')[0]
    f6 = glob.glob(out_dir + 'Total_taufit_result_i*')[0]
    
    f7 = glob.glob(out_dir + 'SDSS+PS1+DES_SF_g*')[0]
    f8 = glob.glob(out_dir + 'SDSS+PS1+DES_SF_r*')[0]
    f9 = glob.glob(out_dir + 'SDSS+PS1+DES_SF_i*')[0]
    
    paths = [f1, f2, f3, f4, f5, f6, f7, f8, f9]
    
    from PyPDF2 import PdfFileReader, PdfFileWriter

    pdf_writer = PdfFileWriter()

    for path in paths:
        file_path = path
        if os.path.exists(file_path):
            pdf_reader = PdfFileReader(file_path)
            for page in range(pdf_reader.getNumPages()):
                # Add each page to the writer object
                pdf_writer.addPage(pdf_reader.getPage(page))

    # Write out the merged PDF
    out_file = os.path.join(out_dir,'diagnostic.pdf')
    with open(out_file, 'wb') as out:
        pdf_writer.write(out)

    # send notification and file to the email if provided
    if email is not None:
        userhome = os.path.expanduser('~')
        username = os.path.split(userhome)[-1]
        text_body = 'To retrieve summary PDF: \n scp %s@mimas.astro.illinois.edu:%s .' % (username, out_file)
        os.system('echo "%s" | mail -s "[TEQUILA SHOTS] Runs at %s completed" %s' % (text_body, out_dir, email)) 

    return 0  
