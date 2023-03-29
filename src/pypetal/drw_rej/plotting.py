import corner
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from matplotlib import gridspec

from pypetal.drw_rej.utils import psd_data, smoothly_broken_power_law

def plot_outcome(x, y, yerr, samples, gp, unit, nsig=0,
                 target=None, show_mean=True,
                 filename=None, jitter=True, show=False):

    """Generate summary image for the DRW fit. Contains a plot of the input light curve along with the DRW fit, a corner plot of the DRW parameters, and a plot of the PSD (data-based and fit-based).


    Parameters
    ----------

    x : list of astropy.units.Quantity
        The times of the light curve.

    y : list of astropy.units.Quantity
        The light curve.

    yerr : list of astro.units.Quantity
        The uncertainty in the light curve.

    samples : array_like
        The output MCMC samples from the DRW fit.

    gp : celerite.GP
        The GP object used to fit the DRW

    unit : str, astropy.units.Unit
        The unit of the light curve.

    nsig : int, optional
        The number of standard deviations away from the mean DRW fit at which to
        consider a point an outlier. Default is 0 (i.e. no points are outliers).

    target : str, optional
        The name of the target. Default is ``None``.

    show_mean : bool, optional
        If True, will plot the mean of the DRW fit. Default is True.

    filename : str, optional
        The name of the file to save the image to. Default is ``None``.

    jitter : bool, optional
        Whether or not the jitter term was included in the DRW fit. Default is True.

    show : bool, optional
        Whether or not to show the image. Default is False.




    Returns
    -------

    fig : matplotlib.figure.Figure
        The figure containing the summary image

    axs : list of matplotlib.axes.Axes
        The axes of the figure

    """

    #Plot probability dists of output params

    #tau = 1/c
    #sig = sqrt(a/2)

    baseline = x[-1] - x[0]
    extra_t = int(baseline.value//10)

    if unit == u.mag:
        unit_label = 'Magnitude'
    elif unit == u.dimensionless_unscaled:
        unit_label = 'Flux'
    else:
        unit_label = 'Flux [' + str(unit) + ']'

    if target is None:
        target_label = ''
    else:
        target_label = 'Light Curve of ' + target

    tau_vals = 1/np.exp(samples[:, 1])
    sig_vals = np.sqrt( np.exp(samples[:, 0])/2 )

    if jitter is True:
        n=3
        jitter_vals = np.exp(samples[:, 2])
        sample_vals = np.vstack((np.log10(sig_vals), np.log10(tau_vals), np.log10(jitter_vals) )).T
        labels = [r'$\log_{10}\ (\sigma_{\rm DRW})$',
                  r'$\log_{10}\ (\tau_{\rm DRW})$',
                  r'$\log_{10}\ (\sigma_n)$']
    else:
        n=2
        sample_vals = np.vstack((np.log10(sig_vals), np.log10(tau_vals) )).T
        labels = [r'$\log_{10}\ (\sigma_{\rm DRW})$',
                  r'$\log_{10}\ (\tau_{\rm DRW})$']


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

    if jitter is True:
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
    if jitter is True:
        axs[2,2].title.set_text(jit_title)
        axs[2,2].title.set_fontsize(17)

    #Red out bad tau regions
    axs[1,1].axvspan(np.log10(0.2*baseline.value), np.log10(100*baseline.value), color= "red", zorder=-5, alpha=0.2)

    if jitter is True:
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

    flag_mask = np.abs(y.value - mu_flag) > nsig*yerr.value
    ax1.errorbar( (x.value - x[0].value)[flag_mask],
                 y.value[flag_mask], yerr.value[flag_mask],
                 fmt='.', color='DodgerBlue', capsize=1., alpha=.4, ms=7)
    ax1.errorbar( (x.value - x[0].value)[~flag_mask],
                y.value[~flag_mask], yerr.value[~flag_mask],
                fmt='.k', capsize=1., alpha=.7, ms=8)


    if show_mean is True:
        ax1.plot(t-x[0].value, mu, color='orange', zorder=-1)
    ax1.fill_between(t-x[0].value, mu+std, mu-std, color='orange', alpha=.3)

    ax1.set_xlabel('Time [' + str(x.unit) + ']', fontsize=20)
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

    sbpl_dat = smoothly_broken_power_law(fLS.value, fit_vals[0], fit_vals[1], fit_vals[2], fit_vals[3])
    ax2.plot(fLS.value, sbpl_dat, color='red', label='SBPL Fit')

    #--------------------------------------------------------------------------
    #Plot values for f_br

    f_br = fit_vals[1]
    f_br_err = fit_err[1]

    f_br_plus = f_br + f_br_err
    f_br_minus = f_br - f_br_err

    if f_br_minus < 0.0 and np.isfinite(f_br_plus) is True:
        ls_diff = np.log10(f_br_plus) - np.log10(f_br)
        f_br_minus = 10**( np.log10(f_br) - ls_diff )


    f_br_val = psd_credint[:, 2][0]

    if np.isfinite(f_br_plus) is True:
        #Dat for f_br line
        t = np.logspace(np.log10(f_br_minus), np.log10(f_br_plus), 100)
        dat = np.full(len(t), f_br_val)
        ax2.plot(t, dat, color='red', linewidth=2)

        arrow_head = (f_br, 10**(np.log10(f_br_val)+.2) )
        arrow_base = (f_br, 10**(np.log10(f_br_val)+.6) )
        ax2.errorbar([f_br], [arrow_base[1]], yerr=[arrow_base[1] - arrow_head[1]],
                     uplims=True, color='red', elinewidth=2, capthick=2)

    #--------------------------------------------------------------------------
    #Red out bad freqs

    ax2.axvspan( np.min(fLS.value) , 1/(2*np.pi*0.2*baseline.value), color='red', alpha=.3)

    fmax = (1./(2*np.pi*np.mean(np.diff(x)[np.diff(x)>0])))
    ax2.axvspan(fmax.value, np.max(fLS.value), color='red', alpha=0.2)


    ax2.set_xlabel(r'Frequency [' + str(x.unit) + '$^{-1}$]', fontsize=18)
    ax2.set_ylabel(r'Power', fontsize=18)

    ax2.set_ylim( np.min(powerLS.value) )
    ax2.set_xlim(np.min(fLS.value), np.max(fLS.value))

    ax2.set_xscale('log', nonpositive='clip')
    ax2.set_yscale('log', nonpositive='clip')

    ax2.tick_params('both', labelsize=15)
    ax2.tick_params('both', which='both', direction='in')
    ax2.tick_params('both', which='major', length=8)
    ax2.tick_params('both', which='minor', length=3)

    fig = plt.gcf()

    if filename is not None:
        fig.savefig(f'{filename}', bbox_inches='tight')

    if show is True:
        plt.show()

    return fig, axs