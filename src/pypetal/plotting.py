import corner
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import palettable
from matplotlib.colors import ListedColormap

import pypetal.utils as utils
from pypetal.petalio import err2str

mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.direction'] = 'in'

mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = 'in'

mpl.rcParams["figure.autolayout"] = False

mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['savefig.format'] = 'pdf'





###################################################################
############################  pyCCF ###############################
###################################################################

def plot_pyccf_results(x1, y1, yerr1, x2, y2, yerr2,
                       ccf_lags, ccf,
                       cccd_lags, ccpd_lags,
                       nbin=50, time_unit='d', lc_unit=['mJy', 'mJy'], lc_names=['', ''],
                       fname=None, show=False):

    """Plot the results of using pyCCF on two light curves.

    Parameters
    ----------

    x1 : array-like
        Time array for first light curve.

    y1 : array-like
        Values for the first light curve.

    yerr1 : array-like
        Uncertainty in the first light curve.

    x2 : array-like
        Time array for second light curve.

    y2 : array-like
        Values for the second light curve.

    yerr2 : array-like
        Uncertainty in the second light curve.

    ccf_lags : array-like
        Lags for the CCF.

    ccf : array-like
        CCF values.

    cccd_lags : array-like
        The output simulated lags from the CCCD.

    ccpd_lags : array-like
        The output simulated lags from the CCPD.

    nbin : int, optional
        Number of bins to use for the histograms. Default is 50.

    time_unit : str, optional
        The unit of time for the light curves for the plot. Default is 'd'.

    lc_unit : str, optional
        The unit of the light curves for the plot. Default is 'mJy'.

    lc_names : list, optional
        The names of the light curves to use for the plot. Default is ['', ''].

    fname : str, optional
        If not ``None``, will save the plot to the given filename. Default is ``None``.

    show : bool, optional
        If ``True``, will show the plot. Default is ``False``.

    """


    cent_med = np.median( cccd_lags )
    cent_hi = np.percentile( cccd_lags, 84 )
    cent_lo = np.percentile( cccd_lags, 16 )

    peak_med = np.median( ccpd_lags )
    peak_hi = np.percentile( ccpd_lags, 84 )
    peak_lo = np.percentile( ccpd_lags, 16 )

    min_lag = np.min(ccf_lags)
    max_lag = np.max(ccf_lags)


    lag_unit = time_unit



    #Setup plot axes
    fig = plt.figure( figsize=(12, 7) )
    gs = gridspec.GridSpec(ncols=3, nrows=6)

    ax1 = fig.add_subplot(gs[:3, :2])
    ax2 = fig.add_subplot(gs[3:6, :2], sharex=ax1)

    ax3 = fig.add_subplot(gs[:2, 2])
    ax4 = fig.add_subplot(gs[2:4, 2], sharex=ax3)
    ax5 = fig.add_subplot(gs[4:6, 2], sharex=ax3)



    #Plot light curve data
    _, caps, bars = ax1.errorbar(x1, y1, yerr1, fmt='.k')
    [bar.set_alpha(.25) for bar in bars]

    _, _, bars = ax2.errorbar(x2, y2, yerr2, fmt='.k')
    [bar.set_alpha(.25) for bar in bars]


    #Increase y-bounds to fit text
    ymin, ymax = ax1.get_ylim()
    ymax = ymin + (ymax - ymin)*1.15
    ax1.set_ylim(ymin, ymax)

    ymin, ymax = ax2.get_ylim()
    ymax = ymin + (ymax - ymin)*1.15
    ax2.set_ylim(ymin, ymax)

    ax1.text( .05, .95, lc_names[0], transform=ax1.transAxes, ha='left', va='top', fontsize=17 )
    ax2.text( .05, .95, lc_names[1], transform=ax2.transAxes, ha='left', va='top', fontsize=17 )

    ax2.set_xlabel('Time [' + str(time_unit) + ']', fontsize=19)

    #----------------------------
    #Plot CCF
    ax3.plot( ccf_lags, ccf )
    ax3.set_ylim(-1, 1)

    #----------------------------
    #Plot CCCD
    bin_vals, bin_edges, _ = ax4.hist( cccd_lags, nbin, range=(min_lag, max_lag), density=False )

    #Increase y-bounds to fit text
    ymin, ymax = ax4.get_ylim()
    ymax = ymin + (ymax - ymin)*1.1
    ax4.set_ylim(ymin, ymax)

    res_txt = r'$' + err2str(cent_med, cent_hi-cent_med, cent_med-cent_lo, dec=2) + r'$ ' + lag_unit
    ax4.text(.05, .8, res_txt, transform=ax4.transAxes, fontsize=14)

    #----------------------------
    #Plot CCPD
    bin_vals, _, _ = ax5.hist( ccpd_lags, bins=bin_edges, range=(min_lag, max_lag), density=False )

    #Increase y-bounds to fit text
    ymin, ymax = ax5.get_ylim()
    ymax = ymin + (ymax - ymin)*1.1
    ax5.set_ylim(ymin, ymax)


    res_txt = r'$' + err2str( peak_med, peak_hi-peak_med, peak_med-peak_lo, dec=2 ) + r'$ ' + lag_unit
    ax5.text(.05, .8, res_txt, transform=ax5.transAxes, fontsize=14)
    ax5.set_xlabel('Lag [' + str(lag_unit) + ']', fontsize=19)


    #Add ylabels
    if lc_unit[0] == 'mag':
        ytxt = 'Magnitude'
        ax1.invert_yaxis()

    elif lc_unit[0] == 'Arbitrary Units':
        ytxt = 'Flux'
    else:
        ytxt = 'Flux [' + str(lc_unit[0]) + ']'
    ax1.set_ylabel( ytxt, fontsize=19, va='center' )


    if lc_unit[1] == 'mag':
        ytxt = 'Magnitude'
        ax2.invert_yaxis()

    elif lc_unit[1] == 'Arbitrary Units':
        ytxt = 'Flux'
    else:
        ytxt = 'Flux [' + str(lc_unit[1]) + ']'
    ax2.set_ylabel( ytxt, fontsize=19, va='center' )



    plt.figtext( .94, .72, 'CCF', rotation=270, fontsize=19 )
    plt.figtext( .94, .45, 'CCCD', rotation=270, fontsize=19 )
    plt.figtext( .94, .2, 'CCPD', rotation=270, fontsize=19 )


    for i, ax in enumerate([ax3, ax4, ax5]):
        ax.yaxis.tick_right()
        ax.yaxis.set_ticks_position('both')

        if i == 2:
            ax.tick_params('both', labelsize=12)
        else:
            ax.tick_params('x', labelsize=0)
            ax.tick_params('y', labelsize=12)

        ax.tick_params('both', which='major', length=6)
        ax.tick_params('both', which='minor', length=3)

        ax.locator_params('x', nbins=5)
        ax.locator_params('y', nbins=3)


    for i, ax in enumerate([ax1, ax2]):

        if i == 1:
            ax.tick_params('both', labelsize=13)
        else:
            ax.tick_params('x', labelsize=0)
            ax.tick_params('y', labelsize=13)

        ax.tick_params('both', which='major', length=7)
        ax.tick_params('both', which='minor', length=3)

        ax.locator_params('y', nbins=6)


    plt.subplots_adjust(hspace=.07, wspace=.1)

    if fname is not None:
        plt.savefig(fname, dpi=300, bbox_inches='tight')

    if show:
        plt.show()

    plt.cla()
    plt.clf()
    plt.close()

    del fig, ax1, ax2, ax3, ax4, ax5, gs




###################################################################
############################# pyZDCF ##############################
###################################################################

def plot_pyzdcf_results(x1, y1, yerr1, x2, y2, yerr2,
                        dcf_df,
                        plike_dict=None,
                        time_unit='d', lc_unit=['mJy', 'mJy'], lc_names=['', ''],
                        fname=None, show=False):

    """Plot the results of pyZDCF.


    Parameters
    ----------

    x1 : list of float
        Time array for the first light curve.

    y1 : list of float
        Values for the first light curve.

    yerr1 : list of float
        Uncertainty in the first light curve.

    x2 : list of float
        Time array for the second light curve.

    y2 : list of float
        Values for the second light curve.

    yerr2 : list of float
        Uncertainty in the second light curve.

    dcf_df : pandas.DataFrame
        The output ``DataFrame`` object from pyZDCF.

    plike_dict : dict, optional
        The output dict from PLIKE, if it is run. Default is ``None``.

    time_unit : str, optional
        The unit of time for the light curves. Default is 'd'.

    lc_unit : str, optional
        The unit of the light curves. Default is 'mJy'.

    lc_names : list of str, optional
        The names of the light curves. Default is ['', ''].

    fname : str, optional
        If not ``None``, will save the plot to the given filename. Default is ``None``.

    show : bool, optional
        If ``True``, will show the plot. Default is ``False``.

    """


    if plike_dict is None:
        include_plike = False
    else:
        include_plike = True

    #Plot the CCF
    fig = plt.figure(  figsize=(12,4))
    gs = gridspec.GridSpec( ncols=3, nrows=2 )


    ax1 = fig.add_subplot( gs[0,:2] )
    ax2 = fig.add_subplot( gs[1,:2], sharex=ax1 )
    ax3 = fig.add_subplot( gs[:,2] )

    lines, caps, bars = ax1.errorbar(x1, y1, yerr1, fmt='.k')
    [bar.set_alpha(.2) for bar in bars]

    ax1.tick_params('x', labelsize=0)


    lines, caps, bars = ax2.errorbar( x2, y2, yerr2, fmt='.k' )
    [bar.set_alpha(.2) for bar in bars]

    ax2.set_xlabel('Time [' + str(time_unit) + ']', fontsize=17)



    lags = dcf_df.tau
    vals = dcf_df.dcf

    err_lo = dcf_df['-err(dcf)']
    err_hi = dcf_df['+err(dcf)']

    tau_err_lo = dcf_df['-sig(tau)']
    tau_err_hi = dcf_df['+sig(tau)']

    lines, caps, bars = ax3.errorbar(lags, vals, yerr=[err_lo, err_hi],
                                     xerr=[tau_err_lo, tau_err_hi], fmt='.k')
    [bar.set_alpha(.2) for bar in bars]


    if include_plike:
        ax3.axvline( plike_dict['ML_lag'], c='b' )
        ax3.axvspan( plike_dict['ML_lag'] - plike_dict['ML_lag_err_lo'],
                    plike_dict['ML_lag'] + plike_dict['ML_lag_err_hi'],
                    color='b', alpha=.1)


    ax3.set_xlabel('Lag [' + str(time_unit) + ']', fontsize=17)

    plt.figtext(.95, .5, 'ZDCF', fontsize=17, rotation=270, va='center')

    ax3.yaxis.tick_right()
    ax3.yaxis.set_ticks_position('both')

    for i, ax in enumerate([ax1, ax2]):
        y1, y2 = ax.get_ylim()
        y2 = y1 + (y2-y1)*1.3
        ax.set_ylim(y1, y2)

        ax.text( .04, .8, lc_names[i], transform=ax.transAxes, fontsize=15 )


    for ax in [ax1, ax2, ax3]:
        ax.tick_params('both', which='major', length=7)
        ax.tick_params('both', which='minor', length=3)


    for ax in [ax1, ax2]:
        ax.tick_params('y', labelsize=13)
        ax.tick_params('x', labelsize=12)

    ax3.tick_params('both', labelsize=12)
    ax3.locator_params('x', nbins=5)
    ax3.locator_params('y', nbins=5)



    #Add ylabels
    if lc_unit[0] == 'mag':
        ytxt = 'Magnitude'
        ax1.invert_yaxis()

    elif lc_unit[0] == 'Arbitrary Units':
        ytxt = 'Flux'
    else:
        ytxt = 'Flux [' + str(lc_unit[0]) + ']'
    plt.figtext( .06, .68, ytxt, fontsize=16, va='center', rotation=90 )


    if lc_unit[1] == 'mag':
        ytxt = 'Magnitude'
        ax2.invert_yaxis()

    elif lc_unit[1] == 'Arbitrary Units':
        ytxt = 'Flux'
    else:
        ytxt = 'Flux [' + str(lc_unit[1]) + ']'
    plt.figtext( .06, .3, ytxt, fontsize=16, va='center', rotation=90 )

    plt.subplots_adjust(hspace=.03, wspace=.06)

    if fname is not None:
        plt.savefig(fname, dpi=300, bbox_inches='tight')

    if show:
        plt.show()

    plt.cla()
    plt.clf()
    plt.close()

