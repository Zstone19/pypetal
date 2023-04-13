import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from pypetal.utils.petalio import err2str

mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.direction'] = 'in'

mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = 'in'

mpl.rcParams["figure.autolayout"] = False

mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['savefig.format'] = 'pdf'





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

        lag_txt = err2str(plike_dict['ML_lag'], plike_dict['ML_lag_err_hi'], plike_dict['ML_lag_err_lo'])
        lag_txt = r'$\tau_{\rm ML} = ' + lag_txt + '$'
        ax2.text( .04, .6, lag_txt, transform=ax2.transAxes, fontsize=15, color='b')



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

    elif (lc_unit[0] == 'Arbitrary Units') or (lc_unit[0] == ''):
        ytxt = 'Flux'
    else:
        ytxt = 'Flux [' + str(lc_unit[0]) + ']'
    plt.figtext( .06, .68, ytxt, fontsize=16, va='center', rotation=90 )


    if lc_unit[1] == 'mag':
        ytxt = 'Magnitude'
        ax2.invert_yaxis()

    elif (lc_unit[1] == 'Arbitrary Units') or (lc_unit[1] == ''):
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

    return
