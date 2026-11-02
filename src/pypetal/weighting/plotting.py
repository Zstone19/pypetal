import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import (BboxConnector,
                                                   TransformedBbox, inset_axes)

from pypetal.utils.petalio import err2str

mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.size'] = 3

mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.size'] = 3

mpl.rcParams["figure.autolayout"] = False

mpl.rcParams['savefig.dpi'] = 200
mpl.rcParams['savefig.format'] = 'pdf'



def plot_weights(output_dir, line_name, res, n0, k, time_unit='d', plot=False):

    """Plot the weight distribution for a given object for a given module's output distribution.

    Parameters
    ----------

    output_dir : str
        Directory to save plots to.

    line_name : str
        Name of the line.

    res : dict
        Dictionary containing the results from the module.

    n0 : int
        The number of overlapping data points at a lag tau = 0.

    k : int, optional
        The power the probability distribution was raised to.

    time_unit : str, optional
        The time unit of the lag distribution. Default is 'd'.

    plot : bool, optional
        Whether to show the plotted results. Default is ``False``.



    Returns
    -------

    fig : matplotlib.figure.Figure
        The figure object containing the plot.

    ax : matplotlib.axes.Axes
        The axes object for the plot.

    """


    lags = res['lags'][-1].copy()
    acf = res['acf'][-1].copy()
    prob_dist = res['weight_dist'][-1].copy()
    ntau = res['ntau'][-1].copy()

    #Get P(tau)
    probs = (ntau/n0)**k

    #Set all ACF values past first negative to zero
        #Left side
    left_inds = np.argwhere( (acf < 0) & (lags < 0) ).T[0]
    if len(left_inds) == 0:
        ind1 = 0
    else:
        ind1 = np.max( left_inds )

        #Right side
    right_inds = np.argwhere( (acf < 0) & (lags > 0) ).T[0]
    if len(right_inds) == 0:
        ind2 = len(acf)-1
    else:
        ind2 = np.min( right_inds )

    acf_new = acf.copy()
    acf_new[:ind1+1] = 0
    acf_new[ind2:] = 0



    fig, ax = plt.subplots()

    ax.plot(lags, probs, c='gray', label=r'P($\tau$) = $[N(\tau) \ / \ N(0)]^{' + '{}'.format(k) + '}$')
    ax.plot(lags, acf_new, c='r', lw=1.5, label='ACF')
    ax.plot( lags, prob_dist/np.max(prob_dist), c='b', lw=2, label = r'P($\tau$) * ACF' )

    ax.set_xlabel('Lag [' + str(time_unit) + ']', fontsize=16)
    ax.set_ylabel('Weights', fontsize=16)

    ax.tick_params(labelsize=11)
    ax.tick_params('both', which='major', length=8)
    ax.tick_params('both', which='minor', length=4)

    plt.figlegend(bbox_to_anchor=(1.3,.9), fontsize=12)
    plt.savefig( output_dir + line_name + r'/weights/' + line_name + '_weights.pdf', dpi=200, bbox_inches='tight' )

    if plot:
        plt.show()

    plt.cla()
    plt.clf()
    plt.close()

    return fig, ax



def plot_weight_output(output_dir, cont_fname, line_fnames, line_names,
                 res, general_kwargs, module, zoom=False, plot=True):

    """Plot the output of the `run_weights` function. Will output histograms for each line
    of the given time lag distribution for the module, after using the weighting procedure
    described by Grier et al. (2019).
    The plot will have a panel for each line. Each panel will contain three subplots within them.
    The top plot will have the weighting function w($\tau$), the ACF of the continuum, and the smoothed
    distribution used to find the peak. The main panel will have the original distribution and the weighted
    distribution after applying w($\tau$). Additionally, there can be an inset plot zooming in on
    the peak of the distribution (see ``zoom``).


    Parameters
    ----------
    output_dir : str
        The directory to output the plots to.
    cont_fname : str
        The name of the continuum file.
    line_fnames : list of str
        The names of the line files.
    line_names : list of str
        The names of the light curves.
    res : dict
        The output of the `run_weights` function for a given module. For example, if ``res``
        is the output of ``run_weights`` when ``run_pyccf=True``, then this input should be
        ``res['pyccf']`` to plot the pyCCF results.
    general_kwargs : dict
        The general keyword arguments used in ``run_pipeline`` function.
    module : str
        The name of the module this function is plotting the results for. Should be either
        'pyccf' or 'javelin'.
    zoom : bool, optional
        If ``True``, the plots will be zoomed in on the peak of the distribution. Default is
        ``False``.
    plot : bool, optional
        If ``True``, the plots will be displayed. Default is ``True``.
    .. note:: Even if ``zoom=False``, the plots will show an inset on the peak of the distribution if the range of the distribution is too large.
    """


    #Read general kwargs
    plot = general_kwargs['plot']
    time_unit = general_kwargs['time_unit']

    nlc = len(res['lags'])

    Ncol = 3
    Nrow = nlc//Ncol + 1

    #--------------------------------------------------------------------------------

    if nlc <= 3:
        Nrow = 1
        Ncol = nlc


    fig = plt.figure(figsize=( 6*Ncol, 5*Nrow ))

    gs = gridspec.GridSpec(Nrow, Ncol, figure=fig)
    ax_tot = np.zeros( (Nrow, Ncol), dtype=object )

    for i in range(Nrow):
        for j in range(Ncol):
            sub_gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[i, j],
                                                    height_ratios=[1, 3], hspace=0)

            ax_bot = plt.subplot(sub_gs[1])
            ax_top = plt.subplot(sub_gs[0], sharex=ax_bot)

            ax_tot[i, j] = [ax_top, ax_bot]


    for i in range(nlc):

        if module == 'pyccf':
            lags = res['lags'][i].copy()
            lag_dist = res['CCCD'][i].copy()

            lag_err_lo = res['centroid'][i][0]
            lag_err_hi = res['centroid'][i][2]
            lag_value = res['centroid'][i][1]

            xlabel = r'{ \rm CCCD }'

        elif module == 'javelin':
            lags = res['lags'][i].copy()
            lag_dist = res['lag_dist'][i].copy()

            lag_err_lo = res['tophat_lag'][i][0]
            lag_err_hi = res['tophat_lag'][i][2]
            lag_value = res['tophat_lag'][i][1]

            xlabel = r't'

        elif module == 'pyroa':
            lags = res['lags'][i].copy()
            lag_dist = res['lag_dist'][i].copy()

            lag_err_lo = res['time_delay'][i][0]
            lag_err_hi = res['time_delay'][i][2]
            lag_value = res['time_delay'][i][1]

            xlabel = r'\tau'

        elif module == 'mica2':
            lags = res['lags'][i].copy()
            lag_dist = res['lag_dist'][i].copy()

            lag_err_lo = res['time_lag'][i][0]
            lag_err_hi = res['time_lag'][i][2]
            lag_value = res['time_lag'][i][1]

            xlabel = r'\tau'


        llim = res['bounds'][i][0]
        rlim = res['bounds'][i][2]
        peak = res['bounds'][i][1]

        weight_dist = res['weight_dist'][i].copy()
        smooth_dist = res['smoothed_dist'][i].copy()
        acf = res['acf'][i].copy()


        #Set ACF=0 when ACF<0
            #Left side
        left_inds = np.argwhere( (acf < 0) & (lags < 0) ).T[0]
        if len(left_inds) == 0:
            ind1 = 0
        else:
            ind1 = np.max( left_inds )

            #Right side
        right_inds = np.argwhere( (acf < 0) & (lags > 0) ).T[0]
        if len(right_inds) == 0:
            ind2 = len(acf)-1
        else:
            ind2 = np.min( right_inds )

        acf[:ind1] = 0
        acf[ind2:] = 0



        if nlc == 1:
            col_ind = 0
            row_ind = 0
        else:
            col_ind = i%Ncol
            row_ind = i//Ncol


        #Plot original lag distribution
        dbin = lags[1] - lags[0]
        bin0 = np.min(lags)

        bin_edges = []
        bin0 = np.min(lags) - dbin/2

        for j in range(len(lags)+1):
            bin_edges.append( bin0 + j*dbin )

        hist, _ = np.histogram(lag_dist, bins=bin_edges)
        hist = np.array(hist, dtype=float)


        #Set top of plot
        good_ind = np.argwhere( ( lags >= llim ) & ( lags <= rlim ) ).T[0]

        if np.argmax(hist) in good_ind:
            ytop = 1.5*np.max(hist)

        elif np.percentile(hist, 99) > np.max(hist):
            ytop = np.percentile(hist, 99)

        else:
            ytop = 1.5*np.max(hist[good_ind])

        #Inset axis?
        lagmax = np.max(lags)
        lagmin = np.min(lags)

        if ( (rlim-llim) < .1*(lagmax-lagmin) ) | zoom:

            #Put inset on the right if peak is on the left and vice-versa
            if peak < (lagmin+lagmax)/2:
                loc=1
            else:
                loc=2

            axins = inset_axes(ax_tot[row_ind, col_ind][1], width='45%', height='40%', loc=loc)

            #Make connectors between the inset and the main plot
            good_ind = np.argwhere( (lags >= llim) & (lags <= rlim) ).T[0]
            if np.max(hist[good_ind]) < .4*ytop:
                loc11 = 3
                loc21 = 2

                loc12 = 4
                loc22 = 1

            else:
                loc11 = 3
                loc21 = 3

                loc12 = 4
                loc22 = 4

            rect = TransformedBbox(axins.viewLim, ax_tot[row_ind, col_ind][1].transData)

            p1 = BboxConnector(axins.bbox, rect, loc1=loc11, loc2=loc21, ec='k')
            axins.add_patch(p1)
            p1.set_clip_on(False)

            p2 = BboxConnector(axins.bbox, rect, loc1=loc12, loc2=loc22, ec='k')
            axins.add_patch(p2)
            p2.set_clip_on(False)

            #Plot the histograms on both the main and inset axes
            for a in [ax_tot[row_ind, col_ind][1], axins]:
                bin1 = a.fill_between(lags, np.zeros_like(hist), hist, step='mid')

                bin2 = a.fill_between(lags, np.zeros_like(hist), weight_dist*hist, step='mid', color='r', alpha=.7)

                a.axvspan( llim, rlim, color='k', alpha=.1 )
                a.set_ylim(0, ytop)

            axins.set_xlim( llim, rlim )
            axins.set_xticks([])


            #If the peak is small, make box around it in the main plot
            if np.max(hist[good_ind]) < .4*ytop:
                y2 = np.max( hist[good_ind] )
                axins.set_ylim(top=1.05*y2)

                x1, x2 = ax_tot[row_ind, col_ind][1].get_xlim()
                xmin = (llim - x1)/(x2-x1)
                xmax = (rlim - x1)/(x2-x1)

                y1_ax, y2_ax = ax_tot[row_ind, col_ind][1].get_ylim()
                ymax = (1.1*y2 - y1_ax)/(y2_ax-y1_ax)

                ax_tot[row_ind, col_ind][1].axhline( 1.05*y2, xmin, xmax, color='k', lw=1.5 )
                ax_tot[row_ind, col_ind][1].axvline( llim, 0, ymax, color='k', lw=1.5 )
                ax_tot[row_ind, col_ind][1].axvline( rlim, 0, ymax, color='k', lw=1.5 )


                yticks = ax_tot[row_ind, col_ind][1].get_yticks()
                m_yticks = ax_tot[row_ind, col_ind][1].get_yticks(minor=True)

                good_ind = np.argwhere( yticks <= 1.05*y2 ).T[0]
                axins.set_yticks( yticks[good_ind] )

                good_ind = np.argwhere( m_yticks <= 1.05*y2 ).T[0]
                axins.set_yticks( m_yticks[good_ind], minor=True )

                axins.tick_params('both', which='major', length=7)
                axins.tick_params('both', which='minor', length=4)

            else:
                y1, y2 = ax_tot[row_ind, col_ind][1].get_ylim()
                axins.set_ylim(y1, y2)

                yticks = ax_tot[row_ind, col_ind][1].get_yticks()
                m_yticks = ax_tot[row_ind, col_ind][1].get_yticks(minor=True)

                axins.set_yticks( yticks )
                axins.set_yticks( m_yticks, minor=True )

                axins.tick_params('both', which='major', length=0)
                axins.tick_params('both', which='minor', length=0)


            axins.set_yticklabels([])



        else:
            loc=1

            bin1 = ax_tot[row_ind, col_ind][1].fill_between(lags, np.zeros_like(hist), hist, step='mid')
            bin2 = ax_tot[row_ind, col_ind][1].fill_between(lags, np.zeros_like(hist), weight_dist*hist, step='mid', color='r', alpha=.7)
            ax_tot[row_ind, col_ind][1].axvspan( llim, rlim, color='k', alpha=.1 )

            ax_tot[row_ind, col_ind][1].set_ylim(0, ytop)


        #Plot ACF
        im3, = ax_tot[row_ind, col_ind][0].plot(lags, acf, c='r')


        #Plot weighting function
        im1, = ax_tot[row_ind, col_ind][0].plot(lags, weight_dist, c='k')
        ax_tot[row_ind, col_ind][0].set_yticks([.5, 1])

        #Plot smoothed distribution
        im2, = ax_tot[row_ind, col_ind][0].plot(lags, smooth_dist/np.max(smooth_dist), c='DodgerBlue')
        ax_tot[row_ind, col_ind][0].axvspan( llim, rlim, color='k', alpha=.1 )

        #Write lag and error
        peak_str = err2str( lag_value, lag_err_hi, lag_err_lo, dec=2 )
        peak_str = r'$' + r'{}'.format(peak_str) + r'$' + ' ' + time_unit

        if loc == 1:
            xtxt = .05
            ha='left'
        else:
            xtxt = .95
            ha='right'

        ax_tot[row_ind, col_ind][1].text( xtxt, .85, peak_str,
                                    ha=ha, transform=ax_tot[row_ind, col_ind][1].transAxes,
                                    fontsize=13 )

        ax_tot[row_ind, col_ind][1].set_xlabel( r'$' + xlabel + r'_{' + r'{}'.format(line_names[i+1]) + r'}$ [' + time_unit + ']', fontsize=15 )


    for i in range(Nrow):
        ax_tot[i, 0][1].set_ylabel('N', fontsize=16)

    for i in range(Nrow):
        for j in range(Ncol):
            ax_tot[i,j][1].tick_params('both', labelsize=11)
            ax_tot[i,j][1].tick_params('both', which='major', length=7)
            ax_tot[i,j][1].tick_params('both', which='minor', length=3)

            ax_tot[i,j][0].tick_params('y', labelsize=11)
            ax_tot[i,j][0].tick_params('x', labelsize=0)
            ax_tot[i,j][0].tick_params('both', which='major', length=7)
            ax_tot[i,j][0].tick_params('both', which='minor', length=3)


    plt.subplots_adjust(wspace=.15)

    #Put legends for the top and bottom plots
    for i in range(Nrow):
        ax_tot[i, -1][0].legend( [im1, im2, im3], [r'w($\tau$)', 'Smoothed Dist', 'ACF'],
                                bbox_to_anchor=(1,1.1), fontsize=11, loc='upper left' )

        ax_tot[i, -1][1].legend( [bin1, bin2], ['Original', 'Weighted'],
                                bbox_to_anchor=(1, 1), fontsize=11, loc='upper left')


    plt.savefig( output_dir + module + '_weights_res.pdf', dpi=200, bbox_inches='tight' )


    if plot:
        plt.show()

    plt.cla()
    plt.clf()
    plt.close()

    return
