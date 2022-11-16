import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

from matplotlib.colors import ListedColormap
import palettable

import corner
import utils

import matplotlib as mpl

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
                       nbin=50, time_unit='d', lc_unit='mJy', lc_names=['', ''],
                       plot_weights=False, fname=None, show=False):

    cent_med = np.median( cccd_lags )
    cent_hi = np.percentile( cccd_lags, 84 )
    cent_lo = np.percentile( cccd_lags, 16 )

    peak_med = np.median( ccpd_lags )
    peak_hi = np.percentile( ccpd_lags, 84 )
    peak_lo = np.percentile( ccpd_lags, 16 )
    
    min_lag = np.min(ccf_lags)
    max_lag = np.max(ccf_lags)
    
    
    
    if (cent_med < .5) & (peak_med < .5) & (time_unit == 'd'):
        cent_med *= 24
        cent_hi *= 24
        cent_lo *= 24
        
        peak_med *= 24
        peak_hi *= 24
        peak_lo *= 24
        
        lag_unit = 'hr'
    else:
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
    
    
    if lc_unit == 'mag':
        ytxt = 'Magnitude'
        
        ax1.invert_yaxis()
        ax2.invert_yaxis()
    elif lc_unit == 'Arbitrary Units':
        ytxt = 'Flux'
    else:
        ytxt = 'Flux [' + str(lc_unit) + ']'
    
    
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

    ax4.text(.05, .8, 
            '${:.2f} ^'.format(cent_med) + '{+' + '{:.2f}'.format(cent_hi-cent_med)  + '}_{-' + '{:.2f}'.format(cent_med-cent_lo) + '}$ ' + lag_unit,
            transform=ax4.transAxes,
            fontsize=14)
    
    if plot_weights:
        bin_cents = np.zeros( len(bin_edges)-1 )
        for i in range(len(bin_cents)):
            bin_cents[i] = (bin_edges[i] + bin_edges[i+1])/2
        
        prob_weights, _ = utils.prob_tau(x1, x2, cccd_lags, lagvals=bin_cents)
    
        ax4.fill_between( bin_cents, np.zeros(len(bin_cents)),
                          bin_vals*prob_weights, step='mid', color='r', alpha=.8 )
        ax4.plot( bin_cents, bin_vals*prob_weights, ds='steps-mid', c='k' )
        
        weighted_lags = utils.make_mc_from_weights(x1, x2, cccd_lags, bin_cents)
        

        weighted_med = np.median( weighted_lags )
        weighted_hi = np.percentile( weighted_lags, 84 )
        weighted_lo = np.percentile( weighted_lags, 16 )
        
        ax4.text(.05, .6, 
                    '${:.2f} ^'.format(weighted_med) + '{+' + '{:.2f}'.format(weighted_hi-weighted_med)  + '}_{-' + '{:.2f}'.format(weighted_med-weighted_lo) + '}$ ' + lag_unit,
                    transform=ax4.transAxes,
                    fontsize=14, color='r')

    #----------------------------
    #Plot CCPD
    bin_vals, _, _ = ax5.hist( ccpd_lags, bins=bin_edges, range=(min_lag, max_lag), density=False )

    #Increase y-bounds to fit text
    ymin, ymax = ax5.get_ylim()
    ymax = ymin + (ymax - ymin)*1.1
    ax5.set_ylim(ymin, ymax)

    ax5.text(.05, .8, 
            '${:.2f} ^'.format(peak_med) + '{+' + '{:.2f}'.format(peak_hi-peak_med)  + '}_{-' + '{:.2f}'.format(peak_med-peak_lo) + '}$ ' + lag_unit,
            transform=ax5.transAxes,
            fontsize=14)
    ax5.set_xlabel('Lag [' + str(lag_unit) + ']', fontsize=19)
    
    
    if plot_weights:
        prob_weights, _ = utils.prob_tau(x1, x2, ccpd_lags, lagvals=bin_cents)
    
        ax5.fill_between( bin_cents, np.zeros(len(bin_cents)),
                          bin_vals*prob_weights, step='mid', color='r', alpha=.8 )
        ax5.plot( bin_cents, bin_vals*prob_weights, ds='steps-mid', c='k' )
        
        weighted_lags = utils.make_mc_from_weights(x1, x2, ccpd_lags, bin_cents)
        
        weighted_med = np.median( weighted_lags )
        weighted_hi = np.percentile( weighted_lags, 84 )
        weighted_lo = np.percentile( weighted_lags, 16 )
        
        ax5.text(.05, .6, 
                    '${:.2f} ^'.format(weighted_med) + '{+' + '{:.2f}'.format(weighted_hi-weighted_med)  + '}_{-' + '{:.2f}'.format(weighted_med-weighted_lo) + '}$ ' + lag_unit,
                    transform=ax5.transAxes,
                    fontsize=14, color='r')


    plt.figtext( .06, .5, ytxt, fontsize=19, rotation=90, va='center' )

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
                        time_unit='d', lc_unit='mJy', lc_names=['', ''],
                        fname=None, show=False):

    if plike_dict is None:
        include_plike = False
    else:
        include_plike = True

    #Plot the CCF
    fig = plt.figure(  figsize=(12,4))
    gs = gridspec.GridSpec( ncols=3, nrows=2 )


    ax1 = fig.add_subplot( gs[0,:2] )
    ax2 = fig.add_subplot( gs[1,:2] )
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


    if lc_unit == 'mag':
        ytxt = 'Magnitude'
        
        ax1.invert_yaxis()
        ax2.invert_yaxis()
    elif lc_unit == 'Arbitrary Units':
        ytxt = 'Flux'
    else:
        ytxt = 'Flux [' + str(lc_unit) + ']'


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

    plt.figtext( .07, .5, ytxt, fontsize=17, rotation=90, va='center' )
    plt.subplots_adjust(hspace=.03, wspace=.06)

    if fname is not None:
        plt.savefig(fname, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    
    plt.cla()
    plt.clf()
    plt.close()




###################################################################
############################ JAVELIN ##############################
###################################################################



def plot_javelin_hist(res, fixed=None, nbin=50, 
                      time_unit='d', lc_unit='mag',
                      plot_weights=False, remove_fixed=True, fname=None):


    if lc_unit == 'Arbitrary Units':
        lc_unit_txt = ''
    else:
        lc_unit_txt = '[' + str(lc_unit) +']'
    
    time_unit_txt = '[' + str(time_unit) + ']'




    Ncol = 3
    Nrow = len(res['tophat_params'])//3 + 1

    names = res['tot_dat'].names[1:]
    x1 = res['tot_dat'].jlist[0]

    if fixed is not None:
        fixed_reshape = np.zeros( (Nrow, Ncol) )
        fixed_reshape[0,0] = fixed[0]
        fixed_reshape[0,1] = fixed[1]
        
        for i in range(Nrow-1):
            for j in range(Ncol):
                fixed_reshape[i+1,j] = fixed[2 + j + i*Ncol ]   

    else:
        fixed_reshape = np.ones( (Nrow, Ncol) )
        fixed_reshape[0, 2] = 0
        
        
    fig, ax = plt.subplots(Nrow, Ncol, figsize=(5*Ncol, 5*Nrow))
        
    ax[0,0].hist( np.log10(res['sigma']), bins=nbin )
    ax[0,1].hist( np.log10(res['tau']), bins=nbin )
    
    
    #Plot tophat params
    for i in range(1, Nrow):
        for j in range(Ncol):
            vals = res['tophat_params'][ (i-1)*Ncol + j ]
            x2 = res['tot_dat'].jlist[i]
             
            if plot_weights:
                if j == 0: 
                    
                    bin_vals, bin_edges, im1 = ax[i,0].hist( vals, bins=nbin )
                    bins = np.zeros(len(bin_edges)-1)
                    for k in range(len(bins)):
                        bins[k] = (bin_edges[k] + bin_edges[k+1])/2
                        
                    prob_weights, _ = utils.prob_tau(x1, x2, lagvals=bins)
                    im2 = ax[i,j].fill_between( bins, np.zeros(len(bins)), prob_weights*bin_vals, 
                                            step='mid', color='r', alpha=.7 )
                    ax[i,j].plot( bins, prob_weights*bin_vals, color='k', ds='steps-mid' )

                else:
                    ax[i,j].hist( vals, bins=nbin )

            else:
                ax[i,j].hist( vals, bins=nbin )


    ax[0,0].set_xlabel(r'$\log_{10}(\sigma_{\rm DRW} \,\, ' + lc_unit_txt + ')$', fontsize=19)
    ax[0,1].set_xlabel(r'$\log_{10}(\tau_{\rm DRW} \,\, ' + time_unit_txt + ')$', fontsize=19)
        
    for i in range(1, Nrow):
        for j in range(Ncol):
            
            if j == 0:
                ax[i,j].set_xlabel(r'$t_{' + names[i-1] + '}$ ' + time_unit_txt, fontsize=22)
            if j == 1:
                ax[i,j].set_xlabel(r'w$_{' + names[i-1] + '}$ ' + time_unit_txt, fontsize=22)
            if j == 2:
                ax[i,j].set_xlabel(r's$_{' + names[i-1] + '}$ ' + lc_unit_txt, fontsize=22)  
                
    for i in range(Nrow):
        ax[i, 0].set_ylabel('N', fontsize=19)
        
    for i in range(2, Ncol):
        ax[0, i].axis('off')
        
        
    for i in range(Nrow):
        for j in range(Ncol):
            ax[i,j].tick_params('both', labelsize=14)
            ax[i,j].tick_params('both', which='major', length=8)
            ax[i,j].tick_params('both', which='minor', length=4)

            if remove_fixed:
                if fixed_reshape[i,j] == 0:
                    ax[i,j].clear()
                    ax[i,j].axis('off')

    if plot_weights:
        plt.figlegend( [im1, im2], ['Original', 'Weighted'], loc=(.81, .88), fontsize=16 )

    plt.subplots_adjust( wspace=.25, hspace=.25 )
    
    if fname is not None:
        plt.savefig( fname, dpi=200, bbox_inches='tight' )

    plt.show()
    
    return fig, ax


def javelin_corner(res, nbin=20, plot_weights=False, fname=None):
    labels = []
    labels.append( r'$\log_{10} (\sigma_{\rm DRW})$' )
    labels.append( r'$\log_{10} (\tau_{\rm DRW})$' )

    for i in range( res['tot_dat'].nlc - 1 ):
        labels.append( r'$t_{' + res['tot_dat'].names[i+1] + '}$' )
        labels.append( r'$w_{' + res['tot_dat'].names[i+1] + '}$' )
        labels.append( r'$s_{' + res['tot_dat'].names[i+1] + '}$' )
    

    #Plot original output with weighted output superposed on histograms
    corner_dat = np.vstack( [np.log10(res['sigma']), np.log10(res['tau']),  res['tophat_params']]).T
    fig = corner.corner( corner_dat,
                labels=labels, show_titles=False, bins=nbin,
                label_kwargs={'fontsize': 20} )
    
    ax = np.array(fig.axes).reshape( (2 + len(res['tophat_params']), 2 + len(res['tophat_params']) ) )
    
    if plot_weights:        
        corner_range = [[x.min(), x.max()] for x in corner_dat.T]

        #For each lag param, make another hist for the prob weighted version
        for i in range(res['tot_dat'].nlc - 1):
            vals = res['tophat_params'][3*i] 
            bin_vals, _ = np.histogram(vals, bins=nbin)
                       
            mc_vals = utils.make_mc_from_weights(res['tot_dat'].jlist[0], 
                                                 res['tot_dat'].jlist[i+1], 
                                                 bin_vals, 20)

            n_bins_1d = int(max(1, np.round(20)))
            bins_1d = np.linspace(min(corner_range[2 + 3*i]), max(corner_range[2 + 3*i]), n_bins_1d + 1)     
            
            #t hist will be on the diagonal            
            for j in range(1,  res['tot_dat'].nlc ):
                ax[3*j - 1, 3*j - 1].hist( mc_vals, bins=bins_1d, color='r', alpha=.5 )
                
            
    if fname is not None:
        plt.savefig( fname, dpi=200, bbox_inches='tight' )
        
    plt.show()
    
    return fig, ax



def plot_javelin_bestfit(res, bestfit_model, time_unit='d', lc_unit='mag', fname=None):
    
    
    cmap = ListedColormap( palettable.cartocolors.qualitative.Vivid_10.mpl_colors )    
    colors = cmap.colors

    tot_dat = res['tot_dat']

    xmin = np.min( [bestfit_model.jlist[i][0] for i in range( bestfit_model.nlc ) ] )
    xmax = np.min( [bestfit_model.jlist[i][-1] for i in range( bestfit_model.nlc ) ] )



    fig, ax = plt.subplots(tot_dat.nlc, 1, sharex=True, figsize=(13, 7))
    for i in range(len(ax)):
        ax[i].set_prop_cycle('color', palettable.cartocolors.qualitative.Bold_10.mpl_colors )


        
    if lc_unit == 'mag':
        ytxt = 'Magnitude'
        
        for i in range(tot_dat.nlc):
            ax[i].invert_yaxis()
            
    elif lc_unit == 'Arbitrary Units':
        ytxt = 'Flux'
    else:
        ytxt = 'Flux [' + str(lc_unit) + ']'
        
        
    for i in range(tot_dat.nlc):
        ax[i].errorbar( tot_dat.jlist[i], tot_dat.mlist[i] + tot_dat.blist[i], 
                    yerr=tot_dat.elist[i], fmt='.', ms=10, mfc=colors[9], mec='k',
                    label='Data', zorder=0 )
        
        ax[i].plot( bestfit_model.jlist[i], bestfit_model.mlist[i]+bestfit_model.blist[i], 
                    c=colors[i % 9], lw=2, label='Best fit', zorder=-1 )
        ax[i].fill_between( bestfit_model.jlist[i], 
                        bestfit_model.mlist[i]+bestfit_model.blist[i]+bestfit_model.elist[i],
                        bestfit_model.mlist[i]+bestfit_model.blist[i]-bestfit_model.elist[i], 
                        color=colors[i % 9], alpha=.3, zorder=-2)
        
        
        y1, y2 = ax[i].get_ylim()
        y2 = y1 + (y2 - y1)*1.25
        ax[i].set_ylim( y1, y2 )
        
        ax[i].set_xlim( xmin, xmax )
        ax[i].text( .02, .8, tot_dat.names[i], transform=ax[i].transAxes, fontsize=16 )

        ax[i].tick_params('both', labelsize=15)
        ax[i].tick_params('both', which='major', length=8)
        ax[i].tick_params('both', which='minor', length=5)
        
        ax[i].legend( bbox_to_anchor=(1.13, 1.03), fontsize=12 )

            
    ax[-1].set_xlabel('Time [' + str(time_unit) + ']', fontsize=20)


    plt.figtext(.05, .5, 'Flux', va='center', rotation=90, fontsize=22)

    plt.subplots_adjust(hspace=0)
    
    if fname is not None:
        plt.savefig( fname, dpi=200, bbox_inches='tight' )
        
    plt.show()
    
    return fig, ax



###################################################################
#########################  DEPRECATED #############################
###################################################################

def plot_pyzdcf_results_peakcent(x1, y1, yerr1, x2, y2, yerr2, 
                       dcf_lags, dcf,
                       dcf_lag_err, dcf_err,
                       cccd_lags, ccpd_lags,
                       time_unit='d', flux_unit='mJy', lc_names=['', ''],
                       fname=None, show=False):

    cent_med = np.median( cccd_lags )
    cent_hi = np.percentile( cccd_lags, 84 )
    cent_lo = np.percentile( cccd_lags, 16 )

    peak_med = np.median( ccpd_lags )
    peak_hi = np.percentile( ccpd_lags, 84 )
    peak_lo = np.percentile( ccpd_lags, 16 )
    
    min_lag = np.min(dcf_lags)
    max_lag = np.max(dcf_lags)
    
    
    
    if (cent_med < .5) & (peak_med < .5) & (time_unit == 'd'):
        cent_med *= 24
        cent_hi *= 24
        cent_lo *= 24
        
        peak_med *= 24
        peak_hi *= 24
        peak_lo *= 24
        
        lag_unit = 'hr'
    else:
        lag_unit = time_unit
    

    
    #Setup plot axes    
    fig = plt.figure( figsize=(12, 7) )
    gs = gridspec.GridSpec(ncols=3, nrows=6)

    ax1 = fig.add_subplot(gs[:3, :2])
    ax2 = fig.add_subplot(gs[3:6, :2], sharex=ax1)

    ax3 = fig.add_subplot(gs[:2, 2])
    ax4 = fig.add_subplot(gs[2:4, 2], sharex=ax3)
    ax5 = fig.add_subplot(gs[4:6, 2], sharex=ax3, sharey=ax4)



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
    #Plot DCF
    _, _, bars = ax3.errorbar( dcf_lags, dcf, yerr=dcf_err, xerr=dcf_lag_err, fmt='.k' )
    [bar.set_alpha(.2) for bar in bars]
    
    ax3.set_ylim(-1, 1)

    #----------------------------
    #Plot CCCD
    _, bins, _ = ax4.hist( cccd_lags, 50, range=(min_lag, max_lag), density=True )

    #Increase y-bounds to fit text
    ymin, ymax = ax4.get_ylim()
    ymax = ymin + (ymax - ymin)*1.1
    ax4.set_ylim(ymin, ymax)

    ax4.text(.05, .8, 
            '${:.2f} ^'.format(cent_med) + '{+' + '{:.2f}'.format(cent_hi-cent_med)  + '}_{-' + '{:.2f}'.format(cent_med-cent_lo) + '}$ ' + lag_unit,
            transform=ax4.transAxes,
            fontsize=14)

    #----------------------------
    #Plot CCPD
    ax5.hist( ccpd_lags, bins=bins, range=(min_lag, max_lag), density=True )

    #Increase y-bounds to fit text
    ymin, ymax = ax5.get_ylim()
    ymax = ymin + (ymax - ymin)*1.1
    ax5.set_ylim(ymin, ymax)

    ax5.text(.05, .8, 
            '${:.2f} ^'.format(peak_med) + '{+' + '{:.2f}'.format(peak_hi-peak_med)  + '}_{-' + '{:.2f}'.format(peak_med-peak_lo) + '}$ ' + lag_unit,
            transform=ax5.transAxes,
            fontsize=14)
    ax5.set_xlabel('Lag [' + str(lag_unit) + ']', fontsize=19)




    plt.figtext( .05, .5, 'Flux [' + str(flux_unit) + ']', fontsize=19, rotation=90, va='center' )

    plt.figtext( .94, .72, 'DCF', rotation=270, fontsize=19 )
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
            ax.tick_params('both', labelsize=12)
        else:
            ax.tick_params('x', labelsize=0)
            ax.tick_params('y', labelsize=12)

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
