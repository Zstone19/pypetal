import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, BboxConnector, TransformedBbox
from matplotlib import gridspec

import numpy as np
from scipy.signal import peak_widths

from pypetal import pyccf
from pypetal.petalio import err2str, write_data
from pypetal import defaults



def find_overlap(x1, x2, gaps):
        
    sort_ind = np.argsort(x1)
    x1 = x1[sort_ind]
    
    sort_ind = np.argsort(x2)
    x2 = x2[sort_ind]
        
    low_lim = np.maximum( x1[0], x2[0] )
    up_lim = np.minimum( x1[-1], x2[-1] )

    Ntot = len( np.argwhere( (x2 >= low_lim) & (x2 <= up_lim) ).T[0] )

    #Remove points in x2 that fall into gaps in x1
    if gaps is not None:
        for i in range(len(gaps[0])):
            lgap = gaps[0][i]
            rgap = gaps[1][i]
                        
            bad_ind = np.argwhere( (x2 > lgap) & (x2 < rgap) ).T[0]
            Ntot -= len(bad_ind)
        
    return Ntot  



def prob_tau(x1, x2, laglim=None, lagvals=None, Nlag=1000, gap_size=30, k=2):
        
    sort_ind = np.argsort(x1)
    x1 = x1[sort_ind]
    
    sort_ind = np.argsort(x2)
    x2 = x2[sort_ind]
    
    baseline = np.max([ x1[-1]-x1[0], x2[-1]-x2[0] ])
    Nlag = int(Nlag)
    
    if (laglim is None) & (lagvals is None):
        laglim = [-baseline, baseline]
        lags = np.linspace(laglim[0], laglim[1], Nlag)
    
    elif (lagvals is None) & (laglim is not None):
        lags = np.linspace(laglim[0], laglim[1], Nlag)
    
    else:
        lags = lagvals
    
    
    #Get where gaps are in data
    gap_ind = np.argwhere(np.diff(x1) > gap_size).T[0]    
    gaps_left = []
    gaps_right = []
    
    for i in range(len(gap_ind)):
        gaps_left.append( x1[gap_ind[i]] )
        gaps_right.append( x1[gap_ind[i]+1] ) 
        
    gaps = np.array([ gaps_left, gaps_right ])
            
    if len(gaps_left) == 0:
        gaps = None
        
    nvals = np.zeros_like(lags)
    N0 = find_overlap( x1, x2, gaps)
    
    for i, tau in enumerate(lags):                
        nvals[i] = find_overlap( x1, x2-tau, gaps )
                    
    probs = (nvals / N0)**k
    
    return probs, nvals, N0, lags



def get_acf(x, y, interp=None, lag_bounds=None, 
            sigmode=0.2, thres=0.8):

    sort_ind = np.argsort(x)
    x = x[sort_ind]
    y = y[sort_ind]

    mean_diff = np.mean(np.diff(x))
    baseline = x[-1] - x[0]

    if lag_bounds is None:
        lag_bounds = [-baseline, baseline]
        
    if interp is None:
        interp = mean_diff/2

    _, _, _, _, ccf_pack, \
        _, _, _ = pyccf.peakcent(x, y, x, y, lag_bounds[0], lag_bounds[1], interp,
                                 thres=thres, sigmode=sigmode)
        
    return ccf_pack[0], ccf_pack[1]





def get_weights(x1, y1, x2, y2, interp=None, lag_bounds=None,
                sigmode=0.2, thres=0.8, gap_size=30, k=2):

    #Get ACF for continuum light curve    
    acf, acf_lags = get_acf(x1, y1, interp=interp, lag_bounds=lag_bounds,
                            sigmode=sigmode, thres=thres)
    

    #Set all values past first negative to zero
    ind1 = np.max( np.argwhere( (acf < 0) & (acf_lags < 0) ).T[0] )
    ind2 = np.min( np.argwhere( (acf < 0) & (acf_lags > 0) ).T[0] )
    acf_new = acf.copy()
    acf_new[:ind1+1] = 0
    acf_new[ind2:] = 0
    
    #Get p(tau)
    probs, ntau, n0, lags = prob_tau(x1, x2, lagvals=acf_lags, gap_size=gap_size, k=k)
    
    #Convolve ACF with p(tau)
    tot_weight = np.convolve(acf_new, probs, mode='same')
    prob_dist = tot_weight / np.max(tot_weight)
    
    return prob_dist, lags, ntau, acf, n0




def gaussian(x, sigma):
    c = np.sqrt(2 * np.pi) * sigma
    return np.exp(-0.5 * (x / sigma)**2) / c

def get_bounds(dist, weights, lags, width=15):
    
    #Bin dist into weight lags
    dbin = lags[1] - lags[0]
    bin0 = np.min(lags)

    bin_edges = []
    bin0 = np.min(lags) - dbin/2

    for i in range(len(lags)+1):
        bin_edges.append( bin0 + i*dbin )

    hist, _ = np.histogram(dist, bins=bin_edges)
    hist = np.array(hist, dtype=float)    
    
    #Weight histogram
    weighted_hist = hist*weights
    
    #Convolve weighted distribution with gaussian
    gauss = gaussian( lags, width )
    smooth_dist = np.convolve( hist, gauss, mode='same')
    smooth_weight_dist = np.convolve( weighted_hist, gauss, mode='same')

    #Find peak
    peak_ind = np.argmax(smooth_weight_dist)
    
    #Find peak bounds
    res = peak_widths( smooth_weight_dist, [peak_ind], rel_height=.99 ) 
    
    peak = lags[peak_ind]
    bound_left = lags[ np.floor(res[2]).astype(int) ] + dbin*( res[2]%1 )
    bound_right = lags[ np.floor(res[3]).astype(int) ] + dbin*( res[3]%1 )

    return bound_left[0], peak, bound_right[0], smooth_dist, smooth_weight_dist
    
    
    


def run_weighting( cont_fname, line_fnames, output_dir, line_names, 
                  run_pyccf, run_javelin, 
                  pyccf_res, javelin_res, 
                  pyccf_params, javelin_params,
                  general_kwargs, kwargs):


    #General kwargs
    verbose = general_kwargs['verbose']
    plot = general_kwargs['plot']
    time_unit = general_kwargs['time_unit']
    lc_unit = general_kwargs['lc_unit']   
    lag_bounds = general_kwargs['lag_bounds']

    #---------------------------
    # pyCCF kwargs
    interp, nsim, mcmode, sigmode, thres, nbin = defaults.set_pyccf( pyccf_params, 
                                                                    np.hstack([ [cont_fname], line_fnames ])
                                                                    )

    #---------------------------
    #JAVELIN kwargs
    _, _, _, _, _, _, _, _, _, _, _, _, _, together, _ = defaults.set_javelin(javelin_params, 
                                                                              np.hstack([ [cont_fname], line_fnames ])
                                                                              )  

    #---------------------------
    #Weighting kwargs    
    gap_size, k, width, zoom = defaults.set_weighting(kwargs)

    #---------------------------
    
    #Make sure number of LCs is the same
    if run_pyccf:
        nlc1 = len(pyccf_res)
        
        if nlc1 != len(line_fnames):
            raise Exception('Number of line LCs does not match number of pyCCF results')
                
        nlc = len(line_fnames)
        
    if run_javelin:
        if together:
            nlc2 = javelin_res['tot_dat'].nlc - 1
        else:
            nlc2 = len(javelin_res)

        if nlc2 != len(line_fnames):
            raise Exception('Number of line LCs does not match number of JAVELIN results')
        
    if run_pyccf & run_javelin:
        if nlc1 != nlc2:
            raise Exception('ERROR: Number of light curves in pyCCF and JAVELIN results do not match.')

    if (not run_pyccf) & (not run_javelin):
        return {}



    nlc = len(line_fnames)
    x_cont, y_cont, yerr_cont = np.loadtxt(cont_fname, usecols=[0,1,2], unpack=True, delimiter=',')
    
    output = {
        'pyccf': {},
        'javelin': {},
        'rmax': []
    }
    
    #Run weighting on pyCCF CCCDs
    if run_pyccf:
        output['pyccf']['centroid'] = []
        output['pyccf']['bounds'] = []
        output['pyccf']['acf'] = []
        output['pyccf']['lags'] = []
        output['pyccf']['weight_dist'] = []
        output['pyccf']['smoothed_dist'] = []
        output['pyccf']['ntau'] = []
        output['pyccf']['CCCD'] = []
        output['pyccf']['downsampled_CCCD'] = []
        output['pyccf']['frac_rejected'] = []
        
        
        for i in range(nlc):
            cccd_lags = pyccf_res[i]['CCCD_lags']
            x_line, y_line, yerr_line = np.loadtxt(line_fnames[i], usecols=[0,1,2], unpack=True, delimiter=',')
            
            prob_dist, lags, ntau, acf, n0 = get_weights(x_cont, y_cont, x_line, y_line, 
                                                interp=interp, lag_bounds=lag_bounds[i],
                                                sigmode=sigmode, thres=thres, 
                                                gap_size=gap_size, k=k)
            
            min_bound, peak, max_bound, smooth_dist, smooth_weight_dist = get_bounds(cccd_lags, prob_dist, lags, width=width)
            downsampled_cccd = cccd_lags[(cccd_lags > min_bound) & (cccd_lags < max_bound)]
            
            med_cent = np.median(downsampled_cccd)
            cent_err_lo = med_cent - np.percentile( downsampled_cccd, 16 )
            cent_err_hi = np.percentile( downsampled_cccd, 84 ) - med_cent
            
            
            #Write diagnostic info
            write_data( [ lags, ntau, prob_dist, acf, smooth_dist, smooth_weight_dist ], 
                        output_dir + line_names[i+1] + '/weights/pyccf_weights.dat',
                        '#lags,ntau,weight_dist,acf,smooth_dist,smooth_weight_dist')
            
            write_data( downsampled_cccd, 
                       output_dir + line_names[i+1] + '/weights/pyccf_weighted_cccd.dat')
            

            #Output results            
            output['pyccf']['centroid'].append( [cent_err_lo, med_cent, cent_err_hi] )
            output['pyccf']['bounds'].append( [min_bound, peak, max_bound] )
            output['pyccf']['acf'].append( acf )
            output['pyccf']['lags'].append( lags )
            output['pyccf']['weight_dist'].append( prob_dist )
            output['pyccf']['smoothed_dist'].append( smooth_weight_dist )
            output['pyccf']['ntau'].append( ntau )
            output['pyccf']['CCCD'].append( cccd_lags )
            output['pyccf']['downsampled_CCCD'].append( downsampled_cccd )
            output['pyccf']['frac_rejected'].append( 1 - len(downsampled_cccd) / len(cccd_lags) )
            
            
            #Plot weights
            plot_weights( output_dir, line_names[i+1], output['pyccf'], n0, k, 
                         time_unit=time_unit, plot=plot )
            
            
            
    #Run weighting on JAVELIN lag distributions
    if run_javelin:
        output['javelin']['tophat_lag'] = []
        output['javelin']['bounds'] = []
        output['javelin']['acf'] = []
        output['javelin']['lags'] = []
        output['javelin']['weight_dist'] = []
        output['javelin']['smoothed_dist'] = []
        output['javelin']['ntau'] = []
        output['javelin']['lag_dist'] = []
        output['javelin']['downsampled_lag_dist'] = []
        output['javelin']['frac_rejected'] = []

        for i in range(nlc):
            
            if together:
                lag_dist = javelin_res['tophat_params'][3*i]
            else:
                lag_dist = javelin_res[i]['tophat_params'][0]

            x_line, y_line, yerr_line = np.loadtxt(line_fnames[i], usecols=[0,1,2], unpack=True, delimiter=',')
            
            prob_dist, lags, ntau, acf, n0 = get_weights(x_cont, y_cont, x_line, y_line, 
                                                interp=interp, lag_bounds=lag_bounds[i],
                                                sigmode=sigmode, thres=thres, 
                                                gap_size=gap_size, k=k)

            
            min_bound, peak, max_bound, smooth_dist, smooth_weight_dist = get_bounds(lag_dist, prob_dist, lags, width=width)
            downsampled_dist = lag_dist[(lag_dist > min_bound) & (lag_dist < max_bound)]
            
            med_lag = np.median(downsampled_dist)
            lag_err_lo = med_lag - np.percentile( downsampled_dist, 16 )
            lag_err_hi = np.percentile( downsampled_dist, 84 ) - med_lag
            
            
            #Write diagnostic info
            if not together:
                write_data( [ lags, ntau, prob_dist, acf, smooth_dist, smooth_weight_dist ], 
                            output_dir + line_names[i+1] + '/weights/javelin_weights.dat',
                            '#lags,ntau,weight_dist,acf,smooth_dist,smooth_weight_dist')
                
                write_data( downsampled_dist, 
                        output_dir + line_names[i+1] + '/weights/javelin_weighted_lag_dist.dat')
            
            else:
                write_data( [ lags, ntau, prob_dist, acf, smooth_dist, smooth_weight_dist ], 
                        output_dir + 'javelin/' + line_names[i+1] + '_javelin_weights.dat',
                        '#lags,ntau,weight_dist,acf,smooth_dist,smooth_weight_dist')
            
                write_data( downsampled_dist, 
                        output_dir + 'javelin/' + line_names[i+1] + '_javelin_weighted_lag_dist.dat')           
        
        
            #Output results
            output['javelin']['tophat_lag'].append( [lag_err_lo, med_lag, lag_err_hi] )
            output['javelin']['bounds'].append( [min_bound, peak, max_bound] )
            output['javelin']['acf'].append( acf )
            output['javelin']['lags'].append( lags )
            output['javelin']['weight_dist'].append( prob_dist )
            output['javelin']['smoothed_dist'].append( smooth_weight_dist )
            output['javelin']['ntau'].append( ntau )
            output['javelin']['lag_dist'].append( lag_dist )
            output['javelin']['downsampled_lag_dist'].append( downsampled_dist )
            output['javelin']['frac_rejected'].append( 1 - len(downsampled_dist) / len(lag_dist) )
        
        
            #Plot weights
            if not run_pyccf:
                plot_weights( output_dir, line_names[i+1], output['javelin'], n0, k, 
                             time_unit=time_unit, plot=plot )

            
            
    if run_pyccf & run_javelin:
        
        for i in range(nlc):
            lag = output['javelin']['tophat_lag'][i][1]
            lag_err_hi = output['javelin']['tophat_lag'][i][2]
            lag_err_lo = output['javelin']['tophat_lag'][i][0]
            
            
            ccf = pyccf_res[i]['CCF']
            ccf_lags = pyccf_res[i]['CCF_lags']
            
            good_ind = np.argwhere( ( ccf_lags >= lag-lag_err_lo ) | ( ccf_lags <= lag+lag_err_hi ) )
            rmax = np.max(ccf[good_ind])
            
            output['rmax'].append(rmax)
            
        write_data( output_dir + 'rmax.dat', output['rmax'] )
    
            
    #Plot results for pyCCF
    if run_pyccf:
        plot_weight_output( output_dir, cont_fname, line_fnames, line_names,
                     output['pyccf'], general_kwargs, module='pyccf', 
                     zoom=zoom, plot=plot)
    
    
    #Plot results for JAVELIN
    if run_javelin:
        plot_weight_output( output_dir, cont_fname, line_fnames, line_names,
                    output['javelin'], general_kwargs, module='javelin', 
                    zoom=zoom, plot=plot)
    
    
    return output


def plot_weights(output_dir, line_name, res, n0, k, time_unit='d', plot=False):
    lags = res['lags'][-1]
    acf = res['acf'][-1]
    prob_dist = res['weight_dist'][-1]
    ntau = res['ntau'][-1]    
    
    #Get P(tau)
    probs = (ntau/n0)**k
    
    #Set all ACF values past first negative to zero
    ind1 = np.max( np.argwhere( (acf < 0) & (lags < 0) ).T[0] )
    ind2 = np.min( np.argwhere( (acf < 0) & (lags > 0) ).T[0] )
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
    verbose = general_kwargs['verbose']
    plot = general_kwargs['plot']
    time_unit = general_kwargs['time_unit']
    


    nlc = len(res['lags'])
        
    Ncol = 3
    Nrow = nlc//Ncol + 1
    
    #--------------------------------------------------------------------------------
        
    if nlc <= 3:
        Nrow = 1
        Ncol = nlc
        
    gs = gridspec.GridSpec(Nrow, Ncol)
    ax_tot = np.zeros( (Nrow, Ncol), dtype=object )
    
    fig, ax = plt.subplots(Nrow, Ncol, figsize=( 6*Ncol, 5*Nrow ))
    for i in range(Nrow):
        for j in range(Ncol):
            sub_gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[i, j], 
                                                    height_ratios=[1, 3], hspace=0)
    
            ax_bot = plt.subplot(sub_gs[1])
            ax_top = plt.subplot(sub_gs[0], sharex=ax_bot)
            
            ax_tot[i, j] = [ax_top, ax_bot]
            
            
    for i in range(nlc):
        
        if module == 'pyccf':
            lags = res['lags'][i]
            lag_dist = res['CCCD'][i]

            llim = res['bounds'][i][0]
            rlim = res['bounds'][i][2]
            peak = res['bounds'][i][1]
            
            weight_dist = res['weight_dist'][i]
            smooth_dist = res['smoothed_dist'][i]
            acf = res['acf'][i]
            
            lag_err_lo = res['centroid'][i][0]
            lag_err_hi = res['centroid'][i][2]
            lag_value = res['centroid'][i][1]
            
            xlabel = r'{ \rm CCCD }'
            
        elif module == 'javelin':
            lags = res['lags'][i]
            lag_dist = res['lag_dist'][i]
            
            llim = res['bounds'][i][0]
            rlim = res['bounds'][i][2]
            peak = res['bounds'][i][1]
            
            weight_dist = res['weight_dist'][i]
            smooth_dist = res['smoothed_dist'][i]
            acf = res['acf'][i]
            
            lag_err_lo = res['tophat_lag'][i][0]
            lag_err_hi = res['tophat_lag'][i][2]
            lag_value = res['tophat_lag'][i][1]
            
            xlabel = r't'
        
        
        
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
        ind1 = np.max( np.argwhere( (acf < 0) & (lags < 0) ).T[0] )
        ind2 = np.min( np.argwhere( (acf < 0) & (lags > 0) ).T[0] )
        acf[:ind1] = 0
        acf[ind2:] = 0
        
        im3, = ax_tot[row_ind, col_ind][0].plot(lags, acf, c='r')                
    
    
        #Plot weighting function
        im1, = ax_tot[row_ind, col_ind][0].plot(lags, weight_dist, c='k')
        ax_tot[row_ind, col_ind][0].set_yticks([.5, 1])
                    
        #Plot smoothed distribution
        im2, = ax_tot[row_ind, col_ind][0].plot(lags, smooth_dist/np.max(smooth_dist), c='DodgerBlue')
        ax_tot[row_ind, col_ind][0].axvspan( llim, rlim, color='k', alpha=.1 )
                    
        #Write lag and error
        peak_str = err2str( lag_value, lag_err_lo, lag_err_hi, dec=2 )
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
    plt.figlegend( [im1, im2, im3], [r'w($\tau$)', 'Smoothed Dist', 'ACF'],
                    bbox_to_anchor=(1,.9), fontsize=11)
    plt.figlegend( [bin1, bin2], ['Original', 'Weighted'],
                    bbox_to_anchor=(.98, .7), fontsize=11)       
     
    plt.savefig( output_dir + module + '_weights_res.pdf', dpi=200, bbox_inches='tight' )
    
    
    if plot:
        plt.show()
        
    plt.cla()
    plt.clf()
    plt.close()
        
    return   
    