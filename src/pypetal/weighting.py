import numpy as np
import os
import glob
import warnings

import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import (BboxConnector,
                                                   TransformedBbox, inset_axes)
from scipy.signal import peak_widths
import pickle

from pypetal import defaults
from pypetal.petalio import err2str, write_data, write_weighting_summary
from pypetal import pyccf
from pypetal.load import get_ordered_line_names
from pypetal.pyroa_funcs import get_samples_chunks


#############################################
############## UTILS ########################
#############################################



def find_overlap(x1, x2, gaps):

    """Find the amount of overlapping data points between two light curves,
    taking into account seasonal gaps in the data.

    Parameters
    ----------

    x1 : list of float
        The time values of the first light curve.

    x2 : list of float
        The time values of the second light curve.

    gaps : list of float
        A list of gaps in the data for the first light curve. The list should be structured
        a a list of 2 element arrays, where the first element is the start of the gap and the
        second element is the end of the gap.




    Returns
    -------

    Ntot : int
        The number of overlapping data points between the two light curves.

    """

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


    """Calculate the probability distribution P(tau) to weight the time lag distribution described in
    Grier et al. (2019).



    Parameters
    ----------

    x1 : list of float
        The time values of the first light curve.

    x2 : list of float
        The time values of the second light curve.

    laglim : list of float, optional
        The limits of lags to search over. If ``None``, the limits will be set to the baseline of the
        light curves. If not ``None``, Nlag must be set. Default is ``None``.

    lagvals : list of float, optional
        The values of lags to search over. If ``None``, the limits will be set to the baseline of the
        light curves. If not ``None``, Nlag must be set. Default is ``None``.

    Nlag : int, optional
        The number of lags to search over. Default is 1000.

    gap_size : float, optional
        The minimum gap size to use when searching for seasonal gaps in the first light curve.
        Default is 30.

    k : int, optional
        The power to raise the probability distribution to. Default is 2.




    Returns
    -------

    probs : list of float
        The probability distribution P(tau).

    nvals : list of float
        The number of overlapping data points for each lag N(tau).

    N0 : int
        The number of overlapping points when tau = 0.

    lags : list of float
        The lags used to calculate the probability distribution.

    """

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


    """Calculate the autocorrelation function of a light curve.


    Parameters
    ----------

    x : list of float
        The time values of the light curve.

    y : list of float
        The flux values of the light curve.

    interp : float, optional
        The interpolation distance to use for pyCCF. Default is ``None``.

    lag_bounds : list of float, optional
        The limits of lags to search over. If ``None``, the limits will be set to the baseline of the
        light curve. Default is ``None``.

    sigmode : float, optional
        The ``sigmode`` option for pyCCF. Default is 0.2.

    thres : float, optional
        The ``thres`` option for pyCCF. Default is 0.8.


    Returns
    -------

    acf : list of float
        The autocorrelation function of the light curve.

    lags : list of float
        The lags used to calculate the ACF.

    """


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


    """Calculate the total weighting distribution w(tau) described in Grier et al. (2019).
    This calculates the probability weights depending on the number of overlapping data points as
    a function of time lag P(tau) and convolves it with the ACF of the continuum light curve.


    Parameters
    ----------

    x1 : list of float
        The time values of the first light curve.

    y1 : list of float
        The flux values of the first light curve.

    x2 : list of float
        The time values of the second light curve.

    y2 : list of float
        The flux values of the second light curve.

    interp : float, optional
        The interpolation distance to use for pyCCF. Default is ``None``.

    lag_bounds : list of float, optional
        The limits of lags to search over. If ``None``, the limits will be set to the baseline of the
        light curve. Default is ``None``.

    sigmode : float, optional
        The ``sigmode`` option for pyCCF. Default is 0.2.

    thres : float, optional
        The ``thres`` option for pyCCF. Default is 0.8.

    gap_size : float, optional
        The minimum gap size to use when searching for seasonal gaps in the first light curve.
        Default is 30.

    k : int, optional
        The power to raise the probability distribution to. Default is 2.



    Returns
    -------

    prob_dist : list of float
        The total weighting distribution w(tau).

    lags : list of float
        The lags used to calculate the total weighting distribution.

    ntau : list of float
        The number of overlapping data points as a function of time lag N(tau).

    acf : list of float
        The ACF of the first light curve.

    n0 : float
        The number of overlapping data points at a lag tau = 0.

    """

    #Get ACF for continuum light curve
    acf, acf_lags = get_acf(x1, y1, interp=interp, lag_bounds=lag_bounds,
                            sigmode=sigmode, thres=thres)


    #Set all values past first negative to zero

        #Left side
    left_inds = np.argwhere( (acf < 0) & (acf_lags < 0) ).T[0]
    if len(left_inds) == 0:
        ind1 = 0
    else:
        ind1 = np.max( left_inds )

        #Right side
    right_inds = np.argwhere( (acf < 0) & (acf_lags > 0) ).T[0]
    if len(right_inds) == 0:
        ind2 = len(acf)-1
    else:
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


    """Find the bounds of the peak of the weighted distribution, given the original distribution and weights.
    This will convolve this weighted distribution with a Gaussian and find the bounds of the primary (i.e tallest)
    peak in the distribution.


    Parameters
    ----------

    dist : list of float
        The original distribution.

    weights : list of float
        The weights to apply to the distribution.

    lags : list of float
        The lags used to calculate the distribution.

    width : float, optional
        The width of the Gaussian to convolve the distribution with. Default is 15.



    Returns
    -------

    peak_lower_bound : float
        The lag of the lower bound of the peak of the distribution.

    peak: int
        The lag of the peak of the distribution.

    peak_upper_bound : float
        The lag of the upper bound of the peak of the distribution.

    smooth_dist : list of float
        The smoothed original distribution.

    smooth_weight_dist : list of float
        The smoothed weighted distribution.

    """

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



#############################################
############## RUNNING ########################
#############################################


def combine_weight_outputs(res_arr, run_pyccf, run_javelin, run_pyroa):

    output = {
        'pyccf': {},
        'javelin': {},
        'pyroa': {},
        'rmax_pyccf': [],
        'rmax_javelin': [],
        'rmax_pyroa': []
    }

    cols_pyccf = ['centroid']
    cols_javelin = ['tophat_lag']
    cols_pyroa = ['time_delay']
    
    cols_common = ['bounds', 'acf', 'lags', 'weight_dist', 'smoothed_dist', 'ntau', 'frac_rejected']
    for col in cols_common:
        cols_pyccf.append(col)
        cols_javelin.append(col)
        cols_pyroa.append(col)
        
    cols_pyccf.append('CCCD')
    cols_javelin.append('lag_dist')
    cols_pyroa.append('lag_dist')
    

    run_arr = [run_pyccf, run_javelin, run_pyroa]
    cols_arr = [cols_pyccf, cols_javelin, cols_pyroa]
    names_arr = ['pyccf', 'javelin', 'pyroa']

    for run, cols, name in zip(run_arr, cols_arr, names_arr):
        if run:
            for col in cols:
                output[name][col] = []
                

    for res in res_arr:
        for run, cols, name in zip(run_arr, cols_arr, names_arr):
            if run:
                for col in cols:
                    output[name][col].append( res[name][col] )


        #rmax
        if run_pyccf:
            output['rmax_pyccf'].append( res['rmax_pyccf'] )
            
            if run_javelin:
                output['rmax_javelin'].append( res['rmax_javelin'] )
            if run_pyroa:
                output['rmax_pyroa'].append( res['rmax_pyroa'] )     
    
    return output






def run_weighting_single( output_dir, cont_fname, line_fname,
                          weighting_kwargs, lag_bounds='baseline',
                          jav_chain_file=None, 
                          pyccf_iccf_file=None, pyccf_dist_file=None,
                          pyroa_sample_file=None,
                          javelin_lag_col=2, pyroa_obj_ind=1, 
                          pyccf_params={}, pyroa_params = {},
                          plot=False, time_unit='d'):

    #---------------------------
    # pyCCF kwargs
    interp, _, _, sigmode, thres, _ = defaults.set_pyccf(pyccf_params)

    #---------------------------
    #Weighting kwargs
    gap_size, k, width, _ = defaults.set_weighting(weighting_kwargs)

    #---------------------------
    #See what modules to run

    if jav_chain_file is not None:
        run_javelin = True
    else:
        run_javelin = False

    if pyccf_dist_file is not None:
        run_pyccf = True
    else:
        run_pyccf = False
        
    if pyroa_sample_file is not None:
        run_pyroa = True
    else:
        run_pyroa = False


    cols = ['n0', 'peak_bounds', 'peak', 'lag', 'lag_err', 'frac_rejected']
    modules = ['pyccf', 'javelin', 'pyroa']

    summary_dict = {}
    summary_dict['k'] = k

    for mod in modules:
        for col in cols:
            summary_dict[ col + '_' + mod ] = None

        summary_dict['rmax_' + mod] = None


    output = {}
    for mod in modules:
        output[mod] = {}
        output['rmax_' + mod] = []


    #---------------------------
    #Get continuum data
    x_cont, y_cont, _ = np.loadtxt(cont_fname, usecols=[0,1,2], unpack=True, delimiter=',')
    x_line, y_line, _ = np.loadtxt(line_fname, unpack=True, delimiter=',', usecols=[0,1,2])

    if lag_bounds == 'baseline':
        baseline = np.max([ x_cont.max(), x_line.max() ]) - np.min([ x_cont.min(), x_line.min() ])
        lag_bounds = [-baseline, baseline]

    #---------------------------
    #pyCCF

    if run_pyccf:
        cccd_lags, ccpd_lags = np.loadtxt(pyccf_dist_file, unpack=True, delimiter=',')
        prob_dist, lags, ntau, acf, n0 = get_weights(x_cont, y_cont, x_line, y_line,
                                    interp=interp, lag_bounds=lag_bounds,
                                    sigmode=sigmode, thres=thres,
                                    gap_size=gap_size, k=k)

        min_bound, peak, max_bound, smooth_dist, smooth_weight_dist = get_bounds(cccd_lags, prob_dist, lags, width=width)
        downsampled_cccd = cccd_lags[(cccd_lags > min_bound) & (cccd_lags < max_bound)]

        med_cent = np.median(downsampled_cccd)
        cent_err_lo = med_cent - np.percentile( downsampled_cccd, 16 )
        cent_err_hi = np.percentile( downsampled_cccd, 84 ) - med_cent 



        #Write diagnostic info
        write_data( [ lags, ntau, prob_dist, acf, smooth_dist, smooth_weight_dist ],
                    output_dir + 'pyccf_weights.dat',
                    '#lags,ntau,weight_dist,acf,smooth_dist,smooth_weight_dist')

        write_data( downsampled_cccd,
                    output_dir + 'pyccf_weighted_cccd.dat')


        #Add to summary dict
        summary_dict['n0_pyccf'] = n0
        summary_dict['peak_bounds_pyccf'] = [min_bound, max_bound]
        summary_dict['peak_pyccf'] = peak
        summary_dict['lag_pyccf'] = med_cent
        summary_dict['lag_err_pyccf'] = [cent_err_lo, cent_err_hi]
        summary_dict['frac_rejected_pyccf'] = 1 - len(downsampled_cccd) / len(cccd_lags)


        #Add to weighting results
        output['pyccf']['centroid'] = [cent_err_lo, med_cent, cent_err_hi]
        output['pyccf']['bounds'] = [min_bound, peak, max_bound]
        output['pyccf']['acf'] = acf
        output['pyccf']['lags'] = lags
        output['pyccf']['weight_dist'] = prob_dist
        output['pyccf']['smoothed_dist'] = smooth_weight_dist
        output['pyccf']['ntau'] = ntau
        output['pyccf']['CCCD'] = cccd_lags
        output['pyccf']['downsampled_CCCD'] = downsampled_cccd
        output['pyccf']['frac_rejected'] = 1 - len(downsampled_cccd) / len(cccd_lags)
        
    else:
        #Add to summary dict
        summary_dict['n0_pyccf'] = np.nan
        summary_dict['peak_bounds_pyccf'] = [np.nan, np.nan]
        summary_dict['peak_pyccf'] = np.nan
        summary_dict['lag_pyccf'] = np.nan
        summary_dict['lag_err_pyccf'] = [np.nan, np.nan]
        summary_dict['frac_rejected_pyccf'] = np.nan


        #Add to weighting results
        output['pyccf']['centroid'] =  [np.nan, np.nan, np.nan]
        output['pyccf']['bounds'] =  [np.nan, np.nan, np.nan]
        output['pyccf']['acf'] = np.nan
        output['pyccf']['lags'] = np.nan
        output['pyccf']['weight_dist'] = np.nan
        output['pyccf']['smoothed_dist'] = np.nan
        output['pyccf']['ntau'] = np.nan
        output['pyccf']['CCCD'] = np.nan
        output['pyccf']['downsampled_CCCD'] = np.nan
        output['pyccf']['frac_rejected'] = np.nan


    #---------------------------
    #JAVELIN

    if run_javelin:
        jav_chains = np.loadtxt(jav_chain_file, unpack=True)
        lag_dist = jav_chains[javelin_lag_col]
        prob_dist, lags, ntau, acf, n0 = get_weights(x_cont, y_cont, x_line, y_line,
                                            interp=interp, lag_bounds=lag_bounds,
                                            sigmode=sigmode, thres=thres,
                                            gap_size=gap_size, k=k)

        min_bound, peak, max_bound, smooth_dist, smooth_weight_dist = get_bounds(lag_dist, prob_dist, lags, width=width)
        downsampled_dist = lag_dist[(lag_dist > min_bound) & (lag_dist < max_bound)]

        med_lag = np.median(downsampled_dist)
        lag_err_lo = med_lag - np.percentile( downsampled_dist, 16 )
        lag_err_hi = np.percentile( downsampled_dist, 84 ) - med_lag


        #Write diagnostic info
        write_data( [ lags, ntau, prob_dist, acf, smooth_dist, smooth_weight_dist ],
                    output_dir + 'javelin_weights.dat',
                    '#lags,ntau,weight_dist,acf,smooth_dist,smooth_weight_dist')

        write_data( downsampled_dist,
                output_dir + 'javelin_weighted_lag_dist.dat')



        #Add to summary dict
        summary_dict['n0_javelin'] = n0
        summary_dict['peak_bounds_javelin'] = [min_bound, max_bound]
        summary_dict['peak_javelin'] = peak
        summary_dict['lag_javelin'] = med_lag
        summary_dict['lag_err_javelin'] = [lag_err_lo, lag_err_hi]
        summary_dict['frac_rejected_javelin'] = 1 - len(downsampled_dist) / len(lag_dist)


        #Add to weighting results
        output['javelin']['tophat_lag'] =  [lag_err_lo, med_lag, lag_err_hi]
        output['javelin']['bounds'] =  [min_bound, peak, max_bound]
        output['javelin']['acf'] = acf
        output['javelin']['lags'] = lags
        output['javelin']['weight_dist'] = prob_dist
        output['javelin']['smoothed_dist'] = smooth_weight_dist
        output['javelin']['ntau'] = ntau
        output['javelin']['lag_dist'] = lag_dist
        output['javelin']['downsampled_lag_dist'] = downsampled_dist
        output['javelin']['frac_rejected'] = 1 - len(downsampled_dist) / len(lag_dist)
    
    else:
        #Add to summary dict
        summary_dict['n0_javelin'] = np.nan
        summary_dict['peak_bounds_javelin'] = [np.nan, np.nan]
        summary_dict['peak_javelin'] = np.nan
        summary_dict['lag_javelin'] = np.nan
        summary_dict['lag_err_javelin'] = [np.nan, np.nan]
        summary_dict['frac_rejected_javelin'] = np.nan


        #Add to weighting results
        output['javelin']['tophat_lag'] =  [np.nan, np.nan, np.nan]
        output['javelin']['bounds'] =  [np.nan, np.nan, np.nan]
        output['javelin']['acf'] = np.nan
        output['javelin']['lags'] = np.nan
        output['javelin']['weight_dist'] = np.nan
        output['javelin']['smoothed_dist'] = np.nan
        output['javelin']['ntau'] = np.nan
        output['javelin']['lag_dist'] = np.nan
        output['javelin']['downsampled_lag_dist'] = np.nan
        output['javelin']['frac_rejected'] = np.nan
        

    #---------------------------
    #PyROA        
        
    if run_pyroa:
        _, nburn, _, _, _, add_var, delay_dist, _, _, _ = defaults.set_pyroa(pyroa_params, 2)
        if isinstance(add_var, list):
            add_var = add_var[pyroa_obj_ind-1]
        if isinstance(delay_dist, list):
            delay_dist = delay_dist[pyroa_obj_ind-1]
        
        pyroa_samples = pickle.load(open(pyroa_sample_file, 'rb'))
        samples_chunks = get_samples_chunks( pyroa_samples, nburn, add_var, delay_dist )
        
        lag_dist = samples_chunks[pyroa_obj_ind][2]
        prob_dist, lags, ntau, acf, n0 = get_weights(x_cont, y_cont, x_line, y_line,
                                            interp=interp, lag_bounds=lag_bounds,
                                            sigmode=sigmode, thres=thres,
                                            gap_size=gap_size, k=k)

        min_bound, peak, max_bound, smooth_dist, smooth_weight_dist = get_bounds(lag_dist, prob_dist, lags, width=width)
        downsampled_dist = lag_dist[(lag_dist > min_bound) & (lag_dist < max_bound)]

        med_lag = np.median(downsampled_dist)
        lag_err_lo = med_lag - np.percentile( downsampled_dist, 16 )
        lag_err_hi = np.percentile( downsampled_dist, 84 ) - med_lag


        #Write diagnostic info
        write_data( [ lags, ntau, prob_dist, acf, smooth_dist, smooth_weight_dist ],
                    output_dir + 'pyroa_weights.dat',
                    '#lags,ntau,weight_dist,acf,smooth_dist,smooth_weight_dist')

        write_data( downsampled_dist,
                output_dir + 'pyroa_weighted_lag_dist.dat')



        #Add to summary dict
        summary_dict['n0_pyroa'] = n0
        summary_dict['peak_bounds_pyroa'] = [min_bound, max_bound]
        summary_dict['peak_pyroa'] = peak
        summary_dict['lag_pyroa'] = med_lag
        summary_dict['lag_err_pyroa'] = [lag_err_lo, lag_err_hi]
        summary_dict['frac_rejected_pyroa'] = 1 - len(downsampled_dist) / len(lag_dist)


        #Add to weighting results
        output['pyroa']['time_delay'] =  [lag_err_lo, med_lag, lag_err_hi]
        output['pyroa']['bounds'] =  [min_bound, peak, max_bound]
        output['pyroa']['acf'] = acf
        output['pyroa']['lags'] = lags
        output['pyroa']['weight_dist'] = prob_dist
        output['pyroa']['smoothed_dist'] = smooth_weight_dist
        output['pyroa']['ntau'] = ntau
        output['pyroa']['lag_dist'] = lag_dist
        output['pyroa']['downsampled_lag_dist'] = downsampled_dist
        output['pyroa']['frac_rejected'] = 1 - len(downsampled_dist) / len(lag_dist)
    
    else:
        #Add to summary dict
        summary_dict['n0_pyroa'] = np.nan
        summary_dict['peak_bounds_pyroa'] = [np.nan, np.nan]
        summary_dict['peak_pyroa'] = np.nan
        summary_dict['lag_pyroa'] = np.nan
        summary_dict['lag_err_pyroa'] = [np.nan, np.nan]
        summary_dict['frac_rejected_pyroa'] = np.nan


        #Add to weighting results
        output['pyroa']['time_delay'] =  [np.nan, np.nan, np.nan]
        output['pyroa']['bounds'] =  [np.nan, np.nan, np.nan]
        output['pyroa']['acf'] = np.nan
        output['pyroa']['lags'] = np.nan
        output['pyroa']['weight_dist'] = np.nan
        output['pyroa']['smoothed_dist'] = np.nan
        output['pyroa']['ntau'] = np.nan
        output['pyroa']['lag_dist'] = np.nan
        output['pyroa']['downsampled_lag_dist'] = np.nan
        output['pyroa']['frac_rejected'] = np.nan
        



    #---------------------------
    #Get rmax

    if run_pyccf:
        ccf_lags, ccf = np.loadtxt(pyccf_iccf_file, unpack=True, delimiter=',')


        #PyCCF lag
        lag = output['pyccf']['centroid'][1]
        lag_err_hi = output['pyccf']['centroid'][2]
        lag_err_lo = output['pyccf']['centroid'][0]
        good_ind = np.argwhere( ( ccf_lags >= lag-lag_err_lo ) & ( ccf_lags <= lag+lag_err_hi ) ).T[0]
        
        if len(good_ind) > 0:
            rmax_pyccf = np.max(ccf[good_ind])
        else:
            rmax_pyccf = ccf[ np.argmin( np.abs(lag - ccf_lags) ) ]

        summary_dict['rmax_pyccf'] = rmax_pyccf
        output['rmax_pyccf'] = rmax_pyccf



        if run_javelin:
            #JAVELIN lag
            lag = output['javelin']['tophat_lag'][1]
            lag_err_hi = output['javelin']['tophat_lag'][2]
            lag_err_lo = output['javelin']['tophat_lag'][0]
            good_ind = np.argwhere( ( ccf_lags >= lag-lag_err_lo ) & ( ccf_lags <= lag+lag_err_hi ) ).T[0]
            
            if len(good_ind) > 0:
                rmax_jav = np.max(ccf[good_ind])
            else:
                rmax_jav = ccf[ np.argmin( np.abs(lag - ccf_lags) ) ]

            summary_dict['rmax_javelin'] = rmax_jav
            output['rmax_javelin'] = rmax_jav
        
        
        if run_pyroa:
            #PyROA lag
            lag = output['pyroa']['time_delay'][1]
            lag_err_hi = output['pyroa']['time_delay'][2]
            lag_err_lo = output['pyroa']['time_delay'][0]            
            good_ind = np.argwhere( ( ccf_lags >= lag-lag_err_lo ) & ( ccf_lags <= lag+lag_err_hi ) ).T[0]

            if len(good_ind) > 0:
                rmax_pyroa = np.max(ccf[good_ind])
            else:
                rmax_pyroa = ccf[ np.argmin( np.abs(lag - ccf_lags) ) ]

            summary_dict['rmax_pyroa'] = rmax_pyroa
            output['rmax_pyroa'] = rmax_pyroa


    #---------------------------
    #Write summary file
    write_weighting_summary(output_dir + 'weight_summary.fits', summary_dict, run_pyccf, run_javelin, run_pyroa)

    return output, summary_dict





def run_weighting_tot(output_dir,
                      jav_chain_fnames=None, pyccf_iccf_fnames=None, pyccf_dist_fnames=None,
                      pyroa_sample_fnames=None,
                      line_names=None, interp=2, together_jav=False,
                      pyroa_obj_inds=None, pyroa_params={},
                      general_kwargs={}, weighting_params={}, share_lag_bounds=True):


    output_dir = os.path.abspath(output_dir) + r'/'
    pyccf_params = {'interp':interp}

    if line_names is None:
        warnings.warn('Assuming that the filenames are in the same order as the line names. Line names will be acquired in chronological order from the given directory, except the first will be the continuum', RuntimeWarning)
        line_names = get_ordered_line_names(output_dir)

    _, _, _, _, _, _, _, _, together_pyroa, _ = defaults.set_pyroa( pyroa_params, len(line_names) )

    #---------------------------
    #Get data fnames

    if output_dir + 'processed_lcs/' in glob.glob( output_dir + '*' ):
        line_fnames = np.array([ output_dir + 'processed_lcs/' + x + '_data.dat' for x in line_names ])
    else:
        line_fnames = np.array([ output_dir + 'light_curves/' + x + '.dat' for x in line_names ])


    general_kwargs = defaults.set_general(general_kwargs, line_fnames)
    _, _, _, zoom = defaults.set_weighting(weighting_params)

    #---------------------------
    #Share lag bounds?
    if share_lag_bounds:

        baselines = []

        x_cont, _, _ = np.loadtxt(line_fnames[0], unpack=True, delimiter=',', usecols=[0,1,2])
        for i in range(len(line_fnames)):
            x_line, _, _ = np.loadtxt(line_fnames[i], unpack=True, delimiter=',', usecols=[0,1,2])

            bl = np.max([ x_cont.max(), x_line.max() ]) - np.min([ x_cont.min(), x_line.min() ])
            baselines.append(bl)



        lag_bounds = []
        lag_bounds_i = [-np.max(baselines), np.max(baselines)]
        for i in range(len(line_fnames)-1):
            lag_bounds.append(lag_bounds_i)


    else:
        lag_bounds = general_kwargs['lag_bounds']


    #---------------------------
    #Account for None inputs
    run_javelin = True
    run_pyccf = True
    run_pyroa = True

    if jav_chain_fnames is None:
        warnings.warn('Assuming JAVELIN was not run.', RuntimeWarning)

    if (pyccf_iccf_fnames is None) or (pyccf_iccf_fnames is None):
        warnings.warn('Assuming PyCCF was not run.', RuntimeWarning)


    if jav_chain_fnames is None:
        jav_chain_fnames = [None] * ( len(line_names)-1)
        run_javelin = False
    elif together_jav:
        jav_chain_fnames = [jav_chain_fnames] * ( len(line_names)-1)
        run_javelin = True
        

    if pyccf_iccf_fnames is None:
        pyccf_iccf_fnames = [None] * ( len(line_names)-1)
        run_pyccf = False

    if pyccf_dist_fnames is None:
        pyccf_dist_fnames = [None] * ( len(line_names)-1)
        run_pyccf = False
        
        
        
    #Look for PyROA files if none input
    if pyroa_sample_fnames is None:
        warnings.warn('Looking for PyROA files, no files input', RuntimeWarning)


        if together_pyroa:
            fnames_i = glob.glob(output_dir + 'pyroa/*.obj')
            if len(fnames_i) == 0:
                warnings.warn('No PyROA files found. Assuming PyROA was not run.', RuntimeWarning)
                run_pyroa = False
            else:
                pyroa_sample_fnames = output_dir + 'pyroa/samples.obj'    
            

        else:

            pyroa_sample_fnames = []        
            for i in range(len(line_fnames[1:])):
                fnames_i = glob.glob(output_dir + r'/' + line_names[i+1] + 'pyroa/*.obj')

                if len(fnames_i) == 0:
                    warnings.warn('No PyROA files found for ' + line_names[i+1], RuntimeWarning)
                    pyroa_sample_fnames.append( None )
                
                else: 
                   pyroa_sample_fnames.append( output_dir + r'/' + line_names[i+1] + 'pyroa/samples.obj' )
    
            if np.all([ x is None for x in pyroa_sample_fnames ]):
                warnings.warn('No PyROA files found. Assuming PyROA was not run.', RuntimeWarning)
                run_pyroa = False
                pyroa_sample_fnames = None
            

    #Check if PyROA files exist
    if isinstance(pyroa_sample_fnames, list):    
        for i in range(len(pyroa_sample_fnames)):
            
            if os.path.isfile(pyroa_sample_fnames[i]):
                continue
            else:
               warnings.warn( pyroa_sample_fnames[i] + ' not found.', RuntimeWarning ) 
               pyroa_sample_fnames[i] = None

        if np.all([ x is None for x in pyroa_sample_fnames ]):
            warnings.warn('No PyROA files found. Assuming PyROA was not run.', RuntimeWarning)
            run_pyroa = False
            pyroa_sample_fnames = None
                     

    elif isinstance(pyroa_sample_fnames, str):
        if not os.path.isfile(pyroa_sample_fnames):
            warnings.warn( pyroa_sample_fnames + ' not found. Assuming PyROA was not run.', RuntimeWarning ) 
            run_pyroa = False
            pyroa_sample_fnames = None

            
    if pyroa_sample_fnames is None:
        pyroa_sample_fnames = [None] * ( len(line_names)-1)
        run_pyroa=False
    elif together_pyroa:
        pyroa_sample_fnames = [pyroa_sample_fnames] * ( len(line_names)-1)
        run_pyroa = True

        
    if run_pyroa:
        if pyroa_obj_inds is None:
            if together_pyroa:
                pyroa_obj_inds = range(1, len(line_names))
            else:
                pyroa_obj_inds = []
                
                for i in range(len(line_names)-1):
                    if pyroa_sample_fnames[i] is None:
                        pyroa_obj_inds.append(None)
                    else:
                        pyroa_obj_inds.append(1)
    else:
        pyroa_obj_inds = [None] * ( len(line_names)-1)

    #---------------------------
    #Make weights directories
    for name in line_names[1:]:
        os.makedirs(output_dir + name + r'/weights/', exist_ok=True)

    #---------------------------
    #Run weighting

    summary_dicts = []
    outputs = []

    summary_fnames = []


    for i in range(len(line_fnames)-1):

        if together_jav:
            javelin_lag_col = 2 + 3*i
        else:
            javelin_lag_col = 2 

        res, summary_dict = run_weighting_single(output_dir + line_names[i+1] + r'/weights/', 
                                                 line_fnames[0], line_fnames[i+1],
                                                 weighting_params, lag_bounds[i],
                                                 jav_chain_fnames[i], pyccf_iccf_fnames[i], pyccf_dist_fnames[i],
                                                 pyroa_sample_fnames[i],
                                                 javelin_lag_col=javelin_lag_col,
                                                 pyroa_obj_ind=pyroa_obj_inds[i],
                                                 pyccf_params=pyccf_params, pyroa_params=pyroa_params)

        summary_dicts.append(summary_dict)
        outputs.append(res)
        summary_fnames.append(output_dir + name + r'/weights/weight_summary.fits')



    #---------------------------
    #Get total results

    res_tot = combine_weight_outputs(outputs, run_pyccf, run_javelin, run_pyroa)

    for i in range(len(line_fnames)-1):
        plot_weights(output_dir, line_names[i+1], res_tot['pyccf'], 
                        summary_dicts[i]['n0_pyccf'], summary_dicts[i]['k'],
                        general_kwargs['time_unit'], general_kwargs['plot'])

    if run_pyccf:
        plot_weight_output(output_dir, line_fnames[0], line_fnames, line_names,
                            res_tot['pyccf'], general_kwargs, 'pyccf', zoom,
                            general_kwargs['plot'])

    if run_javelin:
        plot_weight_output(output_dir, line_fnames[0], line_fnames, line_names,
                            res_tot['javelin'], general_kwargs, 'javelin', zoom,
                            general_kwargs['plot'])
        
    if run_pyroa:
        plot_weight_output(output_dir, line_fnames[0], line_fnames, line_names,
                            res_tot['pyroa'], general_kwargs, 'pyroa', zoom,
                            general_kwargs['plot'])


    return res_tot, summary_dicts



#############################################
############## PLOTTING #####################
#############################################


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
