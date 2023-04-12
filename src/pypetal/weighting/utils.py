import pickle

import numpy as np
from scipy.signal import peak_widths

from pypetal.pyccf.utils import peakcent
from pypetal.pyroa.utils import get_samples_chunks
from pypetal.utils import defaults
from pypetal.utils.petalio import write_data, write_weighting_summary

###########################################################################
########################## ASSIST FUNCTIONS ###############################
###########################################################################

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
        _, _, _ = peakcent(x, y, x, y, lag_bounds[0], lag_bounds[1], interp,
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

def get_bounds(dist, weights, lags, width=15, rel_height=.99):


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
    res = peak_widths( smooth_weight_dist, [peak_ind], rel_height=rel_height )

    peak = lags[peak_ind]
    bound_left = lags[ np.floor(res[2]).astype(int) ] + dbin*( res[2]%1 )
    bound_right = lags[ np.floor(res[3]).astype(int) ] + dbin*( res[3]%1 )

    return bound_left[0], peak, bound_right[0], smooth_dist, smooth_weight_dist

###########################################################################
###########################################################################
###########################################################################

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

    cols_pyccf.append('downsampled_CCCD')
    cols_javelin.append('downsampled_lag_dist')
    cols_pyroa.append('downsampled_lag_dist')


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
    gap_size, k, width, rel_height, _ = defaults.set_weighting(weighting_kwargs)

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

        min_bound, peak, max_bound, smooth_dist, smooth_weight_dist = get_bounds(cccd_lags, prob_dist, lags, width=width, rel_height=rel_height)
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

        else:
            summary_dict['rmax_javelin'] = np.nan
            output['rmax_javelin'] = np.nan


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
        else:
            summary_dict['rmax_pyroa'] = np.nan
            output['rmax_pyroa'] = np.nan


    else:
        summary_dict['rmax_pyccf'] = np.nan
        output['rmax_pyccf'] = np.nan

        summary_dict['rmax_javelin'] = np.nan
        output['rmax_javelin'] = np.nan

        summary_dict['rmax_pyroa'] = np.nan
        output['rmax_pyroa'] = np.nan


    #---------------------------
    #Write summary file
    write_weighting_summary(output_dir + 'weight_summary.fits', summary_dict, run_pyccf, run_javelin, run_pyroa)

    return output, summary_dict
