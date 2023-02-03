import os
import re

import astropy.units as u
import numpy as np
from astropy.io import fits

from pypetal.load import get_modules, get_ordered_line_names


def make_directories(output_dir, fnames, line_names,
                     run_drw_rej, run_pyccf, run_pyzdcf,
                     run_javelin, run_weighting,
                     reject_data, together):

    #Create subdirectories for each line and javelin
    for i in range(len(fnames)):
        os.makedirs( output_dir + line_names[i], exist_ok=True )

        if (run_drw_rej) & (reject_data[i]):
            os.makedirs( output_dir + line_names[i] + '/drw_rej', exist_ok=True )

    for i in range(len(fnames)-1):
        if run_pyccf:
            os.makedirs( output_dir + line_names[i+1] + '/pyccf', exist_ok=True )
        if run_pyzdcf:
            os.makedirs( output_dir + line_names[i+1] + '/pyzdcf', exist_ok=True )

        if run_javelin:
            if together:
                os.makedirs( output_dir + 'javelin/', exist_ok=True )
            else:
                os.makedirs( output_dir + line_names[i+1] + '/javelin', exist_ok=True )

        if ( (run_javelin) & (run_weighting) & (not together) ) | (run_pyccf & run_weighting):
            os.makedirs( output_dir + line_names[i+1] + '/weights', exist_ok=True )


    if run_javelin & together & run_weighting:
        os.makedirs( output_dir + 'javelin/weights', exist_ok=True )


    #Make subdirectories for light curves
    os.makedirs( output_dir + 'light_curves/', exist_ok=True )

    return






def str2unit(x):

    if (x == 'Arbitrary Units'):
        unit = u.dimensionless_unscaled
    else:
        unit = u.Unit(x)

    return unit


#From https://stackoverflow.com/questions/5807952/removing-trailing-zeros-in-python/5808661#5808661
def number_shaver(ch,
                  regx = re.compile('(?<![\d.])0*(?:'
                                    '(\d+)\.?|\.(0)'
                                    '|(\.\d+?)|(\d+\.\d+?)'
                                    ')0*(?![\d.])')  ,
                  repl = lambda mat: mat.group(mat.lastindex)
                                     if mat.lastindex!=3
                                     else '0' + mat.group(3) ):
    return regx.sub(repl,ch)

def err2str(val, up_err, lo_err, dec=2):
    val = number_shaver( str(round(val, dec)) )
    up_err = number_shaver( str(round(up_err, dec)) )
    lo_err = number_shaver( str(round(lo_err, dec)) )

    return val + '^{+' + up_err + '}_{-' + lo_err + '}'


def write_data(arr, fname, header=None):

    arr = np.array(arr, dtype=object)

    ndim = len( arr.shape )
    assert ndim <= 2


    if ndim == 2:
        cols = len(arr)
        rows = len(arr[0])

        with open(fname, "w") as file:

            if header is not None:
                file.write(header + "\n")

            for i in range(rows):

                string = "{},".format(arr[0][i])
                for j in range(1, cols-1):
                    string += "{},".format(arr[j][i])

                string += "{}\n".format(arr[-1][i])

                file.write(string)

    elif ndim == 1:
        rows = len(arr)

        with open(fname, "w") as file:

            if header is not None:
                file.write(header + "\n")

            for i in range(rows):
                file.write( "{}\n".format(arr[i]) )

    return



def write_weight_summary(fname, res):

    k = res['k']

    #pyCCF
    n0_pyccf = res['n0_pyccf']
    peak_bounds_pyccf = res['peak_bounds_pyccf']
    peak_pyccf = res['peak_pyccf']
    lag_pyccf = res['lag_pyccf']
    lag_err_pyccf = res['lag_err_pyccf']
    frac_rejected_pyccf = res['frac_rejected_pyccf']

    #JAVELIN
    n0_javelin = res['n0_javelin']
    peak_bounds_javelin = res['peak_bounds_javelin']
    peak_javelin = res['peak_javelin']
    lag_javelin = res['lag_javelin']
    lag_err_javelin = res['lag_err_javelin']
    frac_rejected_javelin = res['frac_rejected_javelin']

    #Total
    rmax = res['rmax']


    with open(fname, 'w') as file:
        file.write("k = {}\n".format(k))
        file.write("\n")
        file.write("pyCCF\n")
        file.write("----------------\n")
        file.write("n0 = {}\n".format(n0_pyccf))
        file.write("peak bounds = {}\n".format(peak_bounds_pyccf))
        file.write("peak = {}\n".format(peak_pyccf))
        file.write("lag value = {}\n".format(lag_pyccf))
        file.write("lag uncertainty = {}\n".format(lag_err_pyccf))
        file.write("fraction rejected = {}\n".format(frac_rejected_pyccf))
        file.write("\n")
        file.write("JAVELIN\n")
        file.write("----------------\n")
        file.write("n0 = {}\n".format(n0_javelin))
        file.write("peak bounds = {}\n".format(peak_bounds_javelin))
        file.write("peak = {}\n".format(peak_javelin))
        file.write("lag value = {}\n".format(lag_javelin))
        file.write("lag uncertainty = {}\n".format(lag_err_javelin))
        file.write("fraction rejected = {}\n".format(frac_rejected_javelin))
        file.write("\n")
        file.write("Total\n")
        file.write("----------------\n")
        file.write("rmax = {}\n".format(rmax))

    return


def write_to_fits(res, output_dir, line_names=None):

    """Write all output information from a given module into a single FITS file.
    """

    if line_names is None:
        line_names = get_ordered_line_names(line_names)

    run_drw_rej, run_pyccf, run_pyzdcf, run_javelin, run_weighting = get_modules(output_dir)
    nlc = len(line_names)

    if run_drw_rej:

        for i in range(nlc):
            #Create table for MC sample from DRW rejection
            c1 = fits.Column( name='tau_sample', array=res['drw_rej_res']['taus'][i], format='F' )
            c2 = fits.Column( name='sigma_sample', array=res['drw_rej_res']['sigmas'][i], format='F' )

            if res['drw_rej_res']['jitters'][i] is not None:
                c3 = fits.Column( name='jitter_sample', array=res['drw_rej_res']['jitters'][i], format='F' )
            else:
                c3 = fits.Column( name='jitter_sample',
                                 array=np.full( len(res['drw_rej_res']['taus'][i]), np.nan ),
                                 format='F' )


            hdr = fits.Header()
            hdr['TAU_SAMPLE'] = 'MC chains for tau'
            hdr['SIGMA_SAMPLE'] = 'MC chains for sigma'
            hdr['JITTER_SAMPLE'] = 'MC chains for jitter'
            drw_rej_table = fits.BinTableHDU.from_columns( [c1, c2, c3], header=hdr )


            #Create table for masks from DRW rejection
            c1 = fits.Column( name='mask', array=np.array( res['drw_rej_res']['masks'][i], dtype=str), format='5A' )

            hdr = fits.Header()
            hdr['MASK'] = 'Rejection mask where True = rejected'
            drw_mask_table = fits.BinTableHDU.from_columns( [c1], header=hdr )


            #Create table for DRW rejection results
            tau_med = np.median( res['drw_rej_res']['taus'][i] )
            tau_err_hi = np.percentile( res['drw_rej_res']['taus'][i], 84 ) - tau_med
            tau_err_lo = tau_med - np.percentile( res['drw_rej_res']['taus'][i], 16 )
            c1 = fits.Column( name='tau', array=[tau_err_lo, tau_med, tau_err_hi], format='F' )

            sigma_med = np.median( res['drw_rej_res']['sigmas'][i] )
            sigma_err_hi = np.percentile( res['drw_rej_res']['sigmas'][i], 84 ) - sigma_med
            sigma_err_lo = sigma_med - np.percentile( res['drw_rej_res']['sigmas'][i], 16 )
            c2 = fits.Column( name='sigma', array=[sigma_err_lo, sigma_med, sigma_err_hi], format='F' )

            if res['drw_rej_res']['jitters'][i] is not None:
                jitter_med = np.median( res['drw_rej_res']['jitters'][i] )
                jitter_err_hi = np.percentile( res['drw_rej_res']['jitters'][i], 84 ) - jitter_med
                jitter_err_lo = jitter_med - np.percentile( res['drw_rej_res']['jitters'][i], 16 )
                c3 = fits.Column( name='jitter', array=[jitter_err_lo, jitter_med, jitter_err_hi], format='F' )
            else:
                c3 = fits.Column( name='jitter', array=[0, 0, 0], format='F' )


            hdr = fits.Header()
            hdr['TAU'] = 'Median and 1-sigma errors for tau'
            hdr['SIGMA'] = 'Median and 1-sigma errors for sigma'
            hdr['JITTER'] = 'Median and 1-sigma errors for jitter'
            drw_fits_table = fits.BinTableHDU.from_columns( [c1, c2, c3], header=hdr )

            #Create primary HDU
            hdr = fits.Header()
            hdr['NAME'] = line_names[i]
            drw_primary_hdu = fits.PrimaryHDU(header=hdr)

            #Make HDU List
            hdul = fits.HDUList([drw_primary_hdu, drw_rej_table, drw_mask_table, drw_fits_table])

            #Write to file
            hdul.writeto( output_dir + line_names[i] + '/drw_rej.fits', overwrite=True )


    if run_pyccf:

        for i in range(nlc-1):

            #Create table for ICCF
            c1 = fits.Column( name='lags', array=res['pyccf_res'][i]['CCF_lags'], format='F' )
            c2 = fits.Column( name='CCF', array=res['pyccf_res'][i]['CCF'], format='F' )

            hdr = fits.Header()
            hdr['LAGS'] = 'Lags for CCF'
            hdr['CCF'] = 'Cross-correlation function from pyCCF'
            pyccf_ccf_table = fits.BinTableHDU.from_columns( [c1, c2], header=hdr )



            #Create table for peak and centroid
            cent_med = res['pyccf_res'][i]['centroid']
            cent_err_lo = res['pyccf_res'][i]['centroid_err_lo']
            cent_err_hi = res['pyccf_res'][i]['centroid_err_hi']

            peak_med = res['pyccf_res'][i]['peak']
            peak_err_lo = res['pyccf_res'][i]['peak_err_lo']
            peak_err_hi = res['pyccf_res'][i]['peak_err_hi']

            c1 = fits.Column( name='centroid', array=[cent_med, cent_err_lo, cent_err_hi], format='F' )
            c2 = fits.Column( name='peak', array=[peak_med, peak_err_lo, peak_err_hi], format='F' )

            hdr = fits.Header()
            hdr['CENTROID'] = 'Median and 1-sigma errors for centroid'
            hdr['PEAK'] = 'Median and 1-sigma errors for peak'

            pyccf_lag_table = fits.BinTableHDU.from_columns( [c1, c2], header=hdr )


            #Create CCCD/CCPD table
            c1 = fits.Column( name='CCCD', array=res['pyccf_res'][i]['CCCD_lags'], format='F' )
            c2 = fits.Column( name='CCPD', array=res['pyccf_res'][i]['CCPD_lags'], format='F' )

            hdr = fits.Header()
            hdr['CCCD'] = 'The CCCD distribution'
            hdr['CCPD'] = 'The CCPD distribution'

            pyccf_dist_table = fits.BinTableHDU.from_columns( [c1, c2], header=hdr )


            #Create primary HDU
            hdr = fits.Header()
            hdr['NAME'] = line_names[i+1]
            pyccf_primary_hdu = fits.PrimaryHDU(header=hdr)


            #Make HDU List
            hdul = fits.HDUList([pyccf_primary_hdu, pyccf_ccf_table, pyccf_lag_table, pyccf_dist_table])

            #Write to file
            hdul.writeto( output_dir + line_names[i+1] + '/pyccf.fits', overwrite=True )



    if run_pyzdcf:

        for i in range(nlc-1):

            #Create table for ZDCF
            c1 = fits.Column(name='tau', array=res['pyzdcf_res'][i]['tau'], format='F')
            c2 = fits.Column(name='tau_err_lo', array=res['pyzdcf_res'][i]['-sig(tau)'], format='F')
            c3 = fits.Column(name='tau_err_hi', array=res['pyzdcf_res'][i]['+sig(tau)'], format='F')
            c4 = fits.Column(name='ZDCF', array=res['pyzdcf_res'][i]['dcf'], format='F')
            c5 = fits.Column(name='ZDCF_err_lo', array=res['pyzdcf_res'][i]['-err(dcf)'], format='F')
            c6 = fits.Column(name='ZDCF_err_hi', array=res['pyzdcf_res'][i]['+err(dcf)'], format='F')

            hdr = fits.Header()
            hdr['TAU'] = 'Lags for the ZDCF'
            hdr['TAU_ERR_LO'] = 'Lower error on tau'
            hdr['TAU_ERR_HI'] = 'Upper error on tau'
            hdr['ZDCF'] = 'The ZDCF'
            hdr['ZDCF_ERR_LO'] = 'Lower error on ZDCF'
            hdr['ZDCF_ERR_HI'] = 'Upper error on ZDCF'

            pyzdcf_zdcf_table = fits.BinTableHDU.from_columns( [c1, c2, c3, c4, c5, c6], header=hdr )


            #Create primary HDU
            hdr = fits.Header()
            hdr['NAME'] = line_names[i+1]
            pyzdcf_primary_hdu = fits.PrimaryHDU(header=hdr)

            #Create HDU List
            hdul = fits.HDUList([pyzdcf_primary_hdu, pyzdcf_zdcf_table])



            #Create tables for PLIKE
            if run_plike:
                hi
                #NEED TO ADD THIS


            #Save to file
            hdul.writeto( output_dir + line_names[i+1] + '/pyzdcf.fits', overwrite=True )



    if run_javelin:

        if together:
            #NEED TO ADD THIS
            hi

        else:
            for i in range(nlc-1):

                #Create cont HPD table
                if res['javelin_res'][i]['cont_hpd'] is not None:
                    c1 = fits.column(name='sigma_hpd', array=res['javelin_res'][i]['cont_hpd'][0], format='F')
                    c2 = fits.column(name='tau_hpd', array=res['javelin_res'][i]['cont_hpd'][1], format='F')

                    hdr = fits.Header()
                    hdr['SIGMA_HPD'] = 'The median and 1-sigma errors for sigma'
                    hdr['TAU_HPD'] = 'The median and 1-sigma errors for tau'

                    javelin_cont_hpd_table = fits.BinTableHDU.from_columns( [c1, c2], header=hdr )


                #Create MCMC samples table
                c1 = fits.column(name='tau', array=res['javelin_res'][i]['tau'], format='F')
                c2 = fits.column(name='sigma', array=res['javelin_res'][i]['sigma'], format='F')
                c3 = fits.column(name='lag', array=res['javelin_res'][i]['tophat_params'][0], format='F')
                c4 = fits.column(name='width', array=res['javelin_res'][i]['tophat_params'][1], format='F')
                c5 = fits.column(name='scale', array=res['javelin_res'][i]['tophat_params'][2], format='F')

                hdr = fits.Header()
                hdr['TAU'] = 'The MCMC samples for tau'
                hdr['SIGMA'] = 'The MCMC samples for sigma'
                hdr['LAG'] = 'The MCMC samples for the tophat lag'
                hdr['WIDTH'] = 'The MCMC samples for the tophat width'
                hdr['SCALE'] = 'The MCMC samples for the tophat scale'

                javelin_mcmc_table = fits.BinTableHDU.from_columns( [c1, c2, c3, c4, c5], header=hdr )




                #Create rmap HPD table
                c1 = fits.column(name='sigma_hpd', array=res['javelin_res'][i]['rmap_hpd'][0], format='F')
                c2 = fits.column(name='tau_hpd', array=res['javelin_res'][i]['rmap_hpd'][1], format='F')
                c3 = fits.column(name='lag_hpd', array=res['javelin_res'][i]['rmap_hpd'][2], format='F')
                c4 = fits.column(name='width_hpd', array=res['javelin_res'][i]['rmap_hpd'][3], format='F')
                c5 = fits.column(name='scale_hpd', array=res['javelin_res'][i]['rmap_hpd'][4], format='F')

                hdr = fits.Header()
                hdr['SIGMA_HPD'] = 'The median and 1-sigma errors for sigma'
                hdr['TAU_HPD'] = 'The median and 1-sigma errors for tau'
                hdr['LAG_HPD'] = 'The median and 1-sigma errors for the tophat lag'
                hdr['WIDTH_HPD'] = 'The median and 1-sigma errors for the tophat width'
                hdr['SCALE_HPD'] = 'The median and 1-sigma errors for the tophat scale'

                javelin_rmap_hpd_table = fits.BinTableHDU.from_columns( [c1, c2, c3, c4, c5], header=hdr )



                #Create bestfit model table

                #Continuum
                xvals = res['javelin_res'][i]['bestfit_model'].jlist[0]
                yvals = res['javelin_res'][i]['bestfit_model'].mlist[0] + res['javelin_res'][i]['bestfit_model'].blist[0]
                yerr_vals = res['javelin_res'][i]['bestfit_model'].elist[0]

                c1 = fits.column(name='cont_times', array=xvals, format='F')
                c2 = fits.column(name='cont_mean', array=yvals, format='F')
                c3 = fits.column(name='cont_error', array=yerr_vals, format='F')

                hdr = fits.Header()
                hdr['CONT_TIMES'] = 'The times for the bestfit continuum model'
                hdr['CONT_MEAN'] = 'The mean fit for the bestfit continuum model'
                hdr['CONT_ERROR'] = 'The error for the bestfit continuum model'


                #Line
                xvals = res['javelin_res'][i]['bestfit_model'].jlist[1]
                yvals = res['javelin_res'][i]['bestfit_model'].mlist[1] + res['javelin_res'][i]['bestfit_model'].blist[1]
                yerr_vals = res['javelin_res'][i]['bestfit_model'].elist[1]

                c4 = fits.column(name='times', array=xvals, format='F')
                c5 = fits.column(name='mean', array=yvals, format='F')
                c6 = fits.column(name='error', array=yerr_vals, format='F')

                hdr = fits.Header()
                hdr['TIMES'] = 'The times for the bestfit line model'
                hdr['MEAN'] = 'The mean fit for the bestfit line model'
                hdr['ERROR'] = 'The error for the bestfit line model'



                javelin_bestfit_table = fits.BinTableHDU.from_columns( [c1, c2, c3, c4, c5, c6],
                                                                       header=hdr )


                #Create primary HDU
                hdr = fits.Header()
                hdr['NAME'] = line_names[i+1]
                javelin_primary_hdu = fits.PrimaryHDU(header=hdr)

                #Create HDU List
                hdul = fits.HDUList([javelin_primary_hdu, javelin_cont_hpd_table, javelin_mcmc_table, javelin_rmap_hpd_table, javelin_bestfit_table])

                #Write to file
                hdul.writeto(output_dir + line_names[i+1] + '/javelin.fits', overwrite=True)


    if run_weighting:

        for i in range(nlc-1):

            ################
            #    pyCCF
            ################

            #Create weight table
            c1 = fits.Column( name='lags', array=res['weighting_res']['pyccf']['lags'][i], format='F' )
            c2 = fits.Column( name='acf', array=res['weighting_res']['pyccf']['acf'][i], format='F')
            c3 = fits.Column( name='ntau', array=res['weighting_res']['pyccf']['ntau'][i], format='F' )
            c4 = fits.Column( name='weight_dist', array=res['weighting_res']['pyccf']['weight_dist'][i], format='F' )
            c5 = fits.Column( name='smoothed_weight_dist', array=res['weighting_res']['pyccf']['smoothed_dist'][i], format='F' )

            hdr = fits.Header()
            hdr['MODULE'] = 'pyCCF'
            hdr['LAGS'] = 'The lags for the weight distributions'
            hdr['ACF'] = 'The autocorrelation function for the continuum'
            hdr['NTAU'] = 'The non-normalized probability weighting N(tau)'
            hdr['WEIGHT_DIST'] = 'The weight distribution described in Grier et al. (2019)'
            hdr['SMOOTHED_WEIGHT_DIST'] = 'The weight distribution smoothed with a Gaussian kernel'

            pyccf_weight_table = fits.BinTableHDU.from_columns( [c1, c2, c3, c4, c5], header=hdr )


            #Create downsampled dist table
            c1 = fits.Column( name='downsampled_cccd', array=res['weighting_res']['pyccf']['downsampled_CCCD'][i], format='F' )

            hdr = fits.Header()
            hdr['MODULE'] = 'pyCCF'
            hdr['DOWNSAMPLED_CCCD'] = 'The downsampled CCCD after finding the peak'
            hdr['FRAC_REJECTED'] = res['weighting_res']['pyccf']['frac_rejected'][i]

            pyccf_downsample_table = fits.BinTableHDU.from_columns( [c1], header=hdr )


            #Create peak table
            c1 = fits.Column( name='centroid', array=res['weighting_res']['pyccf']['centroid'][i], format='F' )
            c2 = fits.Column( name='primary_peak', array=res['weighting_res']['pyccf']['bounds'][i], format='F' )

            hdr = fits.Header()
            hdr['MODULE'] = 'pyCCF'
            hdr['CENTROID'] = 'The median of the downsampled CCCD, with its uncertainty given as [lower error, value, upper error]'
            hdr['PRIMARY_PEAK'] = 'The primary peak of the weighted and smoothed CCCD, given as [lower bound, peak, upper bound]'

            pyccf_peak_table = fits.BinTableHDU.from_columns( [c1, c2], header=hdr )




            ################
            #    JAVELIN
            ################

            #Create weight table
            c1 = fits.Column( name='lags', array=res['weighting_res']['javelin']['lags'][i], format='F' )
            c2 = fits.Column( name='acf', array=res['weighting_res']['javelin']['acf'][i], format='F')
            c3 = fits.Column( name='ntau', array=res['weighting_res']['javelin']['ntau'][i], format='F' )
            c4 = fits.Column( name='weight_dist', array=res['weighting_res']['javelin']['weight_dist'][i], format='F' )
            c5 = fits.Column( name='smoothed_weight_dist', array=res['weighting_res']['javelin']['smoothed_dist'][i], format='F' )

            hdr = fits.Header()
            hdr['MODULE'] = 'JAVELIN'
            hdr['LAGS'] = 'The lags for the weight distributions'
            hdr['ACF'] = 'The autocorrelation function for the continuum'
            hdr['NTAU'] = 'The non-normalized probability weighting N(tau)'
            hdr['WEIGHT_DIST'] = 'The weight distribution described in Grier et al. (2019)'
            hdr['SMOOTHED_WEIGHT_DIST'] = 'The weight distribution smoothed with a Gaussian kernel'

            javelin_weight_table = fits.BinTableHDU.from_columns( [c1, c2, c3, c4, c5], header=hdr )


            #Create downsampled dist table
            c1 = fits.Column( name='downsampled_lag_dist', array=res['weighting_res']['javelin']['downsampled_lag_dist'][i], format='F' )

            hdr = fits.Header()
            hdr['MODULE'] = 'pyCCF'
            hdr['DOWNSAMPLED_CCCD'] = 'The downsampled CCCD after finding the peak'
            hdr['FRAC_REJECTED'] = res['weighting_res']['javelin']['frac_rejected'][i]

            javelin_downsample_table = fits.BinTableHDU.from_columns( [c1], header=hdr )


            #Create peak table
            c1 = fits.Column( name='tophat_lag', array=res['weighting_res']['javelin']['tophat_lag'][i], format='F' )
            c2 = fits.Column( name='primary_peak', array=res['weighting_res']['javelin']['bounds'][i], format='F' )

            hdr = fits.Header()
            hdr['MODULE'] = 'pyCCF'
            hdr['TOPHAT_LAG'] = 'The median of the downsampled lag distribution, with its uncertainty given as [lower error, value, upper error]'
            hdr['PRIMARY_PEAK'] = 'The primary peak of the weighted and smoothed lag_distribution, given as [lower bound, peak, upper bound]'
            hdr['RMAX'] = res['weighting_res']['rmax'][i]

            javelin_peak_table = fits.BinTableHDU.from_columns( [c1, c2], header=hdr )


            ################
            #    Total
            ################

            #Create primary HDU
            hdr = fits.Header()
            hdr['NAME'] = line_names[i+1]
            primary_hdu = fits.PrimaryHDU(header=hdr)

            hdul = fits.HDUList([primary_hdu,
                                 pyccf_weight_table, pyccf_downsample_table, pyccf_peak_table,
                                 javelin_weight_table, javelin_downsample_table, javelin_peak_table])

    return
