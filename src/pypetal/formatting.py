import numpy as np
import re

import astropy.units as u
from astropy.io import fits





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
            
            
            
def write_to_fits(res, output_dir, line_names, run_drw_rej, run_pyccf, run_pyzdcf, run_javelin):

    """Write all output information from a given module into a single FITS file.
    """

    nlc = len(line_names)

    if run_drw_rej:
        
        for i in range(nlc):
            #Create table for MC sample from DRW rejection
            c1 = fits.Column( name='tau_sample', array=res['drw_rej_res']['taus'][i], format='F' )
            c2 = fits.Column( name='sigma_sample', array=res['drw_rej_res']['sigmas'][i], format='F' ) 
            
            if res['drw_rej_res']['jitters'][i] is not None:      
                c3 = fits.Column( name='jitter_sample', array=res['drw_rej_res']['jitters'][i], format='F' )
            else:
                c3 = fits.Column( name='jitter_sample', array=res['drw_rej_res']['jitters'][i], format='F' )
        
        
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
                xvals = res['javelin_res'][i]['bestfit_model'].jlist[i]
                yvals = res['javelin_res'][i]['bestfit_model'].mlist[i] + res['javelin_res'][i]['bestfit_model'].blist[i]
                yerr_vals = res['javelin_res'][i]['bestfit_model'].elist[i]
                
                c1 = fits.column(name='times', array=xvals, format='F')
                c2 = fits.column(name='mean', array=yvals, format='F')
                c3 = fits.column(name='error', array=yerr_vals, format='F')
                
                hdr = fits.Header()
                hdr['TIMES'] = 'The times for the bestfit model'
                hdr['MEAN'] = 'The mean fit for the bestfit model'
                hdr['ERROR'] = 'The error for the bestfit model'
                
                javelin_bestfit_table = fits.BinTableHDU.from_columns( [c1, c2, c3], header=hdr )
                
                
                #Create primary HDU
                hdr = fits.Header()
                hdr['NAME'] = line_names[i+1]
                javelin_primary_hdu = fits.PrimaryHDU(header=hdr)
                
                #Create HDU List
                hdul = fits.HDUList([javelin_primary_hdu, javelin_cont_hpd_table, javelin_mcmc_table, avelin_rmap_hpd_table, javelin_bestfit_table])
                
                #Write to file
                hdul.writeto(output_dir + line_names[i+1] + '/javelin.fits', overwrite=True)

    return    