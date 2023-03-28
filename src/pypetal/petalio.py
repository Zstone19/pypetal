import os
import re

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.table import Table


def make_directories(output_dir, fnames, line_names,
                     run_drw_rej, run_pyccf, run_pyzdcf, run_pyroa,
                     reject_data, together_pyroa):

    #Create subdirectories for each line
    for i in range(len(fnames)):
        os.makedirs( output_dir + line_names[i], exist_ok=True )

        if (run_drw_rej) & (reject_data[i]):
            os.makedirs( output_dir + line_names[i] + '/drw_rej', exist_ok=True )

    for i in range(len(fnames)-1):
        if run_pyccf:
            os.makedirs( output_dir + line_names[i+1] + '/pyccf', exist_ok=True )
        if run_pyzdcf:
            os.makedirs( output_dir + line_names[i+1] + '/pyzdcf', exist_ok=True )

    if run_pyroa:
        if together_pyroa:
            os.makedirs( output_dir + 'pyroa', exist_ok=True )
        else:
            for i in range(len(fnames)-1):
                os.makedirs( output_dir + line_names[i+1] + '/pyroa', exist_ok=True )


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



def write_weighting_summary(fname, res, run_pyccf, run_javelin, run_pyroa):

    k = res['k']

    #pyCCF
    if run_pyccf:
        n0_pyccf = res['n0_pyccf']
        peak_bounds_pyccf = res['peak_bounds_pyccf']
        peak_pyccf = res['peak_pyccf']
        lag_pyccf = res['lag_pyccf']
        lag_err_pyccf = res['lag_err_pyccf']
        frac_rejected_pyccf = res['frac_rejected_pyccf']
    else:
        n0_pyccf = np.nan
        peak_bounds_pyccf = [np.nan, np.nan]
        peak_pyccf = np.nan
        lag_pyccf = np.nan
        lag_err_pyccf = [np.nan, np.nan]
        frac_rejected_pyccf = np.nan

    #JAVELIN
    if run_javelin:
        n0_javelin = res['n0_javelin']
        peak_bounds_javelin = res['peak_bounds_javelin']
        peak_javelin = res['peak_javelin']
        lag_javelin = res['lag_javelin']
        lag_err_javelin = res['lag_err_javelin']
        frac_rejected_javelin = res['frac_rejected_javelin']
    else:
        n0_javelin = np.nan
        peak_bounds_javelin = [np.nan, np.nan]
        peak_javelin = np.nan
        lag_javelin = np.nan
        lag_err_javelin = [np.nan, np.nan]
        frac_rejected_javelin = np.nan
        
        
    #PyROA
    if run_pyroa:
        n0_pyroa = res['n0_pyroa']
        peak_bounds_pyroa = res['peak_bounds_pyroa']
        peak_pyroa = res['peak_pyroa']
        lag_pyroa = res['lag_pyroa']
        lag_err_pyroa = res['lag_err_pyroa']
        frac_rejected_pyroa = res['frac_rejected_pyroa']
    else:
        n0_pyroa = np.nan
        peak_bounds_pyroa = [np.nan, np.nan]
        peak_pyroa = np.nan
        lag_pyroa = np.nan
        lag_err_pyroa = [np.nan, np.nan]
        frac_rejected_pyroa = np.nan
        


    #Total
    if run_pyccf:
        rmax_pyccf = res['rmax_pyccf']
        
        if run_javelin:
            rmax_jav = res['rmax_javelin']
        if run_pyroa:
            rmax_pyroa = res['rmax_pyroa']
        
    if not run_pyccf:
        rmax_pyccf = np.nan
        rmax_jav = np.nan
        rmax_pyroa = np.nan
    elif not run_javelin:
        rmax_jav = np.nan
    elif not run_pyroa:
        rmax_pyroa = np.nan    
    

    #---------------------------
    #Make table

    k_col = fits.Column(name='k', format='E', array=[k])

    n0_pyccf_col = fits.Column(name='n0_pyccf', format='E', array=[n0_pyccf])
    peak_bounds_pyccf_col = fits.Column(name='peak_bounds_pyccf', format='2E', array=[peak_bounds_pyccf])
    peak_pyccf_col = fits.Column(name='peak_pyccf', format='E', array=[peak_pyccf])
    lag_pyccf_col = fits.Column(name='lag_pyccf', format='E', array=[lag_pyccf])
    lag_err_pyccf_col = fits.Column(name='lag_err_pyccf', format='2E', array=[lag_err_pyccf])
    frac_rejected_pyccf_col = fits.Column(name='frac_rejected_pyccf', format='E', array=[frac_rejected_pyccf])

    n0_javelin_col = fits.Column(name='n0_javelin', format='E', array=[n0_javelin])
    peak_bounds_javelin_col = fits.Column(name='peak_bounds_javelin', format='2E', array=[peak_bounds_javelin])
    peak_javelin_col = fits.Column(name='peak_javelin', format='E', array=[peak_javelin])
    lag_javelin_col = fits.Column(name='lag_javelin', format='E', array=[lag_javelin])
    lag_err_javelin_col = fits.Column(name='lag_err_javelin', format='2E', array=[lag_err_javelin])
    frac_rejected_javelin_col = fits.Column(name='frac_rejected_javelin', format='E', array=[frac_rejected_javelin])
    
    n0_pyroa_col = fits.Column(name='n0_pyroa', format='E', array=[n0_pyroa])
    peak_bounds_pyroa_col = fits.Column(name='peak_bounds_pyroa', format='2E', array=[peak_bounds_pyroa])
    peak_pyroa_col = fits.Column(name='peak_pyroa', format='E', array=[peak_pyroa])
    lag_pyroa_col = fits.Column(name='lag_pyroa', format='E', array=[lag_pyroa])
    lag_err_pyroa_col = fits.Column(name='lag_err_pyroa', format='2E', array=[lag_err_pyroa])
    frac_rejected_pyroa_col = fits.Column(name='frac_rejected_pyroa', format='E', array=[frac_rejected_pyroa])

    rmax_pyccf_col = fits.Column(name='rmax_pyccf', format='E', array=[rmax_pyccf])
    rmax_jav_col = fits.Column(name='rmax_javelin', format='E', array=[rmax_jav])
    rmax_pyroa_col = fits.Column(name='rmax_pyroa', format='E', array=[rmax_pyroa])


    cols = [k_col, 
            n0_pyccf_col, peak_bounds_pyccf_col, peak_pyccf_col, lag_pyccf_col, lag_err_pyccf_col, frac_rejected_pyccf_col, 
            n0_javelin_col, peak_bounds_javelin_col, peak_javelin_col, lag_javelin_col, lag_err_javelin_col, frac_rejected_javelin_col, 
            n0_pyroa_col, peak_bounds_pyroa_col, peak_pyroa_col, lag_pyroa_col, lag_err_pyroa_col, frac_rejected_pyroa_col,
            rmax_pyccf_col, rmax_jav_col, rmax_pyroa_col]
    table = fits.BinTableHDU.from_columns(cols)
    table.writeto(fname, overwrite=True)

    return 


def combine_weight_summary(filenames, output_fname, line_names=None):

    if line_names is None:
        line_names = []

        for i in range(len(filenames)):
            line_names.append('Line {}'.format(i+1))


    #---------------------------
    #Get header
    hdr = fits.open(filenames[0]).header

    #---------------------------
    #Get data

    k = np.zeros(len(filenames))

    n0_pyccf = np.zeros(len(filenames))
    peak_bounds_pyccf = np.zeros( ( len(filenames),2 ) )
    peak_pyccf = np.zeros(len(filenames))
    lag_pyccf = np.zeros(len(filenames))
    lag_err_pyccf = np.zeros( ( len(filenames),2 ) )
    frac_rejected_pyccf = np.zeros(len(filenames))

    n0_javelin = np.zeros(len(filenames))
    peak_bounds_javelin = np.zeros( (len(filenames),2) )
    peak_javelin = np.zeros(len(filenames))
    lag_javelin = np.zeros(len(filenames))
    lag_err_javelin = np.zeros( (len(filenames),2) )
    frac_rejected_javelin = np.zeros(len(filenames))
    
    n0_pyroa = np.zeros(len(filenames))
    peak_bounds_pyroa = np.zeros( (len(filenames),2) )
    peak_pyroa = np.zeros(len(filenames))
    lag_pyroa = np.zeros(len(filenames))
    lag_err_pyroa = np.zeros( (len(filenames),2) )
    frac_rejected_pyroa = np.zeros(len(filenames))

    rmax_pyccf = np.zeros(len(filenames))
    rmax_jav = np.zeros(len(filenames))
    rmax_pyroa = np.zeros(len(filenames))


    for i, fname in enumerate(filenames):
        table = Table.read(fname)

        k[i] = table['k'][0]

        n0_pyccf[i] = table['n0_pyccf'][0]
        peak_bounds_pyccf[i,:] = table['peak_bounds_pyccf'][0]
        peak_pyccf[i] = table['peak_pyccf'][0]
        lag_pyccf[i] = table['lag_pyccf'][0]
        lag_err_pyccf[i,:] = table['lag_err_pyccf'][0]
        frac_rejected_pyccf[i] = table['frac_rejected_pyccf'][0]

        n0_javelin[i] = table['n0_javelin'][0]
        peak_bounds_javelin[i,:] = table['peak_bounds_javelin'][0]
        peak_javelin[i] = table['peak_javelin'][0]
        lag_javelin[i] = table['lag_javelin'][0]
        lag_err_javelin[i,:] = table['lag_err_javelin'][0]
        frac_rejected_javelin[i] = table['frac_rejected_javelin'][0]
        
        n0_pyroa[i] = table['n0_pyroa'][0]
        peak_bounds_pyroa[i,:] = table['peak_bounds_pyroa'][0]
        peak_pyroa[i] = table['peak_pyroa'][0]
        lag_pyroa[i] = table['lag_pyroa'][0]
        lag_err_pyroa[i,:] = table['lag_err_pyroa'][0]
        frac_rejected_pyroa[i] = table['frac_rejected_pyroa'][0]

        rmax_pyccf[i] = table['rmax_pyccf'][0]
        rmax_jav[i] = table['rmax_javelin'][0] 
        rmax_pyroa[i] = table['rmax_pyroa'][0]

    #---------------------------
    #Make columns

    name_col = fits.Column(name='name', format='20A', array=line_names)
    k_col = fits.Column(name='k', format='E', array=k)
    
    n0_pyccf_col = fits.Column(name='n0_pyccf', format='E', array=n0_pyccf)
    peak_bounds_pyccf_col = fits.Column(name='peak_bounds_pyccf', format='2E', array=peak_bounds_pyccf)
    peak_pyccf_col = fits.Column(name='peak_pyccf', format='E', array=peak_pyccf)
    lag_pyccf_col = fits.Column(name='lag_pyccf', format='E', array=lag_pyccf)
    lag_err_pyccf_col = fits.Column(name='lag_err_pyccf', format='2E', array=lag_err_pyccf)
    frac_rejected_pyccf_col = fits.Column(name='frac_rejected_pyccf', format='E', array=frac_rejected_pyccf)

    n0_javelin_col = fits.Column(name='n0_javelin', format='E', array=n0_javelin)
    peak_bounds_javelin_col = fits.Column(name='peak_bounds_javelin', format='2E', array=peak_bounds_javelin)
    peak_javelin_col = fits.Column(name='peak_javelin', format='E', array=peak_javelin)
    lag_javelin_col = fits.Column(name='lag_javelin', format='E', array=lag_javelin)
    lag_err_javelin_col = fits.Column(name='lag_err_javelin', format='2E', array=lag_err_javelin)
    frac_rejected_javelin_col = fits.Column(name='frac_rejected_javelin', format='E', array=frac_rejected_javelin)
    
    n0_pyroa_col = fits.Column(name='n0_pyroa', format='E', array=n0_pyroa)
    peak_bounds_pyroa_col = fits.Column(name='peak_bounds_pyroa', format='2E', array=peak_bounds_pyroa)
    peak_pyroa_col = fits.Column(name='peak_pyroa', format='E', array=peak_pyroa)
    lag_pyroa_col = fits.Column(name='lag_pyroa', format='E', array=lag_pyroa)
    lag_err_pyroa_col = fits.Column(name='lag_err_pyroa', format='2E', array=lag_err_pyroa)
    frac_rejected_pyroa_col = fits.Column(name='frac_rejected_pyroa', format='E', array=frac_rejected_pyroa)

    rmax_pyccf_col = fits.Column(name='rmax_pyccf', format='E', array=rmax_pyccf)    
    rmax_jav_col = fits.Column(name='rmax_javelin', format='E', array=rmax_jav)
    rmax_pyroa_col = fits.Column(name='rmax_pyroa', format='E', array=rmax_pyroa)

    cols = [name_col, k_col,
            n0_pyccf_col, peak_bounds_pyccf_col, peak_pyccf_col, lag_pyccf_col, lag_err_pyccf_col, frac_rejected_pyccf_col,
            n0_javelin_col, peak_bounds_javelin_col, peak_javelin_col, lag_javelin_col, lag_err_javelin_col, frac_rejected_javelin_col,
            n0_pyroa_col, peak_bounds_pyroa_col, peak_pyroa_col, lag_pyroa_col, lag_err_pyroa_col, frac_rejected_pyroa_col,
            rmax_pyccf_col, rmax_jav_col, rmax_pyroa_col]

    #---------------------------
    #Make table

    table_hdu = fits.BinTableHDU.from_columns(cols, header=hdr)
    table_hdu.writeto(output_fname, overwrite=True)

    return
