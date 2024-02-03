import os
import re

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.table import Table

from colorama import init, Fore
init(autoreset=True)

##############################################
# PRINTING

mcolor = Fore.BLUE
acolor = Fore.GREEN
wcolor = Fore.YELLOW
ecolor = Fore.RED

def print_hashstr():
    hashstr = "".join( ['#']*80 )
    print( mcolor + hashstr )
    return

def print_linestr(n):
    linestr = "".join( ['\u25AC']*n )
    print( mcolor + linestr )
    return    


def print_header(header):
    nhead = len(header)
    midpt = int(80/2)

    hstart = midpt - int(nhead/2) - 2
    hend = midpt + (nhead-int(nhead/2)) + 2

    nstart = hstart
    nend = 80 - hend
    assert nstart + 2 + nhead + 2 + nend == 80

    hashstr1 = "".join( ['#']*nstart )
    hashstr2 = "".join( ['#']*nend )

    print('')
    print_hashstr()
    print(mcolor + hashstr1, end='')
    print('  ', end='')
    print( acolor + header, end='')
    print('  ', end='')
    print(mcolor + hashstr2)
    print_hashstr()
    print('')
    return


def print_subheader(subheader, n, d):
    if n < len(subheader) + 5:
        n = len(subheader) + 5
    
    print_linestr(n)
    print(acolor + subheader)
    print_linestr(n)
    print_subheader_dict(d)
    print_linestr(n)
    print('')
    return


def print_subheader_dict(d):
    n = 15 + 5
    
    keys = list(d.keys())
    vals = list(d.values())
    for i in range(len(keys)):
        print(keys[i] + ':', end='')
        print( ''.join([' ']*(n-len(keys[i])-1)), end=''  )
        print('{}'.format(vals[i]))
        
    return

def print_warning(warning):
    print(wcolor + '*** ', end='')
    print(wcolor + warning)
    return

def print_error(error):
    print(ecolor + '*** ', end='')
    print(ecolor + error)
    return

##############################################

def make_directories(output_dir, fnames, line_names,
                     run_drw_rej, run_pyccf, run_pyzdcf, run_pyroa, run_mica2,
                     reject_data, together_pyroa, together_mica2):

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


    if run_mica2:
        if together_mica2:
            os.makedirs( output_dir + 'mica2', exist_ok=True )
        else:
            for i in range(len(fnames)-1):
                os.makedirs( output_dir + line_names[i+1] + '/mica2', exist_ok=True )


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



def write_weighting_summary(fname, res, run_pyccf, run_javelin, run_pyroa, run_mica2):

    run_arr = [run_pyccf, run_javelin, run_pyroa, run_mica2]
    names = ['pyccf', 'javelin', 'pyroa', 'mica2']
    mod_ncols = 6

    colnames = ['n0', 'peak_bounds', 'peak', 'lag', 'lag_err', 'frac_rejected']
    colnames_tot = []
    for name in names:
        for col in colnames:
            colnames_tot.append( col + '_' + name )

    weighting_dat = []
    weighting_dat.append( res['k'] )
    for i, run in enumerate(run_arr):
        if run:
            for i in range( i*mod_ncols, (i+1)*mod_ncols ):
                weighting_dat.append( res[ colnames_tot[i] ] )

        else:
            for i in range(mod_ncols):
                if i in [1, 4]:
                    weighting_dat.append( [np.nan, np.nan] )
                else:
                    weighting_dat.append( np.nan )

    if run_pyccf:
        weighting_dat.append( res['rmax_pyccf'] )
        
        for i in range(1, len(run_arr)):
            if run_arr[i]:
                weighting_dat.append( res['rmax_' + names[i]] )
            else:
                weighting_dat.append( np.nan )
        
    else:
        weighting_dat.append(np.nan)
        weighting_dat.append(np.nan)
        weighting_dat.append(np.nan)
        weighting_dat.append(np.nan)


    tot_fits_cols = []

    fits_colnames = []
    fits_colnames.append('k')
    for i in range(len(colnames_tot)):
        fits_colnames.append( colnames_tot[i] )

    fits_colnames.append('rmax_pyccf')
    fits_colnames.append('rmax_javelin')
    fits_colnames.append('rmax_pyroa')
    fits_colnames.append('rmax_mica2')



    col_fmts = []
    col_fmts.append('E')
    for _ in range(len(run_arr)):
        for j in range(mod_ncols):
            if j in [1,4]:
                col_fmts.append('2E')
            else:
                col_fmts.append('E')

    for _ in range(len(run_arr)):
        col_fmts.append('E')


    for i in range(len(weighting_dat)):
        fits_col_i = fits.Column(name=fits_colnames[i], format=col_fmts[i], array=[ weighting_dat[i] ])
        tot_fits_cols.append(fits_col_i)

    table = fits.BinTableHDU.from_columns(tot_fits_cols)
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

    n0_mica2 = np.zeros(len(filenames))
    peak_bounds_mica2 = np.zeros( (len(filenames),2) )
    peak_mica2 = np.zeros(len(filenames))
    lag_mica2 = np.zeros(len(filenames))
    lag_err_mica2 = np.zeros( (len(filenames),2) )
    frac_rejected_mica2 = np.zeros(len(filenames))

    rmax_pyccf = np.zeros(len(filenames))
    rmax_jav = np.zeros(len(filenames))
    rmax_pyroa = np.zeros(len(filenames))
    rmax_mica2 = np.zeros(len(filenames))


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
        
        n0_mica2[i] = table['n0_mica2'][0]
        peak_bounds_mica2[i,:] = table['peak_bounds_mica2'][0]
        peak_mica2[i] = table['peak_mica2'][0]
        lag_mica2[i] = table['lag_mica2'][0]
        lag_err_mica2[i,:] = table['lag_err_mica2'][0]
        frac_rejected_mica2[i] = table['frac_rejected_mica2'][0]
        

        rmax_pyccf[i] = table['rmax_pyccf'][0]
        rmax_jav[i] = table['rmax_javelin'][0]
        rmax_pyroa[i] = table['rmax_pyroa'][0]
        rmax_mica2[i] = table['rmax_mica2'][0]

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

    n0_mica2_col = fits.Column(name='n0_mica2', format='E', array=n0_mica2)
    peak_bounds_mica2_col = fits.Column(name='peak_bounds_mica2', format='2E', array=peak_bounds_mica2)
    peak_mica2_col = fits.Column(name='peak_mica2', format='E', array=peak_mica2)
    lag_mica2_col = fits.Column(name='lag_mica2', format='E', array=lag_mica2)
    lag_err_mica2_col = fits.Column(name='lag_err_mica2', format='2E', array=lag_err_mica2)
    frac_rejected_mica2_col = fits.Column(name='frac_rejected_mica2', format='E', array=frac_rejected_mica2)

    rmax_pyccf_col = fits.Column(name='rmax_pyccf', format='E', array=rmax_pyccf)
    rmax_jav_col = fits.Column(name='rmax_javelin', format='E', array=rmax_jav)
    rmax_pyroa_col = fits.Column(name='rmax_pyroa', format='E', array=rmax_pyroa)
    rmax_mica2_col = fits.Column(name='rmax_mica2', format='E', array=rmax_mica2)

    cols = [name_col, k_col,
            n0_pyccf_col, peak_bounds_pyccf_col, peak_pyccf_col, lag_pyccf_col, lag_err_pyccf_col, frac_rejected_pyccf_col,
            n0_javelin_col, peak_bounds_javelin_col, peak_javelin_col, lag_javelin_col, lag_err_javelin_col, frac_rejected_javelin_col,
            n0_pyroa_col, peak_bounds_pyroa_col, peak_pyroa_col, lag_pyroa_col, lag_err_pyroa_col, frac_rejected_pyroa_col,
            n0_mica2_col, peak_bounds_mica2_col, peak_mica2_col, lag_mica2_col, lag_err_mica2_col, frac_rejected_mica2_col,
            rmax_pyccf_col, rmax_jav_col, rmax_pyroa_col, rmax_mica2_col]

    #---------------------------
    #Make table

    table_hdu = fits.BinTableHDU.from_columns(cols, header=hdr)
    table_hdu.writeto(output_fname, overwrite=True)

    return
