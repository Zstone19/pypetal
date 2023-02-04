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
