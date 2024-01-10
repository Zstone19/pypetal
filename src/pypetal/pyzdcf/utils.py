import os
import subprocess

import numpy as np
from pyzdcf import pyzdcf
from pypetal.utils.petalio import print_error


def run_plike(dcf_fname, lag_bounds, plike_dir, verbose=False):

    """
    Runs the PLIKE algorithm to compute the maximum likelihood peak of the ZDCF.
    Must have PLIKE (https://ui.adsabs.harvard.edu/abs/2013arXiv1302.1508A/abstract).
    Will output a file containing the PLIKE results ('plike.out') in ``plike_dir``.


    Parameters
    ----------

    plike_dir : str
            Path to the directory with the PLIKE executable.

    lag_bounds : (2,) array_like
            Lower and upper bounds of lags to search for ML peak.

    dcf_name :str
            Path to the ZDCF file, usually ouptput by pyZDCF.

    verbose : bool, optional
            If True, will read output of PLIKE. Default is False.


    Returns
    -------

    """

    #Make sure plike dir exists
    assert os.path.exists( plike_dir )
    plike_dir = os.path.abspath(plike_dir) + r'/'

    cwd = os.getcwd()
    os.chdir(plike_dir)

    #Delete old plike.out if it exists
    if os.path.exists( plike_dir + 'plike.out' ):
        os.remove( plike_dir + 'plike.out' )

    #Make sure dcf file exists
    assert os.path.exists(dcf_fname)
    dcf_fname = os.path.abspath(dcf_fname)


    #Make sure there are two lag bounds
    if len(lag_bounds) != 2:
        print_error('ERROR: Must provide two lag bounds.')
        print_error('Lag bounds: ', lag_bounds)
        raise ValueError('Must provide two lag bounds.')


    #Make file with the arguments to pass to plike
    with open('args.txt', 'w') as f:
        f.write(dcf_fname + '\n')
        f.write(str(lag_bounds[0]) + '\n')
        f.write(str(lag_bounds[1]))


    if verbose:
        print('Executing PLIKE')

    exec_str = './plike < args.txt'
    res = subprocess.Popen(exec_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output, error = res.communicate()
    if res.returncode != 0:
        raise Exception("File handling failed %d %s %s" % (res.returncode, output, error))


    #Delete args.txt
    os.remove('args.txt')


    os.chdir(cwd)

    return


def get_zdcf(input_dir, fname1, fname2, out_dir, prefix='zdcf', num_MC=500, minpts=0,
             uniform_sampling=False, autocf=False, omit_zero_lags=True,
             sparse='auto', sep=',', verbose=False):

    """
    Runs the pyZDCF algorithm to compute the Z-Transformed Discrete Correlation Function (ZDCF).
    For more information on the algorithm, see pyZDCF (https://pyzdcf.readthedocs.io/).
    The algorithm will take in two light curves and output a file containing the ZDCF ('zdcf.dcf')
    in the specified output directory.


    Parameters
    ----------

    input_dir : str
            Path to the directory containing the light curves.

    fname1 : str
            Name of the first light curve file.

    fname2 : str
            Name of the second light curve file.

    out_dir : str
            Path to the directory to place the output ZDCF file.

    num_MC : float, optional
            The number of Monte Carlo simulations to run. Default is 500.

    minpts : int, optional
            The minimum number of points to use in each bin when computing the ZDCF.
            Must be larger than 11. If set to 0, it will be set to 11. Default is 0.

    uniform_sampling: bool, optional
            If True, the light curves will be assumed to be uniformly sampled.
            Default is ``False``.

    autocf : bool, optional
            If True, the auto-correlation function for the first light curve will be computed.
            If False, the ZDCF will be computed between the light curves.
            Default is ``False``.

    omit_zero_lags : bool, optional
            If True, will omit the points with zero lags when computing the ZDCF.
            Default is ``True``.

    sparse : (bool, str), optional
            Determines whether to use a sparse matrix implementation for reduced RAM usage.
            This feature is suitable for longer light curves (> 3000 data points). If True, will
            use sparse matrix implementation. If set to 'auto', will use sparse matrix implementation
            if there are more than 3000 data points per light curve. Default is 'auto'.

    sep : str, optional
            The delimiter used in the light curve files. Default is ',' for CSV files.

    verbose : bool, optional
            If ``True``, will output progress of pyZDCF. Default is ``False``.



    Returns
    -------

    dcf_df : pandas.DataFrame
        A pandas DataFrame object containing the ZDCF and its errors. The
        columns within the DataFrame match the order of columns within the
        output ZDCF file 'zdcf.dcf'.

    """


    #Files need to be csv
    params = dict(
        autocf = autocf,
        prefix = prefix,
        uniform_sampling = uniform_sampling,
        omit_zero_lags = True,
        minpts = 0,
        num_MC = num_MC,
        lc1_name = fname1,
        lc2_name = fname2
    )

    dcf_df = pyzdcf(input_dir=input_dir, output_dir=out_dir, intr=False,
                    verbose=verbose, parameters=params, sep=sep, sparse=sparse )

    return dcf_df
