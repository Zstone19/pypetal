import itertools
import multiprocessing as mp
import os
import shutil
from functools import partial

import numpy as np
from astropy.table import Table

from pypetal.utils import defaults
import pypetal.pyzdcf.utils as utils
from pypetal.pyzdcf.plotting import plot_pyzdcf_results


#For multiprocessing
def mp_map(func, args, threads):
    n_inputs = len(args)

    if (threads > 1) & (n_inputs > 1):
        pool = mp.Pool(threads)
        res = pool.starmap(func, args)

        pool.close()
        pool.join()

    else:
        res = list(itertools.starmap(func, args))

    return res



def pyzdcf_tot(cont_fname, line_fnames, line_names, output_dir,
               general_kwargs, kwargs):


    if line_fnames is str:
        line_fnames = [line_fnames]
        line_names = [line_names]

    #--------------------------------------------------
    #Read general kwargs

    verbose = general_kwargs['verbose']
    plot = general_kwargs['plot']
    time_unit = general_kwargs['time_unit']
    lc_unit = general_kwargs['lc_unit']
    lag_bounds = general_kwargs['lag_bounds']
    threads = general_kwargs['threads']

    #--------------------------------------------------
    #Read kwargs

    nsim, minpts, uniform_sampling, omit_zero_lags, \
        sparse, prefix, run_plike, plike_dir = defaults.set_pyzdcf(kwargs,
                                                                   np.hstack([ [cont_fname], line_fnames ])
                                                                   )

    #-------------------------------------------
    if (run_plike) & (plike_dir is None):
        print('Error: plike_dir must be specified if run_plike=True')
        print('Skipping PLIKE')
        run_plike = False

    input_dir = os.path.dirname( os.path.realpath(cont_fname) )

    for i in range(len(line_fnames)):
        if input_dir != os.path.dirname( os.path.realpath(line_fnames[i]) ):
            print('ERROR: All light curve files must be in the same directory')
            return {}


    if verbose:

        txt_str = """
Running pyZDCF
----------------------
nsim: {}
minpts: {}
uniform_sampling: {}
omit_zero_lags: {}
sparse: {}
prefix: {}
run_plike: {}
plike_dir: {}
----------------------
        """.format( nsim, minpts, uniform_sampling, omit_zero_lags,
                    sparse, prefix, run_plike, plike_dir)

        print(txt_str)



    input_dir += r'/'
    cont_fname_short = os.path.basename(cont_fname)
    line_fnames_short = [ os.path.basename(x) for x in line_fnames ]


    pyzdcf_func = partial( utils.get_zdcf, num_MC=nsim, minpts=minpts,
                           uniform_sampling=uniform_sampling, omit_zero_lags=omit_zero_lags,
                           sparse=sparse, sep=',', verbose=False)

    arg1 = np.full( len(line_fnames), input_dir )
    arg2 = np.full( len(line_fnames), cont_fname_short )
    arg3 = line_fnames_short
    arg4 = [ output_dir + x + r'/pyzdcf/' for x in line_names[1:] ]
    arg5 = [ x + '_' + prefix for x in line_names[1:] ]

    args = list(zip(arg1, arg2, arg3, arg4, arg5))
    res_tot = mp_map(pyzdcf_func, args, threads)

    plike_tot = []
    for i in range(len(line_fnames)):

        x1, y1, yerr1 = np.loadtxt( cont_fname, delimiter=',', unpack=True, usecols=[0,1,2] )
        x2, y2, yerr2 = np.loadtxt( line_fnames[i], delimiter=',', unpack=True, usecols=[0,1,2] )

        plike_dict = None
        if run_plike:

            utils.run_plike( arg4[i] + arg5[i] + '.dcf', lag_bounds[i], plike_dir,
                            verbose=verbose)

            plike_fname = plike_dir + 'plike.out'
            output_fname = output_dir + line_names[i+1] + r'/pyzdcf/' + line_names[i+1] + '_plike.out'
            assert os.path.exists( os.path.dirname(output_fname) )

            shutil.move( plike_fname, output_fname )
            plike_dat = Table.read( output_fname,
                                    format='ascii',
                                    names=['num', 'lag', 'r', '-dr', '+dr', 'likelihood'])


            file = open(output_fname, 'r')
            output_str = list(file)[-3:]

            ml_lag = float( output_str[1].split()[7] )
            ml_lag_err_hi = np.abs( float( output_str[1].split()[8] )  )
            ml_lag_err_lo = np.abs( float( output_str[1].split()[9] )  )

            plike_dict = {
                'output': plike_dat,
                'ML_lag': ml_lag,
                'ML_lag_err_hi': ml_lag_err_hi,
                'ML_lag_err_lo': ml_lag_err_lo
            }

            plike_tot.append( plike_dict )

        plot_pyzdcf_results(x1, y1, yerr1, x2, y2, yerr2,
                                     res_tot[i], plike_dict,
                                     time_unit=time_unit, lc_unit=[lc_unit[0], lc_unit[i+1]],
                                     lc_names=[line_names[0], line_names[i+1]],
                                     fname=output_dir + line_names[i+1] + r'/pyzdcf/' + line_names[i+1] + '_zdcf.pdf',
                                     show=plot)

    return res_tot, plike_tot
