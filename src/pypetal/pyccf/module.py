import itertools
import multiprocessing as mp

import numpy as np

from pypetal.pyccf.plotting import plot_pyccf_results
from pypetal.pyccf.utils import get_pyccf_lags
from pypetal.utils import defaults
from pypetal.utils.petalio import write_data, print_subheader


#For multiprocessing
def mp_map(func, args, threads):
    n_inputs = len(args)

    if (threads > 1) & (n_inputs > 1):
        pool = mp.get_context('fork').Pool(threads)
        res = pool.starmap(func, args)

        pool.close()
        pool.join()

    else:
        res = list(itertools.starmap(func, args))

    return res



def pyccf_tot(cont_fname, line_fnames, line_names, output_dir,
              general_kwargs, kwargs):

    if line_fnames is str:
        line_fnames = [line_fnames]

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
    interp, nsim, mcmode, sigmode, thres, nbin = defaults.set_pyccf(kwargs)

    #-------------------------------------------
    #Run pyCCF for each line

    if verbose:

        if len(lag_bounds) > 2:
            lag_bounds_str = 'array'
        else:
            lag_bounds_str = lag_bounds
            
        print_dict = {
            'lag_bounds': lag_bounds_str,
            'interp': interp,
            'nsim': nsim,
            'mcmode': mcmode,
            'sigmode': sigmode,
            'thres': thres,
            'nbin': nbin
        }
        
        print_subheader('Running pyCCF', 35, print_dict)


    res_tot = []
    for i in range(len(line_fnames)):
        res = get_pyccf_lags( cont_fname, line_fnames[i], lag_bounds[i],
                                   interp=interp, nsim=nsim,
                                   mcmode=mcmode, sigmode=sigmode, thres=thres,
                                   threads=threads, verbose=verbose)

        #Write CCF to file
        dat_fname = output_dir + line_names[i+1] + r'/pyccf/' + line_names[i+1] + '_ccf.dat'
        header = '#Lag,CCF'
        write_data( [ res['CCF_lags'], res['CCF'] ], dat_fname, header )

        #Write CCCD, CCPD to file
        dat_fname = output_dir + line_names[i+1] + r'/pyccf/' + line_names[i+1] + '_ccf_dists.dat'
        header = '#CCCD,CCPD'
        write_data( [ res['CCCD_lags'], res['CCPD_lags'] ], dat_fname, header )


        res['name'] = line_names[i+1]

        x1, y1, yerr1 = np.loadtxt( cont_fname, delimiter=',', unpack=True, usecols=[0,1,2] )
        x2, y2, yerr2 = np.loadtxt( line_fnames[i], delimiter=',', unpack=True, usecols=[0,1,2] )

        plot_pyccf_results(x1, y1, yerr1, x2, y2, yerr2,
                                    res['CCF_lags'], res['CCF'],
                                    res['CCCD_lags'], res['CCPD_lags'],
                                    nbin=nbin, time_unit=time_unit, lc_unit=[lc_unit[0], lc_unit[i+1]],
                                    lc_names=[line_names[0], line_names[i+1]],
                                    fname = output_dir + line_names[i+1] + r'/pyccf/' + line_names[i+1] + '_ccf.pdf',
                                    show=plot)


        res_tot.append(res)


    return res_tot
