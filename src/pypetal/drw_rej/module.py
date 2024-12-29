import itertools
import multiprocessing as mp
from functools import partial

import astropy.units as u
import numpy as np

from pypetal.drw_rej.utils import drw_flag
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



def str2unit(x):

    if (x == 'Arbitrary Units'):
        unit = u.dimensionless_unscaled
    else:
        unit = u.Unit(x)

    return unit


def drw_rej_tot(cont_fname, line_fnames, line_names, output_dir,
                general_kwargs, kwargs):


    #--------------------------------------------------
    #Read general kwargs

    verbose = general_kwargs['verbose']
    plot = general_kwargs['plot']
    time_unit = general_kwargs['time_unit']
    lc_unit = general_kwargs['lc_unit']
    threads = general_kwargs['threads']

    #--------------------------------------------------
    #Read kwargs

    jitter, nsig, nwalker, nburn, nchain, clip, \
        reject_data, use_for_javelin, solver = defaults.set_drw_rej(kwargs,
                                                            np.hstack([ [cont_fname], line_fnames ])
                                                            )

    #--------------------------------------------------
    #Get units

    time_unit = str2unit(time_unit)

    for i in range(len(lc_unit)):
        lc_unit[i] = str2unit(lc_unit[i])

    #--------------------------------------------------
    if verbose:

        if len(clip) > 2:
            clip_str = 'array'
        else:
            clip_str = clip

        print_dict = {
            'jitter': jitter,
            'nsig': nsig,
            'nwalker': nwalker,
            'nburn': nburn,
            'nchain': nchain,
            'clip': clip_str,
            'reject_data': reject_data,
            'use_for_javelin': use_for_javelin,
            'solver': solver
        }

        print_subheader('Performing DRW Rejection', 35, print_dict)


    sigmas = []
    taus = []

    if jitter:
        jitters = []

    rejecting = np.any(reject_data)
    fnames = np.hstack([ [cont_fname], line_fnames ])


    drw_rej_func = partial( drw_flag, nwalkers=nwalker, nburn=nburn, nsamp=nchain,
                            nsig=nsig, jitter=jitter, solver=solver, plot=plot)



    arg1 = []
    arg2 = []
    arg3 = []
    arg4 = line_names
    arg5 = []

    for i in range(len(fnames)):
        if reject_data[i]:
            x, y, yerr = np.loadtxt( fnames[i], delimiter=',', unpack=True, usecols=[0,1,2] )

            arg1.append( x*time_unit )
            arg2.append( y*lc_unit[i] )
            arg3.append( yerr*lc_unit[i] )
            arg5.append( output_dir + line_names[i] + '/drw_rej/' + line_names[i] + '_drw_fit.pdf' )


    args = list(zip(arg1, arg2, arg3, arg4, arg5))
    res_tot = mp_map(drw_rej_func, args, threads)

    masks_tot = []
    n = 0
    for i in range(len(fnames)):

        if reject_data[i]:
            res = res_tot[n]

            sigmas.append( res['sigma'])
            taus.append( res['tau'])

            if jitter:
                jitters.append( res['jitter'] )

            mask = res['mask']

            #Write mask
            dat_fname = output_dir + line_names[i] + '/drw_rej/' + line_names[i] + '_mask.dat'
            write_data(mask, dat_fname)

            #Write LC fit
            dat_fname = output_dir + line_names[i] + '/drw_rej/' + line_names[i] + '_drw_fit.dat'
            header = '#time,light curve,error'
            write_data( [ res['fit_x'], res['fit_y'], res['fit_err'] ], dat_fname, header )

            #Write the MCMC chains
            dat_fname = output_dir + line_names[i] + '/drw_rej/' + line_names[i] + '_chain.dat'

            if jitter:
                header = '#sigma,tau,jitter'
                write_data( [ res['sigma'], res['tau'], res['jitter'] ], dat_fname, header )
            else:
                header = '#sigma,tau'
                write_data( [ res['sigma'], res['tau'] ], dat_fname, header )


            x_new = x[~mask]
            y_new = y[~mask]
            yerr_new = yerr[~mask]

            dat_fname = output_dir + 'processed_lcs/' + line_names[i] + '_data.dat'
            write_data( [x_new, y_new, yerr_new], dat_fname )


            n += 1

        else:
            x, y, yerr = np.loadtxt( fnames[i], delimiter=',', unpack=True, usecols=[0,1,2] )

            mask = np.full( len(x), False )
            taus.append(None)
            sigmas.append(None)

            if jitter:
                jitters.append(None)

            if rejecting:
                dat_fname = output_dir + 'processed_lcs/' + line_names[i] + '_data.dat'
                write_data([x, y, yerr], dat_fname)


        masks_tot.append(mask)

    output = {
        'masks': masks_tot,
        'reject_data': reject_data,
        'taus': taus,
        'sigmas': sigmas
    }

    if jitter:
        output['jitters'] = jitters

    return output
