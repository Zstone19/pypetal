import glob
import os

import numpy as np
import pandas as pd
from astropy.table import Table
from javelin.zylc import get_data


def get_line_names(main_dir):
    dirs = glob.glob(main_dir + '*/')

    if main_dir + 'light_curves/' in dirs:
        dirs.remove(main_dir + 'light_curves/')

    if main_dir + 'processed_lcs/' in dirs:
        dirs.remove(main_dir + 'processed_lcs/')

    if main_dir + 'javelin/' in dirs:
        if len(dirs) > 2:
            dirs.remove(main_dir + 'javelin/' )
            line_names = [ os.path.basename(os.path.dirname(x)) for x in dirs ]
        else:
            fits_files = glob.glob(main_dir + 'javelin/*_lc_fits.dat')
            line_names = [ os.path.basename(x[:-12]) for x in fits_files ]

    else:
        line_names = [ os.path.basename(os.path.dirname(x)) for x in dirs ]

    return line_names


def str2bool(string):
    if string == 'False':
        return False
    elif string == 'True':
        return True


def get_modules(main_dir):

    line_names = get_line_names(main_dir)
    line_dirs = np.array([ main_dir + x + r'/' for x in line_names ])

    ###########################################################
    #DRW Rejection
    has_drw = np.zeros( len(line_names), dtype=bool )
    for i, dir_i in enumerate(line_dirs):
        subdirs = glob.glob(dir_i + '*/')
        if dir_i + 'drw_rej/' in subdirs:
            has_drw[i] = True

    run_drw_rej = np.any( has_drw )

    ###########################################################
    #Detrend

    has_detrend = np.zeros( len(line_names), dtype=bool )
    for i, dir_i in enumerate(line_dirs):
        subdirs = glob.glob(dir_i + '*.pdf')
        if dir_i + 'detrend.pdf' in subdirs:
            has_detrend[i] = True

    run_detrend = np.any( has_detrend )


    #Get number of lines (not counting continuum)
    n_lnc = len(line_names) - 1

    ###########################################################
    #pyCCF
    has_pyccf = np.zeros( len(line_names), dtype=bool )
    for i, dir_i in enumerate(line_dirs):
        subdirs = glob.glob(dir_i + '*/')
        if dir_i + 'pyccf/' in subdirs:
            has_pyccf[i] = True


    n_pyccf = len( np.argwhere(has_pyccf).T[0] )
    if ( n_pyccf != n_lnc ) & ( n_pyccf > 0 ):
        print('pyCCF was not completed for all lines, so will assume run_pyccf=False')
        run_pyccf = False
    else:
        run_pyccf = np.any( has_pyccf )


    ###########################################################
    #pyZDCF
    has_pyzdcf = np.zeros( len(line_names), dtype=bool )
    for i, dir_i in enumerate(line_dirs):
        subdirs = glob.glob(dir_i + '*/')
        if dir_i + 'pyzdcf/' in subdirs:
            has_pyzdcf[i] = True


    n_pyzdcf = len( np.argwhere(has_pyzdcf).T[0] )
    if ( n_pyzdcf != n_lnc ) & ( n_pyzdcf > 0 ):
        print('pyZDCF was not completed for all lines, so will assume run_pyzdcf=False')
        run_pyzdcf = False
    else:
        run_pyzdcf = np.any( has_pyzdcf )



    ###########################################################
    #JAVELIN
    run_javelin = False

    dirs_tot = glob.glob(main_dir + '*/')
    if main_dir + 'javelin/' in dirs_tot:
        run_javelin = True

    else:
        has_javelin = np.zeros( len(line_names), dtype=bool )
        for i,dir_i in enumerate(line_dirs):
            subdirs = glob.glob(dir_i + '*/')

            if dir_i + 'javelin/' in subdirs:
                has_javelin[i] = True

        n_javelin = len( np.argwhere(has_javelin).T[0] )
        if ( n_javelin != n_lnc ) & ( n_javelin > 0 ):
            print('JAVELIN was not completed for all lines, so will assume run_javelin=False')
        else:
            run_javelin = np.any( has_javelin )



    ###########################################################
    #Weighting
    has_weighting = np.zeros( len(line_names), dtype=bool )
    for i, dir_i in enumerate(line_dirs):
        subdirs = glob.glob(dir_i + '*/')

        if dir_i + 'weights/' in subdirs:
            has_weighting[i] = True

    n_weighting = len( np.argwhere(has_weighting).T[0] )
    if ( n_weighting != n_lnc ) & ( n_weighting > 0 ):
        print('Weighting was not completed for all lines, so will assume run_weighting=False')
        run_weighting = False
    else:
        run_weighting = np.any( has_weighting )



    return run_drw_rej, run_pyccf, run_pyzdcf, run_javelin, run_weighting



def get_javelin_together(main_dir):


    dirs_tot = glob.glob( main_dir + '*/')

    dirs_tot.remove( main_dir + 'light_curves/')
    if main_dir + 'processed_lcs/' in dirs_tot:
        dirs_tot.remove( main_dir + 'processed_lcs/')

    if main_dir + 'javelin/' in dirs_tot:
        together = True
    else:
        together = False

    return together



def get_cont_name(main_dir):

    run_drw_rej, run_pyccf, run_pyzdcf, run_javelin, _ = get_modules(main_dir)

    if run_pyccf:

        dirs_tot = glob.glob( main_dir + '*/')

        dirs_tot.remove( main_dir + 'light_curves/')
        if main_dir + 'processed_lcs/' in dirs_tot:
            dirs_tot.remove( main_dir + 'processed_lcs/')

        has_pyccf = np.zeros( len(dirs_tot), dtype=bool )
        for i, dir_i in enumerate(dirs_tot):
            subdirs = glob.glob(dir_i + '*/')
            if dir_i + 'pyccf/' in subdirs:
                has_pyccf[i] = True

        ind = np.argwhere( ~has_pyccf ).T[0][0]

        cont_name = os.path.basename( os.path.dirname( dirs_tot[ind] ) )


    elif run_pyzdcf:

        dirs_tot = glob.glob( main_dir + '*/')

        dirs_tot.remove( main_dir + 'light_curves/')
        if main_dir + 'processed_lcs/' in dirs_tot:
            dirs_tot.remove( main_dir + 'processed_lcs/')

        has_pyzdcf = np.zeros( len(dirs_tot), dtype=bool )
        for i, dir_i in enumerate(dirs_tot):
            subdirs = glob.glob(dir_i + '*/')
            if dir_i + 'pyzdcf/' in subdirs:
                has_pyzdcf[i] = True

        ind = np.argwhere( ~has_pyzdcf ).T[0][0]

        cont_name = os.path.basename(  os.path.dirname( dirs_tot[ind] )  )


    elif run_javelin:
        together = get_javelin_together(main_dir)

        dirs_tot = glob.glob( main_dir + '*/')

        processed = False
        if main_dir + 'processed_lcs/' in dirs_tot:
            dirs_tot.remove( main_dir + 'processed_lcs/')
            processed = True

        if together:
            dirs_tot.remove( main_dir + 'javelin/' )

            lc_files = glob.glob( main_dir + 'light_curves/*.dat' )
            cont_dat = get_data( main_dir + 'javelin/cont_lcfile.dat' )

            xcont = cont_dat.jlist[0]
            ycont = cont_dat.mlist[0] + cont_dat.blist[0]
            yerrcont = cont_dat.elist[0]
            for f in lc_files:
                xline, yline, yerrline = np.loadtxt(f, unpack=True, delimiter=',', usecols=[0,1,2])

                if np.all( xline == xcont ) & np.all( yline == ycont ) & np.all( yerrline == yerrcont ):
                    cont_name = os.path.basename(f)[:-4]
                    break


        else:

            has_javelin = np.zeros( len(dirs_tot), dtype=bool )
            for i, dir_i in enumerate(dirs_tot):
                subdirs = glob.glob(dir_i + '*/')

                if dir_i + 'javelin/' in subdirs:
                    has_javelin[i] = True

            ind = np.argwhere( ~has_javelin ).T[0][0]
            cont_name = os.path.basename( os.path.dirname( dirs_tot[ind] ) )



    elif run_drw_rej:
        line_names, reject_data = get_reject_data(main_dir, ordered=False)

        if np.all(~reject_data):
            cont_name = line_names[0]
        else:
            possible_names = np.array(line_names)[reject_data]

            cont_name = possible_names[0]

    return cont_name



def get_ordered_line_names(main_dir):
    line_names = get_line_names(main_dir)
    cont_name = get_cont_name(main_dir)

    line_names.remove(cont_name)

    line_names_tot = np.hstack([ [cont_name], line_names ])

    return line_names_tot


#######################################################################
#                            DRW REJECTION
#######################################################################


def get_reject_data(main_dir, ordered=True):

    if ordered:
        line_names = get_ordered_line_names(main_dir)
    else:
        line_names = get_line_names(main_dir)

    reject_data = np.zeros( len(line_names), dtype=bool )
    for i, name in enumerate(line_names):
        dir_i = main_dir + name + r'/'

        subdirs = glob.glob(dir_i + '*/')
        if dir_i + 'drw_rej/' in subdirs:
            reject_data[i] = True

    return line_names, reject_data


def load_drw_rej(dir_loc):

    line_names = get_ordered_line_names(dir_loc)
    cont_name = get_cont_name(dir_loc)

    line_dirs = np.array([ dir_loc + x + r'/' for x in line_names ])
    cont_dir = dir_loc + cont_name


    res_dict = {}
    res_dict['taus'] = []
    res_dict['sigmas'] = []
    res_dict['jitters'] = []

    res_dict['fit_x'] = []
    res_dict['fit_y'] = []
    res_dict['fit_err'] = []

    res_dict['masks'] = []
    res_dict['reject_data'] = get_reject_data(dir_loc)[1]
    res_dict['names'] = []

    for i, dir_i in enumerate(line_dirs):

        if res_dict['reject_data'][i]:
            chain_dat = Table.read( dir_i + 'drw_rej/' + line_names[i] + '_chain.dat',
                                   format='ascii.csv')

            res_dict['taus'].append( np.array(chain_dat['tau']) )
            res_dict['sigmas'].append( np.array(chain_dat['#sigma']) )

            if 'jitter' in chain_dat.colnames:
                res_dict['jitters'].append( np.array(chain_dat['jitter']) )
            else:
                res_dict['jitters'].append(None)


            x, y, yerr = np.loadtxt( dir_i + 'drw_rej/' + line_names[i] + '_drw_fit.dat', unpack=True, delimiter=',', usecols=[0,1,2] )
            res_dict['fit_x'].append(x)
            res_dict['fit_y'].append(y)
            res_dict['fit_err'].append(yerr)

            mask = np.loadtxt( dir_i + 'drw_rej/' + line_names[i] + '_mask.dat', dtype=str)
            mask = list( map(str2bool, mask) )
            res_dict['masks'].append(mask)

            res_dict['names'].append(line_names[i])

        else:
            res_dict['taus'].append(None)
            res_dict['sigmas'].append(None)
            res_dict['jitters'].append(None)

            res_dict['fit_x'].append(None)
            res_dict['fit_y'].append(None)
            res_dict['fit_err'].append(None)


            x, y, yerr = np.loadtxt( dir_loc + 'light_curves/' + line_names[i] + '.dat', unpack=True, delimiter=',', usecols=[0,1,2] )
            res_dict['masks'].append( np.zeros( len(x), dtype=bool ) )
            res_dict['names'].append(line_names[i])

    return res_dict


#######################################################################
#                              PYCCF
#######################################################################


def load_pyccf(dir_loc):

    line_names = get_ordered_line_names(dir_loc)
    dirs_tot = np.array([ dir_loc + x + r'/' for x in line_names ])

    line_dirs = dirs_tot[1:]
    cont_dir = dirs_tot[0]

    res_dict_tot = []
    for i, dir_i in enumerate(line_dirs):
        pyccf_dir = dir_i + 'pyccf/'
        dict_i = {}

        cccd_dist, ccpd_dist = np.loadtxt( pyccf_dir + line_names[i+1] + '_ccf_dists.dat', delimiter=',', unpack=True, usecols=[0,1])
        dict_i['CCCD_lags'] = cccd_dist
        dict_i['CCPD_lags'] = ccpd_dist

        lags, ccf = np.loadtxt( pyccf_dir + line_names[i+1] + '_ccf.dat', delimiter=',', unpack=True, usecols=[0,1])
        dict_i['CCF_lags'] = lags
        dict_i['CCF'] = ccf

        dict_i['name'] = line_names[i+1]

        res_dict_tot.append(dict_i)

    return res_dict_tot


#######################################################################
#                             PYZDCF
#######################################################################

def get_run_plike(line_dir):
    pyzdcf_dir = line_dir + 'pyzdcf/'

    line_name = os.path.basename( os.path.dirname(line_dir) )

    if pyzdcf_dir + line_name + '_plike.out' in glob.glob( pyzdcf_dir + '*' ):
        return True
    else:
        return False


def load_pyzdcf(dir_loc):

    line_names = get_ordered_line_names(dir_loc)
    dirs_tot = np.array([ dir_loc + x + r'/' for x in line_names ])

    line_dirs = dirs_tot[1:]
    cont_dir = dirs_tot[0]

    pyzdcf_dict_tot = []
    plike_dict_tot = []
    for i, dir_i in enumerate(line_dirs):
        pyzdcf_dir = dir_i + 'pyzdcf/'
        dict_i = {}

        dcf_file = glob.glob(pyzdcf_dir + '*.dcf')[0]
        zdcf_df = Table.read(dcf_file, format='ascii',
                             names=['tau', '-sig(tau)', '+sig(tau)', 'dcf', '-err(dcf)', '+err(dcf)', '#bin']).to_pandas()


        dict_i['output'] = zdcf_df
        dict_i['name'] = line_names[i+1]

        pyzdcf_dict_tot.append(dict_i)


        run_plike = get_run_plike(dir_i)
        if run_plike:
            dict_i = {}

            plike_file = glob.glob(pyzdcf_dir + '*.out')[0]
            plike_df = Table.read(plike_file, format='ascii',
                                  names=['num', 'lag', 'r', '-dr', '+dr', 'likelihood'])


            with open(plike_file, 'r') as file:
                output_str = list(file)[-3:]

                ml_lag = float( output_str[1].split()[7] )
                ml_lag_err_hi = np.abs( float( output_str[1].split()[8] )  )
                ml_lag_err_lo = np.abs( float( output_str[1].split()[9] )  )


            dict_i = {
                'output': plike_df,
                'ML_lag': ml_lag,
                'ML_lag_err_hi': ml_lag_err_hi,
                'ML_lag_err_lo': ml_lag_err_lo,
                'name': line_names[i+1]
            }
            plike_dict_tot.append(dict_i)

        else:
            plike_dict_tot.append({})

    return pyzdcf_dict_tot, plike_dict_tot


#######################################################################
#                             JAVELIN
#######################################################################

def load_javelin(dir_loc):

    line_names = get_ordered_line_names(dir_loc)
    dirs_tot = np.array([ dir_loc + x + r'/' for x in line_names ])
    together = get_javelin_together(dir_loc)

    line_dirs = dirs_tot[1:]
    cont_dir = dirs_tot[0]


    if together:
        res_dict = {}

        javelin_dir = dir_loc + 'javelin/'

        chains = np.loadtxt( javelin_dir + 'chain_rmap.txt', unpack=True )
        res_dict['tau'] = np.exp(chains[1])
        res_dict['sigma'] = np.exp(chains[0])
        res_dict['tophat_params'] = chains[2:]

        x, y, yerr = np.loadtxt(javelin_dir + line_names[0] + '_lc_fits.dat', unpack=True, delimiter=',', usecols=[0,1,2])
        res_dict['cont_fit_x'] = x
        res_dict['cont_fit_y'] = y
        res_dict['cont_fit_yerr'] = yerr

        for i in range(1, len(line_names)):
            x, y, yerr = np.loadtxt(javelin_dir + line_names[i] + '_lc_fits.dat', unpack=True, delimiter=',', usecols=[0,1,2])
            res_dict[line_names[i] + '_fit_x'] = x
            res_dict[line_names[i] + '_fit_y'] = y
            res_dict[line_names[i] + '_fit_yerr'] = yerr

        return res_dict


    else:

        res_dict_tot = []
        for i, dir_i in enumerate(line_dirs):
            javelin_dir = dir_i + 'javelin/'
            dict_i = {}


            chains = np.loadtxt( javelin_dir + 'chain_rmap.txt', unpack=True )
            dict_i['tau'] = np.exp(chains[1])
            dict_i['sigma'] = np.exp(chains[0])
            dict_i['tophat_params'] = chains[2:]

            x, y, yerr = np.loadtxt(javelin_dir + line_names[0] + '_lc_fits.dat', unpack=True, delimiter=',', usecols=[0,1,2])
            dict_i['cont_fit_x'] = x
            dict_i['cont_fit_y'] = y
            dict_i['cont_fit_yerr'] = yerr

            x, y, yerr = np.loadtxt(javelin_dir + line_names[i+1] + '_lc_fits.dat', unpack=True, delimiter=',', usecols=[0,1,2])
            dict_i['fit_x'] = x
            dict_i['fit_y'] = y
            dict_i['fit_yerr'] = yerr

            dict_i['name'] = line_names[i+1]
            res_dict_tot.append(dict_i)


        return res_dict_tot


#######################################################################
#                             Weighting
#######################################################################

def read_weighting_summary(fname):

    res_dict = {}
    with open(fname) as f:

        for i, line in enumerate( f.readlines() ):
            line = line.strip('\n')

            #Total
            if i == 0:
                res_dict['k'] = float(line.strip('k = '))


            #pyCCF
            if i == 4:
                res_dict['pyccf_n0'] = int(line.strip('n0 = '))

            if i == 5:
                line = line.strip('peak bounds = [ ]')
                line = line.split(',')

                bounds = list(map(float, line))
                res_dict['pyccf_lag_bounds'] = bounds

            if i == 6:
                res_dict['pyccf_peak'] = float(line.strip('peak = '))

            if i == 7:
                res_dict['pyccf_lag_value'] = float(line.strip('lag value = '))

            if i == 8:
                line = line.strip('lag uncertainty = [ ]')
                line = line.split(',')

                bounds = list(map(float, line))
                res_dict['pyccf_lag_uncertainty'] = bounds

            if i == 9:
                res_dict['pyccf_frac_rejected'] = float(line.strip('fraction rejected = '))



            #JAVELIN
            if i == 13:
                res_dict['javelin_n0'] = int(line.strip('n0 = '))

            if i == 14:
                line = line.strip('peak bounds = [ ]')
                line = line.split(',')

                bounds = list(map(float, line))
                res_dict['javelin_lag_bounds'] = bounds

            if i == 15:
                res_dict['javelin_peak'] = float(line.strip('peak = '))

            if i == 16:
                res_dict['javelin_lag_value'] = float(line.strip('lag value = '))

            if i == 17:
                line = line.strip('lag uncertainty = [ ]')
                line = line.split(',')

                bounds = list(map(float, line))
                res_dict['javelin_lag_uncertainty'] = bounds

            if i == 18:
                res_dict['javelin_frac_rejected'] = float(line.strip('fraction rejected = '))



            #Total
            if i == 22:
                res_dict['rmax'] = float(line.strip('rmax = '))


    return res_dict




def load_weighting(main_dir):

    line_names = get_ordered_line_names(main_dir)
    dirs_tot = np.array([ main_dir + x + r'/' for x in line_names ])

    line_dirs = dirs_tot[1:]
    cont_dir = dirs_tot[0]



    output = {
        'pyccf': {},
        'javelin': {},
        'rmax': [],
        'names': []
    }

    output['pyccf']['centroid'] = []
    output['pyccf']['bounds'] = []
    output['pyccf']['acf'] = []
    output['pyccf']['lags'] = []
    output['pyccf']['weight_dist'] = []
    output['pyccf']['smoothed_dist'] = []
    output['pyccf']['ntau'] = []
    output['pyccf']['downsampled_CCCD'] = []
    output['pyccf']['frac_rejected'] = []


    output['javelin']['tophat_lag'] = []
    output['javelin']['bounds'] = []
    output['javelin']['acf'] = []
    output['javelin']['lags'] = []
    output['javelin']['weight_dist'] = []
    output['javelin']['smoothed_dist'] = []
    output['javelin']['ntau'] = []
    output['javelin']['downsampled_lag_dist'] = []
    output['javelin']['frac_rejected'] = []


    res_dicts_tot = []
    for i, dir_i in enumerate(line_dirs):
        weight_dir = dir_i + 'weights/'
        summary_dict = read_weighting_summary(weight_dir + 'weight_summary.txt')


        javelin_files = glob.glob( weight_dir + 'javelin*' )
        pyccf_files = glob.glob( weight_dir + 'pyccf*' )

        jweight = False
        pweight = False

        if len(javelin_files) > 0:
            jweight = True
        if len(pyccf_files) > 0:
            pweight = True


        if jweight:
            lag_err_lo = summary_dict['javelin_lag_uncertainty'][0]
            lag_err_hi = summary_dict['javelin_lag_uncertainty'][1]
            med_lag = summary_dict['javelin_lag_value']
            output['javelin']['tophat_lag'].append([lag_err_lo, med_lag, lag_err_hi])

            min_bound = summary_dict['javelin_lag_bounds'][0]
            max_bound = summary_dict['javelin_lag_bounds'][1]
            peak = summary_dict['javelin_peak']
            output['javelin']['bounds'].append([min_bound, peak, max_bound])

            lags, ntau, weight_dist, acf, smooth_dist, smooth_weight_dist = np.loadtxt( weight_dir + 'javelin_weights.dat', delimiter=',', unpack=True )
            output['javelin']['lags'].append(lags)
            output['javelin']['ntau'].append(ntau)
            output['javelin']['weight_dist'].append(weight_dist)
            output['javelin']['acf'].append(acf)
            output['javelin']['smoothed_dist'].append(smooth_weight_dist)

            weighted_lag_dist = np.loadtxt( weight_dir + 'javelin_weighted_lag_dist.dat', delimiter=',', unpack=True )
            output['javelin']['downsampled_lag_dist'].append(weighted_lag_dist)

            output['javelin']['frac_rejected'].append( summary_dict['javelin_frac_rejected'] )

        else:
            output['javelin']['tophat_lag'].append(None)
            output['javelin']['bounds'].append(None)
            output['javelin']['lags'].append(None)
            output['javelin']['ntau'].append(None)
            output['javelin']['weight_dist'].append(None)
            output['javelin']['acf'].append(None)
            output['javelin']['smoothed_dist'].append(None)
            output['javelin']['downsampled_lag_dist'].append(None)
            output['javelin']['frac_rejected'].append(None)



        if pweight:
            cent_err_lo = summary_dict['pyccf_lag_uncertainty'][0]
            cent_err_hi = summary_dict['pyccf_lag_uncertainty'][1]
            med_cent = summary_dict['pyccf_lag_value']
            output['pyccf']['centroid'].append([cent_err_lo, med_cent, cent_err_hi])

            min_bound = summary_dict['pyccf_lag_bounds'][0]
            max_bound = summary_dict['pyccf_lag_bounds'][1]
            peak = summary_dict['pyccf_peak']
            output['pyccf']['bounds'].append([min_bound, peak, max_bound])

            lags, ntau, weight_dist, acf, smooth_dist, smooth_weight_dist = np.loadtxt( weight_dir + 'pyccf_weights.dat', delimiter=',', unpack=True )
            output['pyccf']['lags'].append(lags)
            output['pyccf']['ntau'].append(ntau)
            output['pyccf']['weight_dist'].append(weight_dist)
            output['pyccf']['acf'].append(acf)
            output['pyccf']['smoothed_dist'].append(smooth_weight_dist)

            weighted_cccd = np.loadtxt( weight_dir + 'pyccf_weighted_cccd.dat', delimiter=',', unpack=True )
            output['pyccf']['downsampled_CCCD'].append(weighted_cccd)

            output['pyccf']['frac_rejected'].append( summary_dict['pyccf_frac_rejected'] )

        else:
            output['pyccf']['centroid'].append(None)
            output['pyccf']['bounds'].append(None)
            output['pyccf']['lags'].append(None)
            output['pyccf']['ntau'].append(None)
            output['pyccf']['weight_dist'].append(None)
            output['pyccf']['acf'].append(None)
            output['pyccf']['smoothed_dist'].append(None)
            output['pyccf']['downsampled_CCCD'].append(None)
            output['pyccf']['frac_rejected'].append(None)


        if jweight & pweight:
            output['rmax'].append( summary_dict['rmax'] )

        output['names'].append( line_names[i+1] )


    return output


#######################################################################
#                             All Modules
#######################################################################

def load(main_dir, verbose=False):

    """Load the results from a previous run of pyPetal.

    Parameters
    ----------

    main_dir : str
        The main output directory of the pyPetal run.

    verbose : bool, optional
        Whether or not to print the summary of the run (i.e. what modules were run).




    Returns
    -------

    output : dict
        A dictionary containing the results of the run.

    """

    main_dir = os.path.abspath(main_dir) + r'/'

    run_drw_rej, run_pyccf, run_pyzdcf, run_javelin, run_weighting = get_modules(main_dir)

    txt = """
Prior pyPetal run
---------------------
DRW Rejection: {}
pyCCF: {}
pyZDCF: {}
JAVELIN: {}
Weighting: {}
---------------------
""".format(run_drw_rej, run_pyccf, run_pyzdcf, run_javelin, run_weighting)

    if verbose:
        print(txt)


    drw_rej_res = {}
    pyccf_res = {}
    pyzdcf_res = {}
    plike_res = {}
    javelin_res = {}
    weighting_res = {}

    if run_drw_rej:
        drw_rej_res = load_drw_rej(main_dir)

    if run_pyccf:
        pyccf_res = load_pyccf(main_dir)

    if run_pyzdcf:
        pyzdcf_res, plike_res = load_pyzdcf(main_dir)

    if run_javelin:
        javelin_res = load_javelin(main_dir)

    if run_weighting:
        weighting_res = load_weighting(main_dir)

    res_tot = {}
    res_tot['drw_rej_res'] = drw_rej_res
    res_tot['pyccf_res'] = pyccf_res
    res_tot['pyzdcf_res'] = pyzdcf_res
    res_tot['plike_res'] = plike_res
    res_tot['javelin_res'] = javelin_res
    res_tot['weighting_res'] = weighting_res

    return res_tot
