import glob
import os
import shutil
import unittest

import javelin.lcmodel
import javelin.zylc
import numpy as np
import pandas as pd
from javelin.zylc import get_data

import pypetal.pipeline as pl
from pypetal.load import read_weighting_summary
from pypetal.weighting import prob_tau


def format_float(x):
    return float("%11.3e" % x)

def format_int(x):
    return int("%5.1i" % x)

def str2bool(string_, default='raise'):

    true = ['true', 't', '1', 'y', 'yes', 'enabled', 'enable', 'on']
    false = ['false', 'f', '0', 'n', 'no', 'disabled', 'disable', 'off']
    if string_.lower() in true:
        return True
    elif string_.lower() in false or (not default):
        return False
    else:
        raise ValueError('The value \'{}\' cannot be mapped to boolean.'
                         .format(string_))






class TestAll(unittest.TestCase):

    def setUp(self):

        main_dir = 'examples/dat/javelin_'
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)


        drw_rej_params = {
            'nchain': 100,
            'nburn': 50,
            'nwalker': 20,
            'reject_data': [True, False, True],
            'nsig': 1.0,
            'use_for_javelin': True
        }

        detrend_params = {
            'K': 2
        }

        pyccf_params = {
            'nsim': 100,
            'interp': 1.75,
            'sigmode': .25,
            'thres': 0.79
        }

        pyzdcf_params = {
            'nsim': 100,
            'minpts': 12,
            'prefix': 'myname',
            'run_plike': True,
            'plike_dir': 'plike_v4/'
        }

        javelin_params = {
            'nchain': 40,
            'nburn': 20,
            'nwalker': 20,
            'rm_type': "phot"
        }

        weighting_params = {
            'k': 3,
            'gap_size': 20
        }

        lag_bounds = [ [-1000, 1000], 'baseline' ]

        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                            run_drw_rej=True, drw_rej_params=drw_rej_params,
                            run_detrend=True, detrend_params=detrend_params,
                            run_pyccf=True, pyccf_params=pyccf_params,
                            run_pyzdcf=True, pyzdcf_params=pyzdcf_params,
                            run_javelin=True, javelin_params=javelin_params,
                            run_weighting=True, weighting_params=weighting_params,
                            lag_bounds=lag_bounds,
                            file_fmt='ascii',
                            time_unit='d',
                            lc_unit=['Jy', 'mJy', 'mJy'])

        self.filenames = filenames
        self.line_names = line_names
        self.res = res

        self.drw_mc = drw_rej_params['nchain'] * drw_rej_params['nwalker']
        self.drw_burn = drw_rej_params['nburn'] * drw_rej_params['nwalker']
        self.reject_data = drw_rej_params['reject_data']

        self.pyccf_mc = pyccf_params['nsim']
        self.pyzdcf_mc = pyzdcf_params['nsim']

        self.javelin_mc = javelin_params['nchain'] * javelin_params['nwalker']
        self.javelin_burn = javelin_params['nburn'] * javelin_params['nwalker']

        self.lag_bounds = lag_bounds



    def test_all(self):

        #################################################################################################
        # FILE LAYOUT
        #################################################################################################

        main_directories = glob.glob('.tmp/*/')
        subdirs = ['.tmp/continuum/', '.tmp/yelm/', '.tmp/zing/']

        #Make sure each line has a subdirectory
        self.assertIn( '.tmp/continuum/' , main_directories )
        self.assertIn( '.tmp/yelm/' , main_directories )
        self.assertIn( '.tmp/zing/' , main_directories )

        #Make sure eachline has a subdir for all modules
        for fdir, name in zip(subdirs[1:], self.line_names[1:]):
            for subdir in ['pyccf', 'pyzdcf', 'javelin', 'weights']:
                self.assertIn( '.tmp/' + name + '/' + subdir, glob.glob( fdir + '*' ) )

        self.assertIn( '.tmp/continuum/drw_rej', glob.glob( subdirs[0] + '*' ) )
        self.assertNotIn( '.tmp/yelm/drw_rej', glob.glob( subdirs[1] + '*' ) )
        self.assertIn( '.tmp/zing/drw_rej', glob.glob( subdirs[2] + '*' )  )


        #Make sure the light curves get saved
        self.assertIn( '.tmp/light_curves/', main_directories )

        lc_files = glob.glob('.tmp/light_curves/*')
        for name in self.line_names:
            self.assertIn( '.tmp/light_curves/' + name + '.dat', lc_files )


        #Make sure processed light curves get saved
        self.assertIn( '.tmp/processed_lcs/', main_directories )

        lc_files = glob.glob('.tmp/processed_lcs/*')
        for name in self.line_names:
            self.assertIn( '.tmp/processed_lcs/' + name + '_data.dat', lc_files )



        #################################################################################################
        # INDIVIDUAL FILES EXIST
        #################################################################################################

        #DRW rejection
        drw_suffix = ['_chain.dat', '_drw_fit.dat', '_mask.dat', '_drw_fit.pdf']
        for i, name in enumerate(self.line_names):
            subdir = '.tmp/' + name + r'/'
            drw_dir = subdir + 'drw_rej/'

            if not self.reject_data[i]:
                self.assertNotIn( drw_dir, glob.glob( subdir + '*/' ) )
                continue

            self.assertIn(drw_dir, glob.glob( subdir + '*/' ) )

            for suffix in drw_suffix:
                self.assertIn( drw_dir + name + suffix, glob.glob( drw_dir + '*' ) )


        #Detrend
        for name in self.line_names:
            subdir = '.tmp/' + name + r'/'
            self.assertIn( subdir + 'detrend.pdf', glob.glob( subdir + '*' ) )


        #PyCCF
        pyccf_suffix = ['_ccf.dat', '_ccf_dists.dat', '_ccf.pdf']
        for name in self.line_names[1:]:
            subdir = '.tmp/' + name + r'/'
            pyccf_dir = subdir + 'pyccf/'
            self.assertIn(pyccf_dir, glob.glob( subdir + '*/' ) )

            for suffix in pyccf_suffix:
                self.assertIn( pyccf_dir + name + suffix, glob.glob( pyccf_dir + '*' ) )


        #PyZDCF
        pyzdcf_suffix = ['_myname.dcf', '_zdcf.pdf']
        for name in self.line_names[1:]:
            subdir = '.tmp/' + name + r'/'
            pyzdcf_dir = subdir + 'pyzdcf/'
            self.assertIn(pyzdcf_dir, glob.glob( subdir + '*/' ) )

            for suffix in pyzdcf_suffix:
                self.assertIn( pyzdcf_dir + name + suffix, glob.glob( pyzdcf_dir + '*' ) )

        #PLIKE
        plike_suffix = ['_plike.out']
        for name in self.line_names[1:]:
            subdir = '.tmp/' + name + r'/'
            pyzdcf_dir = subdir + 'pyzdcf/'
            self.assertIn(plike_dir, glob.glob( subdir + '*/' ) )

            for suffix in plike_suffix:
                self.assertIn( plike_dir + name + suffix, glob.glob( pyzdcf_dir + '*' ) )


        #Javelin
        types = ['cont', 'rmap']
        javelin_suffix1 = ['burn_', 'chain_', 'logp_']
        javelin_suffix2 = ['_lc_fits.dat']
        javelin_files = ['cont_lcfile.dat', 'tot_lcfile.dat', 'javelin_bestfit.pdf', 'javelin_histogram.pdf']

        for name in self.line_names[1:]:
            subdir = '.tmp/' + name + r'/'
            javelin_dir = subdir + 'javelin/'
            self.assertIn(javelin_dir, glob.glob( subdir + '*/' ) )

            for suffix in javelin_suffix1:
                for t in ['cont']:
                    self.assertNotIn( javelin_dir + suffix + t + '.txt', glob.glob( javelin_dir + '*' ) )
                for t in ['rmap']:
                    self.assertIn( javelin_dir + suffix + t + '.txt', glob.glob( javelin_dir + '*' ) )

            for suffix in javelin_suffix2:
                self.assertIn( javelin_dir + name + suffix, glob.glob( javelin_dir + '*' ) )

            for f in javelin_files:
                self.assertIn( javelin_dir + f, glob.glob( javelin_dir + '*' ) )

            self.assertNotIn( javelin_dir + 'javelin_corner.pdf', glob.glob( javelin_dir + '*' ) )


        #Weighting
        files1 = ['javelin_weighted_lag_dist.dat', 'javelin_weights.dat', 'pyccf_weighted_cccd.dat', 'pyccf_weights.dat', 'weight_summary.txt']
        files2 = ['javelin_weights_res.pdf', 'pyccf_weights_res.pdf']
        for name in self.line_names[1:]:
            subdir = '.tmp/' + name + r'/'
            weight_dir = subdir + 'weights/'
            self.assertIn( weight_dir, glob.glob( subdir + '*/' ) )

            self.assertIn( weight_dir + name + '_weights.pdf', glob.glob(weight_dir + '*') )

            for f in files1:
                self.assertIn( weight_dir + f, glob.glob( weight_dir + '*' ) )

            for f in files2:
                self.assertIn( '.tmp/' + f, glob.glob( '.tmp/*' ) )



        #Light Curves
        self.assertIn( '.tmp/light_curves', glob.glob('.tmp/*') )
        for name in self.line_names:
            self.assertIn( '.tmp/light_curves/' + name + '.dat', glob.glob('.tmp/light_curves/*') )

        #Processed Light Curves
        self.assertIn( '.tmp/processed_lcs', glob.glob('.tmp/*') )
        for name in self.line_names:
            self.assertIn( '.tmp/processed_lcs/' + name + '_data.dat', glob.glob('.tmp/processed_lcs/*') )



        #################################################################################################
        # RES KEYS
        #################################################################################################

        general_keys = ['drw_rej_res', 'pyccf_res', 'pyzdcf_res', 'plike_res', 'javelin_res', 'weighting_res']

        drw_rej_keys = ['masks', 'reject_data', 'taus', 'sigmas', 'jitters']

        pyccf_keys = ['CCF', 'CCF_lags', 'centroid', 'centroid_err_lo', 'centroid_err_hi',
                        'peak', 'peak_err_lo', 'peak_err_hi', 'CCCD_lags', 'CCPD_lags', 'name']

        pyzdcf_cols = ['tau', '-sig(tau)', '+sig(tau)', 'dcf', '-err(dcf)', '+err(dcf)', '#bin' ]

        plike_keys = ['output', 'ML_lag', 'ML_lag_err_lo', 'ML_lag_err_hi']
        plike_cols = ['num', 'lag', 'r', '-dr', '+dr', 'likelihood']

        javelin_keys = ['cont_hpd', 'tau', 'sigma', 'tophat_params', 'hpd',
                        'cont_model', 'rmap_model', 'cont_dat', 'tot_dat', 'bestfit_model']

        weight_keys = ['pyccf', 'javelin']
        weight_keys1 = ['centroid', 'bounds', 'acf', 'lags', 'weight_dist',
                        'smoothed_dist', 'ntau', 'downsampled_CCCD', 'frac_rejected',
                        'CCCD']
        weight_keys2 = ['tophat_lag', 'bounds', 'acf', 'lags', 'weight_dist',
                        'smoothed_dist', 'ntau', 'downsampled_lag_dist', 'frac_rejected',
                        'lag_dist']


        keys_tot = [drw_rej_keys, pyccf_keys, None, plike_keys, javelin_keys, weight_keys]


        #DRW Rej, pyCCF, PLIKE, JAVELIN
        for i, key in enumerate(general_keys):
            self.assertIn(key, self.res.keys())

            if keys_tot[i] is None:
                continue

            for k in keys_tot[i]:

                if i in [1, 3, 4]:
                    for j in range(len(self.line_names[1:])):
                        self.assertIn(k, list(self.res[key][j].keys()) )
                else:
                    self.assertIn(k, list(self.res[key].keys()) )


        #pyZDCF
        for i in range(len(self.line_names[1:])):
            self.assertListEqual( pyzdcf_cols, list(self.res['pyzdcf_res'][i].columns) )

        #PLIKE
        for i in range(len(self.line_names[1:])):
            self.assertListEqual( plike_cols, list(self.res['plike_res'][i]['output'].colnames) )

        #Weighting
        for mod, keys in zip(weight_keys, [weight_keys1, weight_keys2]):
            for k in keys:
                self.assertIn(k, self.res['weighting_res'][mod].keys())




        #################################################################################################
        # RES VALUES
        #################################################################################################
        #Test lengths of arrays and their types

        #DRW rejection
        for key in ['taus', 'sigmas', 'jitters']:
            self.assertEqual( len(self.res['drw_rej_res'][key]), len(self.line_names) )

        for i in range(len(self.line_names)):
            x, _, _ = np.loadtxt( self.filenames[i], unpack=True, usecols=[0,1,2] )
            self.assertEqual( len(self.res['drw_rej_res']['masks'][i]), len(x) )
            self.assertIs( self.res['drw_rej_res']['masks'][i].dtype.type, np.bool_ )

            for key in ['taus', 'sigmas', 'jitters']:
                if self.reject_data[i]:
                    self.assertEqual( len(self.res['drw_rej_res'][key][i]), self.drw_mc )
                    self.assertIs( np.array(self.res['drw_rej_res'][key][i]).dtype.type, np.float64 )
                else:
                    self.assertIs( self.res['drw_rej_res'][key][i], None )


        #PyCCF
        self.assertEqual( len(self.res['pyccf_res']), len(self.line_names)-1 )
        for i in range(len(self.line_names)-1):

            for key in [ 'CCF', 'CCF_lags', 'centroid', 'centroid_err_lo', 'centroid_err_hi',
                        'peak', 'peak_err_lo', 'peak_err_hi', 'CCCD_lags', 'CCPD_lags' ]:
                self.assertIs( self.res['pyccf_res'][i][key].dtype.type, np.float64 )

            self.assertIs( type(self.res['pyccf_res'][i]['name']), str )

            self.assertEqual( len(self.res['pyccf_res'][i]['CCCD_lags']), self.pyccf_mc )
            self.assertEqual( len(self.res['pyccf_res'][i]['CCPD_lags']), self.pyccf_mc )


        #pyZDCF
        for i in range(len(self.line_names[1:])):
            self.assertIs( type(self.res['pyzdcf_res'][i]), pd.DataFrame )

        #PLIKE
        for i in range(len(self.line_names[1:])):
            self.assertIs( type(self.res['plike_res'][i]['output']), Table )

            for key in ['ML_lag', 'ML_lag_err_lo', 'ML_lag_err_hi']:
                self.assertIn( type(self.res['plike_res'][i][key]), [float, np.float64] )

            self.assertGreaterEqual( np.min(self.res['plike_res'][i]['output']['lag']), self.lag_bounds[i][0] )
            self.assertLessEqual( np.max(self.res['plike_res'][i]['output']['lag']), self.lag_bounds[i][1] )


        #Javelin
        self.assertEqual( len(self.res['javelin_res']), len(self.line_names)-1 )
        for i in range(len(self.line_names)-1):

            for key in [ 'tau', 'sigma' ]:
                self.assertIs( type(self.res['javelin_res'][i][key][0]), np.float64 )

            for key in ['tophat_params', 'hpd']:
                self.assertIs( type(self.res['javelin_res'][i][key][0][0]), np.float64 )


            self.assertEqual( len(self.res['javelin_res'][i]['tau']), self.javelin_mc )
            self.assertEqual( len(self.res['javelin_res'][i]['sigma']), self.javelin_mc )

            self.assertEqual( len(self.res['javelin_res'][i]['tophat_params']), 4 )
            for j in range(4):
                self.assertEqual( len(self.res['javelin_res'][i]['tophat_params'][j]), self.javelin_mc)

            self.assertIs( self.res['javelin_res'][i]['cont_hpd'], None )
            self.assertEqual( len(self.res['javelin_res'][i]['hpd']), 3 )
            for j in range(3):
                self.assertEqual( len(self.res['javelin_res'][i]['hpd'][j]), 6)


            self.assertIs( self.res['javelin_res'][i]['cont_model'], None )
            self.assertIs( type(self.res['javelin_res'][i]['rmap_model']), javelin.lcmodel.Pmap_Model )

            for key in ['cont_dat', 'tot_dat', 'bestfit_model']:
                self.assertIs( type(self.res['javelin_res'][i][key]), javelin.zylc.LightCurve )

            self.assertGreaterEqual( np.min(self.res['javelin_res'][i]['tophat_params'][0]), self.lag_bounds[i][0] )
            self.assertLessEqual( np.max(self.res['javelin_res'][i]['tophat_params'][0]), self.lag_bounds[i][1] )



        #Weighting
        modules = ['pyccf', 'javelin']
        for mod in modules:

            for i in range(len(self.line_names[1:])):
                self.assertEqual( len(self.res['weighting_res'][mod]['bounds'][i]), 3 )

                for key in ['bounds', 'acf', 'lags', 'weight_dist', 'smoothed_dist']:
                    self.assertIs( type(self.res['weighting_res'][mod][key][i][0]), np.float64 )

                self.assertIs( type(self.res['weighting_res'][mod]['ntau'][i][0]), np.float64 )
                self.assertIs( type(self.res['weighting_res'][mod]['frac_rejected'][i]), float )

                if mod == 'pyccf':
                    self.assertIs( type(self.res['weighting_res']['pyccf']['downsampled_CCCD'][i][0]), np.float64 )

                    self.assertEqual( len(self.res['weighting_res'][mod]['centroid'][i]), 3 )
                    self.assertIs( type(self.res['weighting_res'][mod]['centroid'][i][0]), np.float64 )

                    self.assertAlmostEqual( len(self.res['pyccf_res'][i]['CCCD_lags'])*(1 - self.res['weighting_res'][mod]['frac_rejected'][i]), len(self.res['weighting_res'][mod]['downsampled_CCCD'][i]), places=5 )

                    #Make sure the downsampled data is within the peaks bounds
                    for j in range(len(self.res['weighting_res'][mod]['downsampled_CCCD'][i])):
                        self.assertGreaterEqual( self.res['weighting_res'][mod]['downsampled_CCCD'][i][j], self.res['weighting_res'][mod]['bounds'][i][0] )
                        self.assertLessEqual( self.res['weighting_res'][mod]['downsampled_CCCD'][i][j], self.res['weighting_res'][mod]['bounds'][i][2] )


                    #Make sure the peak is the same as the peak of the downsampled data
                    peak = np.median(self.res['weighting_res'][mod]['downsampled_CCCD'][i])
                    err_lo = peak - np.percentile(self.res['weighting_res'][mod]['downsampled_CCCD'][i], 16)
                    err_hi = np.percentile(self.res['weighting_res'][mod]['downsampled_CCCD'][i], 84) - peak

                    self.assertEqual( peak, self.res['weighting_res'][mod]['centroid'][i][1] )
                    self.assertEqual( err_lo, self.res['weighting_res'][mod]['centroid'][i][0] )
                    self.assertEqual( err_hi, self.res['weighting_res'][mod]['centroid'][i][2] )


                    dist_ds = self.res['weighting_res'][mod]['downsampled_CCCD'][i]
                    dist = self.res['pyccf_res'][i]['CCCD_lags']

                else:
                    self.assertIs( type(self.res['weighting_res']['javelin']['downsampled_lag_dist'][i][0]), np.float64 )

                    self.assertEqual( len(self.res['weighting_res'][mod]['tophat_lag'][i]), 3 )
                    self.assertIs( type(self.res['weighting_res'][mod]['tophat_lag'][i][0]), np.float64 )

                    self.assertAlmostEqual( len(self.res['javelin_res'][i]['tophat_params'][0])*(1 - self.res['weighting_res'][mod]['frac_rejected'][i]), len(self.res['weighting_res'][mod]['downsampled_lag_dist'][i]), places=5 )

                    for j in range(len(self.res['weighting_res'][mod]['downsampled_lag_dist'][i])):
                        self.assertGreaterEqual( self.res['weighting_res'][mod]['downsampled_lag_dist'][i][j], self.res['weighting_res'][mod]['bounds'][i][0] )
                        self.assertLessEqual( self.res['weighting_res'][mod]['downsampled_lag_dist'][i][j], self.res['weighting_res'][mod]['bounds'][i][2] )


                    peak = np.median(self.res['weighting_res'][mod]['downsampled_lag_dist'][i])
                    err_lo = peak - np.percentile(self.res['weighting_res'][mod]['downsampled_lag_dist'][i], 16)
                    err_hi = np.percentile(self.res['weighting_res'][mod]['downsampled_lag_dist'][i], 84) - peak

                    self.assertEqual( peak, self.res['weighting_res'][mod]['tophat_lag'][i][1] )
                    self.assertEqual( err_lo, self.res['weighting_res'][mod]['tophat_lag'][i][0] )
                    self.assertEqual( err_hi, self.res['weighting_res'][mod]['tophat_lag'][i][2] )


                    dist_ds = self.res['weighting_res'][mod]['downsampled_lag_dist'][i]
                    dist = self.res['javelin_res'][i]['tophat_params'][0]




                self.assertListEqual( list(self.res['weighting_res']['pyccf']['CCCD'][i]), list(self.res['pyccf_res'][i]['CCCD_lags']) )
                self.assertListEqual( list(self.res['weighting_res']['javelin']['lag_dist'][i]), list(self.res['javelin_res'][i]['tophat_params'][0]) )

                #Make sure the downsampled data points are in the original sample
                for j in range(len(dist_ds)):
                    self.assertIn( dist_ds[j], dist )




        #################################################################################################
        # INDIVIDUAL FILE CONTENT
        #################################################################################################

        #DRW rejection
        for i in range(len(self.line_names)):

            if not self.reject_data[i]:
                x, _, _ = np.loadtxt( '.tmp/light_curves/' + self.line_names[i] + '.dat', unpack=True, delimiter=',', usecols=[0,1,2] )
                self.assertListEqual( list(np.zeros(len(x), dtype=bool)), list( self.res['drw_rej_res']['masks'][i] ) )
                continue

            #Chains
            file_sig, file_tau, file_jit = np.loadtxt( '.tmp/' + self.line_names[i] + '/drw_rej/' + self.line_names[i] + '_chain.dat', unpack=True, delimiter=',' )
            self.assertListEqual( list(file_sig), list(self.res['drw_rej_res']['sigmas'][i]) )
            self.assertListEqual( list(file_tau), list(self.res['drw_rej_res']['taus'][i]) )
            self.assertListEqual( list(file_jit), list(self.res['drw_rej_res']['jitters'][i]) )

            #Masks
            file_mask = np.loadtxt( '.tmp/' + self.line_names[i] + '/drw_rej/' + self.line_names[i] + '_mask.dat', dtype=str)
            file_mask = list(map(str2bool, file_mask))
            self.assertListEqual( list(file_mask), list( self.res['drw_rej_res']['masks'][i] ) )


        #PyCCF
        for i in range(len(self.line_names[1:])):

            #Dists
            file_cccd, file_ccpd = np.loadtxt( '.tmp/' + self.line_names[i+1] + '/pyccf/' + self.line_names[i+1] + '_ccf_dists.dat', unpack=True, delimiter=',' )
            self.assertListEqual( list(file_cccd), list(self.res['pyccf_res'][i]['CCCD_lags']) )
            self.assertListEqual( list(file_ccpd), list(self.res['pyccf_res'][i]['CCPD_lags']) )

            #CCF
            file_lags, file_ccf = np.loadtxt( '.tmp/' + self.line_names[i+1] + '/pyccf/' + self.line_names[i+1] + '_ccf.dat', unpack=True, delimiter=',' )
            self.assertListEqual( list(file_lags), list(self.res['pyccf_res'][i]['CCF_lags']) )
            self.assertListEqual( list(file_ccf), list(self.res['pyccf_res'][i]['CCF']) )


        #PyZDCF
        df_cols = ['tau', '-sig(tau)', '+sig(tau)', 'dcf', '-err(dcf)', '+err(dcf)', '#bin' ]
        for i in range(len(self.line_names[1:])):
            file_df = np.loadtxt('.tmp/' + self.line_names[i+1] + '/pyzdcf/' + self.line_names[i+1] + '_myname.dcf')
            file_df = pd.DataFrame(file_df, columns=df_cols)

            #Read the res DF with the proper formatting
            res_df = self.res['pyzdcf_res'][i].copy(deep=True)
            for col in df_cols:
                if col != '#bin':
                    res_df[col] = list(map( format_float, res_df[col] ))
                else:
                    res_df[col] = list(map( format_int, res_df[col] ))

            for col in df_cols:
                self.assertListEqual( list(file_df[col]), list(res_df[col]) )



        #PLIKE
        df_cols = ['num', 'lag', 'r', '-dr', '+dr', 'likelihood']
        for i in range(len(self.line_names[1:])):
            file_df = Table.read('.tmp/' + self.line_names[i+1] + '/pyzdcf/' + self.line_names[i+1] + '_plike.out',
                                 format='ascii',
                                 names=['num', 'lag', 'r', '-dr', '+dr', 'likelihood'])

            res_df = self.res['plike_res'][i]['output']

            for col in df_cols:
                self.assertListEqual( list(res_df[col]), list(file_df[col]) )


            with open('.tmp/' + self.line_names[i+1] + '/pyzdcf/' + self.line_names[i+1] + '_plike.out', 'r') as f:
                output_str = list(f)[-3:]

                ml_lag = float( output_str[1].split()[7] )
                ml_lag_err_hi = np.abs( float( output_str[1].split()[8] )  )
                ml_lag_err_lo = np.abs( float( output_str[1].split()[9] )  )


            self.assertEqual( ml_lag, self.res['plike_res'][i]['ML_lag'] )
            self.assertEqual( ml_lag_err_hi, self.res['plike_res'][i]['ML_lag_err_hi'] )
            self.assertEqual( ml_lag_err_lo, self.res['plike_res'][i]['ML_lag_err_lo'] )





        #Javelin
        for i in range(len(self.line_names[1:])):
            #Chains
            file_sig, file_tau, t, w, s, a = np.loadtxt( '.tmp/' + self.line_names[i+1] + '/javelin/chain_rmap.txt', unpack=True )
            file_tophat = [t.tolist(), w.tolist(), s.tolist(), a.tolist()]

            self.assertListEqual( list(np.exp(file_sig)), list(self.res['javelin_res'][i]['sigma']) )
            self.assertListEqual( list(np.exp(file_tau)), list(self.res['javelin_res'][i]['tau']) )
            self.assertListEqual( file_tophat, self.res['javelin_res'][i]['tophat_params'].tolist() )


            #Use for javelin
            drw_tau = np.median(self.res['drw_rej_res']['taus'][0])
            drw_sig = np.median(self.res['drw_rej_res']['sigmas'][0])

            for j in range(self.javelin_mc):
                self.assertAlmostEqual( np.exp(file_sig[j]), drw_sig, places=5 )
                self.assertAlmostEqual( np.exp(file_tau[j]), drw_tau, places=5 )


            #Light curve files
            file_cont = get_data('.tmp/' + self.line_names[i+1] + '/javelin/cont_lcfile.dat')
            file_tot = get_data('.tmp/' + self.line_names[i+1] + '/javelin/tot_lcfile.dat')

                #Continuum
            file_ydat = file_cont.mlist[0] + file_cont.blist[0]
            res_ydat = self.res['javelin_res'][i]['cont_dat'].mlist[0] + self.res['javelin_res'][i]['cont_dat'].blist[0]
            self.assertEqual(file_cont.nlc, self.res['javelin_res'][i]['cont_dat'].nlc)
            self.assertListEqual(file_cont.jlist[0].tolist(), self.res['javelin_res'][i]['cont_dat'].jlist[0].tolist() )
            self.assertListEqual(file_ydat.tolist(), res_ydat.tolist() )
            self.assertListEqual(file_cont.elist[0].tolist(), self.res['javelin_res'][i]['cont_dat'].elist[0].tolist() )
            self.assertListEqual(file_cont.ilist[0].tolist(), self.res['javelin_res'][i]['cont_dat'].ilist[0].tolist() )

            for j in range(len(file_cont.jlist[0])):
                file_ydat = file_cont.mlist[0][j] + file_cont.blist[0]
                res_ydat = self.res['javelin_res'][i]['cont_dat'].mlist[0][j] + self.res['javelin_res'][i]['cont_dat'].blist[0]
                self.assertAlmostEqual(file_cont.jlist[0][j], self.res['javelin_res'][i]['tot_dat'].jlist[0][j], places=4 )
                self.assertAlmostEqual(file_ydat, res_ydat, places=4 )
                self.assertAlmostEqual(file_cont.elist[0][j], self.res['javelin_res'][i]['tot_dat'].elist[0][j], places=4 )
                self.assertAlmostEqual(file_cont.ilist[0][j], self.res['javelin_res'][i]['tot_dat'].ilist[0][j], places=4 )

                #Tot
            self.assertEqual(file_tot.nlc, self.res['javelin_res'][i]['tot_dat'].nlc)
            for j in range(len(self.filenames[1:])):
                for k in range(len(file_tot.jlist[j])):
                    file_ydat = file_tot.mlist[j][k] + file_tot.blist[j]
                    res_ydat = self.res['javelin_res'][i]['tot_dat'].mlist[j][k] + self.res['javelin_res'][i]['tot_dat'].blist[j]
                    self.assertAlmostEqual(file_tot.jlist[j][k], self.res['javelin_res'][i]['tot_dat'].jlist[j][k], places=3 )
                    self.assertAlmostEqual(file_ydat, res_ydat, places=3 )
                    self.assertAlmostEqual(file_tot.elist[j][k], self.res['javelin_res'][i]['tot_dat'].elist[j][k], places=3 )
                    self.assertAlmostEqual(file_tot.ilist[j][k], self.res['javelin_res'][i]['tot_dat'].ilist[j][k], places=3 )



        #Weighting

            #Dists
        for i, name in enumerate(self.line_names[1:]):
            lag_dist_ds = np.loadtxt('.tmp/' + name + '/weights/javelin_weighted_lag_dist.dat')
            cccd_ds = np.loadtxt('.tmp/' + name + '/weights/pyccf_weighted_cccd.dat')

            self.assertEqual( list(self.res['weighting_res']['javelin']['downsampled_lag_dist'][i]), list(lag_dist_ds) )
            self.assertEqual( list(self.res['weighting_res']['pyccf']['downsampled_CCCD'][i]), list(cccd_ds) )

            #Make sure downsampled dists are in the original ones
            for j in range(len(lag_dist_ds)):
                self.assertIn( lag_dist_ds[j], self.res['javelin_res'][i]['tophat_params'][0] )

            for j in range(len(cccd_ds)):
                self.assertIn( cccd_ds[j], self.res['pyccf_res'][i]['CCCD_lags'] )


            lags, ntau, weight_dist, acf, smooth_dist, smooth_weight_dist = np.loadtxt( '.tmp/' + name + '/weights/javelin_weights.dat', unpack=True, delimiter=',' )
            self.assertListEqual( list(self.res['weighting_res']['javelin']['lags'][i]), list(lags) )
            self.assertListEqual( list(self.res['weighting_res']['javelin']['ntau'][i]), list(ntau) )
            self.assertListEqual( list(self.res['weighting_res']['javelin']['weight_dist'][i]), list(weight_dist) )
            self.assertListEqual( list(self.res['weighting_res']['javelin']['acf'][i]), list(acf) )
            self.assertListEqual( list(self.res['weighting_res']['javelin']['smoothed_dist'][i]), list(smooth_weight_dist) )

            lags, ntau, weight_dist, acf, smooth_dist, smooth_weight_dist = np.loadtxt( '.tmp/' + name + '/weights/pyccf_weights.dat', unpack=True, delimiter=',' )
            self.assertListEqual( list(self.res['weighting_res']['pyccf']['lags'][i]), list(lags) )
            self.assertListEqual( list(self.res['weighting_res']['pyccf']['ntau'][i]), list(ntau) )
            self.assertListEqual( list(self.res['weighting_res']['pyccf']['weight_dist'][i]), list(weight_dist) )
            self.assertListEqual( list(self.res['weighting_res']['pyccf']['acf'][i]), list(acf) )
            self.assertListEqual( list(self.res['weighting_res']['pyccf']['smoothed_dist'][i]), list(smooth_weight_dist) )



            #Summary files
        xcont, _, _ = np.loadtxt(self.filenames[0], unpack=True)
        for i, name in enumerate(self.line_names[1:]):
            xline, _, _ = np.loadtxt(self.filenames[i+1], unpack=True)
            weight_dict = read_weighting_summary( '.tmp/' + name + '/weights/weight_summary.txt' )

            #k
            self.assertEqual( weight_dict['k'], 3 )


            #N0
            pyccf_lags = self.res['weighting_res']['pyccf']['lags'][i]
            jav_lags = self.res['weighting_res']['javelin']['lags'][i]

            _, _, pyccf_n0, _ = prob_tau(xcont, xline, lagvals=pyccf_lags, gap_size=20, k=3)
            _, _, jav_n0, _ = prob_tau(xcont, xline, lagvals=jav_lags, gap_size=20, k=3)
            self.assertLessEqual( np.abs(weight_dict['pyccf_n0'] - pyccf_n0), 1 )
            self.assertLessEqual( np.abs(weight_dict['javelin_n0'] - jav_n0), 1 )



            #Bounds
            pyccf_bounds = [ self.res['weighting_res']['pyccf']['bounds'][i][0], self.res['weighting_res']['pyccf']['bounds'][i][2] ]
            pyccf_peak = self.res['weighting_res']['pyccf']['bounds'][i][1]

            jav_bounds = [ self.res['weighting_res']['javelin']['bounds'][i][0], self.res['weighting_res']['javelin']['bounds'][i][2] ]
            jav_peak = self.res['weighting_res']['javelin']['bounds'][i][1]

            self.assertListEqual( list(weight_dict['pyccf_lag_bounds']), list(pyccf_bounds) )
            self.assertEqual( weight_dict['pyccf_peak'], pyccf_peak )
            self.assertListEqual( list(weight_dict['javelin_lag_bounds']), list(jav_bounds) )
            self.assertEqual( weight_dict['javelin_peak'], jav_peak )



            #Lag
            pyccf_lag_err = [ self.res['weighting_res']['pyccf']['centroid'][i][0], self.res['weighting_res']['pyccf']['centroid'][i][2] ]
            pyccf_lag = self.res['weighting_res']['pyccf']['centroid'][i][1]

            jav_lag_err = [ self.res['weighting_res']['javelin']['tophat_lag'][i][0], self.res['weighting_res']['javelin']['tophat_lag'][i][2] ]
            jav_lag = self.res['weighting_res']['javelin']['tophat_lag'][i][1]

            self.assertListEqual( list(weight_dict['pyccf_lag_uncertainty']), list(pyccf_lag_err) )
            self.assertEqual( weight_dict['pyccf_lag_value'], pyccf_lag )
            self.assertListEqual( list(weight_dict['javelin_lag_uncertainty']), list(jav_lag_err) )
            self.assertEqual( weight_dict['javelin_lag_value'], jav_lag )


            #Rmax
            self.assertEqual( weight_dict['rmax'], self.res['weighting_res']['rmax'][i] )


            #Frac rejected
            cccd = self.res['weighting_res']['pyccf']['CCCD'][i]
            lag_dist = self.res['weighting_res']['javelin']['lag_dist'][i]

            cccd_ds = self.res['weighting_res']['pyccf']['downsampled_CCCD'][i]
            lag_dist_ds = self.res['weighting_res']['javelin']['downsampled_lag_dist'][i]

            self.assertEqual( weight_dict['pyccf_frac_rejected'], 1- len(cccd_ds)/len(cccd) )
            self.assertEqual( weight_dict['javelin_frac_rejected'], 1- len(lag_dist_ds)/len(lag_dist) )





    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
