import glob
import os
import shutil
import unittest

import numpy as np

import pypetal.load as load
import pypetal.pipeline as pl


class TestLoad(unittest.TestCase):

    def setUp(self):

        main_dir = 'examples/dat/javelin_'
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']

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
            'rm_type': "spec",
            'together': True
        }

        weighting_params = {
            'k': 3,
            'gap_size': 20
        }

        lag_bounds = 'baseline'

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

        self.loaded_res = load('.tmp/')
        self.reject_data = drw_rej_params['reject_data']





    def test_match(self):

        ###############################################################################################
        # DRW REJECTION
        ###############################################################################################

        res_keys = list(self.res['drw_rej_res'].keys())
        load_keys = list(self.loaded_res['drw_rej_res'].keys())

        for key in res_keys:
            self.assertIn(key, load_keys)
        for key in ['fit_x', 'fit_y', 'fit_err', 'names']:
            self.assertIn(key, load_keys)


        #Names
        self.assertListEqual( self.loaded_res['drw_rej_res']['names'], self.line_names )

        #Reject_data
        self.assertListEqual( self.loaded_res['drw_rej_res']['reject_data'].tolist(), self.res['drw_rej_res']['reject_data'] )

        #Chains and masks
        for i in range(len(self.filenames)):

            if not self.reject_data[i]:
                self.assertIsNone( self.loaded_res['drw_rej_res']['taus'][i] )
                self.assertIsNone( self.loaded_res['drw_rej_res']['sigmas'][i] )
                self.assertIsNone( self.loaded_res['drw_rej_res']['jitters'][i] )

                self.assertIsNone( self.res['drw_rej_res']['taus'][i] )
                self.assertIsNone( self.res['drw_rej_res']['sigmas'][i] )
                self.assertIsNone( self.res['drw_rej_res']['jitters'][i] )
                continue

            self.assertListEqual( list(self.res['drw_rej_res']['taus'][i]), list(self.loaded_res['drw_rej_res']['taus'][i]) )
            self.assertListEqual( list(self.res['drw_rej_res']['sigmas'][i]), list(self.loaded_res['drw_rej_res']['sigmas'][i]) )
            self.assertListEqual( list(self.res['drw_rej_res']['jitters'][i]), list(self.loaded_res['drw_rej_res']['jitters'][i]) )
            self.assertListEqual( list(self.res['drw_rej_res']['masks'][i]), list(self.loaded_res['drw_rej_res']['masks'][i]) )


        #Fits
        for i, name in enumerate(self.line_names):

            if not self.reject_data[i]:
                self.assertIsNone( self.loaded_res['drw_rej_res']['fit_x'][i] )
                self.assertIsNone( self.loaded_res['drw_rej_res']['fit_y'][i] )
                self.assertIsNone( self.loaded_res['drw_rej_res']['fit_err'][i] )
                continue

            xfit, yfit, yerrfit = np.loadtxt( '.tmp/' + name + '/drw_rej/' + name + '_drw_fit.dat', unpack=True, delimiter=',' )
            self.assertListEqual( list(xfit), list(self.loaded_res['drw_rej_res']['fit_x'][i]) )
            self.assertListEqual( list(yfit), list(self.loaded_res['drw_rej_res']['fit_y'][i]) )
            self.assertListEqual( list(yerrfit), list(self.loaded_res['drw_rej_res']['fit_err'][i]) )


        ###############################################################################################
        # PYCCF
        ###############################################################################################

        for i, name in enumerate(self.line_names[1:]):
            res_keys = list(self.res['pyccf_res'][i].keys())
            load_keys = list(self.loaded_res['pyccf_res'][i].keys())


            for key in load_keys:
                self.assertIn(key, res_keys)

                if key == 'name':
                    self.assertEqual(self.res['pyccf_res'][i][key], self.loaded_res['pyccf_res'][i][key])
                else:
                    self.assertListEqual(self.res['pyccf_res'][i][key].tolist(), self.loaded_res['pyccf_res'][i][key].tolist())


        ###############################################################################################
        # PYZDCF
        ###############################################################################################

        for i in range(len(self.filenames[1:])):
            load_keys = list(self.loaded_res['pyzdcf_res'][i].keys())
            self.assertListEqual(load_keys, ['output', 'name'])

            res_cols = list(self.res['pyzdcf_res'][i].columns)
            load_cols = list(self.loaded_res['pyzdcf_res'][i]['output'].columns)
            self.assertListEqual(res_cols, load_cols)

            for col in res_cols:
                if col == 'tau':
                    for j in range(len(self.res['pyzdcf_res'][i][col])):
                        res_int = int(self.res['pyzdcf_res'][i][col].iloc[j])
                        load_int = int(self.loaded_res['pyzdcf_res'][i]['output'][col].iloc[j])
                        self.assertLessEqual( np.abs( res_int - load_int ), 1 )

                elif col in ['-sig(tau)', '+sig(tau)']:
                    for j in range(len(self.res['pyzdcf_res'][i][col])):
                        res_val = self.res['pyzdcf_res'][i][col].iloc[j]
                        load_val = self.loaded_res['pyzdcf_res'][i]['output'][col].iloc[j]
                        self.assertAlmostEqual( res_val, load_val, places=2 )

                elif col in ['dcf', '-err(dcf)', '+err(dcf)']:
                    for j in range(len(self.res['pyzdcf_res'][i][col])):
                        res_val = self.res['pyzdcf_res'][i][col].iloc[j]
                        load_val = self.loaded_res['pyzdcf_res'][i]['output'][col].iloc[j]
                        self.assertAlmostEqual( res_val, load_val, places=4 )

                else:
                    self.assertListEqual( list(self.res['pyzdcf_res'][i][col]), list(self.loaded_res['pyzdcf_res'][i]['output'][col]) )


        ###############################################################################################
        # PLIKE
        ###############################################################################################

        for i in range(len(self.filenames[1:])):
            plike_cols = ['num', 'lag', 'r', '-dr', '+dr', 'likelihood']

            plike_keys = ['output', 'ML_lag', 'ML_lag_err_lo', 'ML_lag_err_hi', 'name']
            res_keys = list(self.res['plike_res'][i].keys())
            load_keys = list(self.loaded_res['plike_res'][i].keys())

            self.assertIn( 'name', load_keys )
            for key in plike_keys[:-1]:
                self.assertIn(key, load_keys)
                self.assertIn(key, res_keys)

            self.assertEqual( self.line_names[i+1], self.loaded_res['plike_res'][i]['name'] )

            for key in plike_keys[1:-1]:
                self.assertEqual( self.res['plike_res'][i][key], self.loaded_res['plike_res'][i][key] )

            self.assertListEqual( list(self.res['plike_res'][i]['output'].colnames), plike_cols )
            self.assertListEqual( list(self.loaded_res['plike_res'][i]['output'].colnames), plike_cols )
            for col in plike_cols:
                self.assertListEqual( list(self.res['plike_res'][i]['output'][col]), list(self.loaded_res['plike_res'][i]['output'][col]) )


        ###############################################################################################
        # JAVELIN
        ###############################################################################################

        jav_keys = ['tau', 'sigma', 'tophat_params',
                    'cont_fit_x', 'cont_fit_y', 'cont_fit_yerr']
        res_keys = list(self.res['javelin_res'].keys())
        load_keys = list(self.loaded_res['javelin_res'].keys())

        for key in jav_keys:
            self.assertIn(key, load_keys)
            for name in self.line_names[1:]:
                self.assertIn( name + '_fit_x', load_keys )
                self.assertIn( name + '_fit_y', load_keys )
                self.assertIn( name + '_fit_yerr', load_keys )


        for key in ['tau', 'sigma', 'tophat_params']:

            if key in ['tau','sigma']:
                self.assertListEqual(self.res['javelin_res'][key].tolist(), self.loaded_res['javelin_res'][key].tolist())
            elif key == 'tophat_params':
                for j in range(6):
                    self.assertListEqual(self.res['javelin_res'][key][j].tolist(), self.loaded_res['javelin_res'][key][j].tolist())


        c_xfit, c_yfit, c_yerrfit = np.loadtxt( '.tmp/javelin/continuum_lc_fits.dat', unpack=True, delimiter=',' )
        self.assertListEqual(c_xfit.tolist(), self.loaded_res['javelin_res']['cont_fit_x'].tolist() )
        self.assertListEqual(c_yfit.tolist(), self.loaded_res['javelin_res']['cont_fit_y'].tolist() )
        self.assertListEqual(c_yerrfit.tolist(), self.loaded_res['javelin_res']['cont_fit_yerr'].tolist() )


        for i in range(1, len(self.line_names)):
            xfit, yfit, yerrfit = np.loadtxt( '.tmp/javelin/' + self.line_names[i] + '_lc_fits.dat', unpack=True, delimiter=',' )
            self.assertListEqual(xfit.tolist(), self.loaded_res['javelin_res'][ self.line_names[i] + '_fit_x'].tolist() )
            self.assertListEqual(yfit.tolist(), self.loaded_res['javelin_res'][ self.line_names[i] + '_fit_y'].tolist() )
            self.assertListEqual(yerrfit.tolist(), self.loaded_res['javelin_res'][ self.line_names[i] + '_fit_yerr'].tolist() )




        ###############################################################################################
        # WEIGHTING
        ###############################################################################################

        res_keys = list(self.res['weighting_res'].keys())
        load_keys = list(self.loaded_res['weighting_res'].keys())

        self.assertIn( 'names', load_keys )
        for key in res_keys:
            self.assertIn(key, load_keys)


        weight_keys1 = ['centroid', 'bounds', 'acf', 'lags', 'weight_dist',
                        'smoothed_dist', 'ntau', 'downsampled_CCCD', 'frac_rejected']
        weight_keys2 = ['tophat_lag', 'bounds', 'acf', 'lags', 'weight_dist',
                        'smoothed_dist', 'ntau', 'downsampled_lag_dist', 'frac_rejected']


        for key in weight_keys1:
            self.assertIn(key, self.loaded_res['weighting_res']['pyccf'].keys())
        for key in weight_keys2:
            self.assertIn(key, self.loaded_res['weighting_res']['javelin'].keys())

        self.assertNotIn( 'CCCD', self.loaded_res['weighting_res']['pyccf'].keys() )
        self.assertNotIn( 'lag_dist', self.loaded_res['weighting_res']['javelin'].keys() )

        for key in ['bounds', 'acf', 'lags', 'weight_dist', 'smoothed_dist', 'ntau']:
            for mod in ['pyccf', 'javelin']:
                for i in range(len(self.line_names[1:])):

                    self.assertListEqual( list(self.res['weighting_res'][mod][key][i]), list(self.loaded_res['weighting_res'][mod][key][i]) )



        for i in range(len(self.line_names[1:])):
            self.assertListEqual( list(self.res['weighting_res']['pyccf']['downsampled_CCCD'][i]), list(self.loaded_res['weighting_res']['pyccf']['downsampled_CCCD'][i]) )
            self.assertListEqual( list(self.res['weighting_res']['javelin']['downsampled_lag_dist'][i]), list(self.loaded_res['weighting_res']['javelin']['downsampled_lag_dist'][i]) )

            self.assertListEqual( list(self.res['weighting_res']['pyccf']['centroid'][i]), list(self.loaded_res['weighting_res']['pyccf']['centroid'][i]) )
            self.assertListEqual( list(self.res['weighting_res']['javelin']['tophat_lag'][i]), list(self.loaded_res['weighting_res']['javelin']['tophat_lag'][i]) )


        for mod in ['pyccf', 'javelin']:
            self.assertListEqual( list(self.res['weighting_res'][mod]['frac_rejected']), list(self.loaded_res['weighting_res'][mod]['frac_rejected']) )



    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
