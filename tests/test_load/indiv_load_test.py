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
        self.filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        self.output_dir = '.tmp/'
        self.line_names = ['continuum', 'yelm', 'zing']

        self.plike_dir = 'plike_v4/'



    def test_drw_match(self):

        ###############################################################################################
        # RUN PYPETAL
        ###############################################################################################

        drw_rej_params = {
            'nchain': 100,
            'nburn': 50,
            'nwalker': 20,
            'reject_data': [True, False, True],
            'nsig': 1.0,
            'use_for_javelin': True
        }

        #Run pypetal
        res = pl.run_pipeline(self.output_dir, self.filenames, self.line_names,
                            run_drw_rej=True, drw_rej_params=drw_rej_params,
                            file_fmt='ascii',
                            time_unit='d',
                            lc_unit=['Jy', 'mJy', 'mJy'])

        reject_data = drw_rej_params['reject_data']
        loaded_res = load('.tmp/')


        ###############################################################################################
        # MATCH
        ###############################################################################################

        res_keys = list(res['drw_rej_res'].keys())
        load_keys = list(loaded_res['drw_rej_res'].keys())

        for key in res_keys:
            self.assertIn(key, load_keys)
        for key in ['fit_x', 'fit_y', 'fit_err', 'names']:
            self.assertIn(key, load_keys)


        #Names
        self.assertListEqual( loaded_res['drw_rej_res']['names'], self.line_names )

        #Reject_data
        self.assertListEqual( loaded_res['drw_rej_res']['reject_data'].tolist(), res['drw_rej_res']['reject_data'] )

        #Chains and masks
        for i in range(len(self.filenames)):

            if not reject_data[i]:
                self.assertIsNone( loaded_res['drw_rej_res']['taus'][i] )
                self.assertIsNone( loaded_res['drw_rej_res']['sigmas'][i] )
                self.assertIsNone( loaded_res['drw_rej_res']['jitters'][i] )

                self.assertIsNone( res['drw_rej_res']['taus'][i] )
                self.assertIsNone( res['drw_rej_res']['sigmas'][i] )
                self.assertIsNone( res['drw_rej_res']['jitters'][i] )
                continue

            self.assertListEqual( list(res['drw_rej_res']['taus'][i]), list(loaded_res['drw_rej_res']['taus'][i]) )
            self.assertListEqual( list(res['drw_rej_res']['sigmas'][i]), list(loaded_res['drw_rej_res']['sigmas'][i]) )
            self.assertListEqual( list(res['drw_rej_res']['jitters'][i]), list(loaded_res['drw_rej_res']['jitters'][i]) )
            self.assertListEqual( list(res['drw_rej_res']['masks'][i]), list(loaded_res['drw_rej_res']['masks'][i]) )


        #Fits
        for i, name in enumerate(self.line_names):

            if not reject_data[i]:
                self.assertIsNone( loaded_res['drw_rej_res']['fit_x'][i] )
                self.assertIsNone( loaded_res['drw_rej_res']['fit_y'][i] )
                self.assertIsNone( loaded_res['drw_rej_res']['fit_err'][i] )
                continue

            xfit, yfit, yerrfit = np.loadtxt( '.tmp/' + name + '/drw_rej/' + name + '_drw_fit.dat', unpack=True, delimiter=',' )
            self.assertListEqual( list(xfit), list(loaded_res['drw_rej_res']['fit_x'][i]) )
            self.assertListEqual( list(yfit), list(loaded_res['drw_rej_res']['fit_y'][i]) )
            self.assertListEqual( list(yerrfit), list(loaded_res['drw_rej_res']['fit_err'][i]) )



    def test_pyccf_match(self):

        ###############################################################################################
        # RUN PYPETAL
        ###############################################################################################

        pyccf_params = {
            'nsim': 100,
            'interp': 1.75,
            'sigmode': .25,
            'thres': 0.79
        }

        lag_bounds = [ [-1000,1000], 'baseline' ]

        #Run pypetal
        res = pl.run_pipeline(self.output_dir, self.filenames, self.line_names,
                            run_pyccf=True, pyccf_params=pyccf_params,
                            lag_bounds=lag_bounds,
                            file_fmt='ascii',
                            time_unit='d',
                            lc_unit=['Jy', 'mJy', 'mJy'])

        loaded_res = load('.tmp/')

        ###############################################################################################
        # MATCH
        ###############################################################################################

        for i, name in enumerate(self.line_names[1:]):
            res_keys = list(res['pyccf_res'][i].keys())
            load_keys = list(loaded_res['pyccf_res'][i].keys())


            for key in load_keys:
                self.assertIn(key, res_keys)

                if key == 'name':
                    self.assertEqual(res['pyccf_res'][i][key], loaded_res['pyccf_res'][i][key])
                else:
                    self.assertListEqual(res['pyccf_res'][i][key].tolist(), loaded_res['pyccf_res'][i][key].tolist())






    def test_pyzdcf_match(self):

        ###############################################################################################
        # RUN PYPETAL
        ###############################################################################################

        pyzdcf_params = {
            'nsim': 100,
            'minpts': 12,
            'prefix': 'myname',
            'run_plike': True,
            'plike_dir': self.plike_dir
        }

        lag_bounds = [0,500]

        #Run pypetal
        res = pl.run_pipeline(self.output_dir, self.filenames, self.line_names,
                            run_pyzdcf=True, pyzdcf_params=pyzdcf_params,
                            lag_bounds=lag_bounds,
                            file_fmt='ascii',
                            time_unit='d',
                            lc_unit=['Jy', 'mJy', 'mJy'])

        loaded_res = load('.tmp/')

        ###############################################################################################
        # MATCH PYZDCF
        ###############################################################################################

        for i in range(len(self.filenames[1:])):
            load_keys = list(loaded_res['pyzdcf_res'][i].keys())
            self.assertListEqual(load_keys, ['output', 'name'])

            res_cols = list(res['pyzdcf_res'][i].columns)
            load_cols = list(loaded_res['pyzdcf_res'][i]['output'].columns)
            self.assertListEqual(res_cols, load_cols)

            for col in res_cols:
                if col == 'tau':
                    for j in range(len(res['pyzdcf_res'][i][col])):
                        res_int = int(res['pyzdcf_res'][i][col].iloc[j])
                        load_int = int(loaded_res['pyzdcf_res'][i]['output'][col].iloc[j])
                        self.assertLessEqual( np.abs( res_int - load_int ), 1 )

                elif col in ['-sig(tau)', '+sig(tau)']:
                    for j in range(len(res['pyzdcf_res'][i][col])):
                        res_val = res['pyzdcf_res'][i][col].iloc[j]
                        load_val = loaded_res['pyzdcf_res'][i]['output'][col].iloc[j]
                        self.assertAlmostEqual( res_val, load_val, places=2 )

                elif col in ['dcf', '-err(dcf)', '+err(dcf)']:
                    for j in range(len(res['pyzdcf_res'][i][col])):
                        res_val = res['pyzdcf_res'][i][col].iloc[j]
                        load_val = loaded_res['pyzdcf_res'][i]['output'][col].iloc[j]
                        self.assertAlmostEqual( res_val, load_val, places=4 )

                else:
                    self.assertListEqual( list(res['pyzdcf_res'][i][col]), list(loaded_res['pyzdcf_res'][i]['output'][col]) )


        ###############################################################################################
        # MATCH PLIKE
        ###############################################################################################

        for i in range(len(self.filenames[1:])):
            plike_cols = ['num', 'lag', 'r', '-dr', '+dr', 'likelihood']

            plike_keys = ['output', 'ML_lag', 'ML_lag_err_lo', 'ML_lag_err_hi', 'name']
            res_keys = list(res['plike_res'][i].keys())
            load_keys = list(loaded_res['plike_res'][i].keys())

            self.assertIn( 'name', load_keys )
            for key in plike_keys[:-1]:
                self.assertIn(key, load_keys)
                self.assertIn(key, res_keys)

            self.assertEqual( self.line_names[i+1], loaded_res['plike_res'][i]['name'] )

            for key in plike_keys[1:-1]:
                self.assertEqual( res['plike_res'][i][key], loaded_res['plike_res'][i][key] )

            self.assertListEqual( list(res['plike_res'][i]['output'].colnames), plike_cols )
            self.assertListEqual( list(loaded_res['plike_res'][i]['output'].colnames), plike_cols )
            for col in plike_cols:
                self.assertListEqual( list(res['plike_res'][i]['output'][col]), list(loaded_res['plike_res'][i]['output'][col]) )





    def test_javelin_match(self):

        ###############################################################################################
        # RUN JAVELIN
        ###############################################################################################

        javelin_params = {
            'nchain': 40,
            'nburn': 20,
            'nwalker': 20,
            'rm_type': "phot"
        }

        lag_bounds = [0,500]

        #Run pypetal
        res = pl.run_pipeline(self.output_dir, self.filenames, self.line_names,
                            run_javelin=True, javelin_params=javelin_params,
                            lag_bounds=lag_bounds,
                            file_fmt='ascii',
                            time_unit='d',
                            lc_unit=['Jy', 'mJy', 'mJy'])

        loaded_res = load('.tmp/')

        ###############################################################################################
        # MATCH JAVELIN
        ###############################################################################################

        for i in range(len(self.filenames[1:])):

            jav_keys = ['tau', 'sigma', 'tophat_params',
                        'cont_fit_x', 'cont_fit_y', 'cont_fit_yerr',
                        'fit_x', 'fit_y', 'fit_yerr', 'name']
            res_keys = list(res['javelin_res'][i].keys())
            load_keys = list(loaded_res['javelin_res'][i].keys())

            for key in jav_keys:
                self.assertIn(key, load_keys)

            self.assertEqual( self.line_names[i+1], loaded_res['javelin_res'][i]['name'] )

            for key in ['tau', 'sigma', 'tophat_params']:

                if key in ['tau','sigma']:
                    self.assertListEqual(res['javelin_res'][i][key].tolist(), loaded_res['javelin_res'][i][key].tolist())
                elif key == 'tophat_params':
                    for j in range(4):
                        self.assertListEqual(res['javelin_res'][i][key][j].tolist(), loaded_res['javelin_res'][i][key][j].tolist())


            c_xfit, c_yfit, c_yerrfit = np.loadtxt( '.tmp/' + self.line_names[i+1] + '/javelin/continuum_lc_fits.dat', unpack=True, delimiter=',' )
            xfit, yfit, yerrfit = np.loadtxt( '.tmp/' + self.line_names[i+1] + '/javelin/' + self.line_names[i+1] + '_lc_fits.dat', unpack=True, delimiter=',' )

            self.assertListEqual(c_xfit.tolist(), loaded_res['javelin_res'][i]['cont_fit_x'].tolist() )
            self.assertListEqual(c_yfit.tolist(), loaded_res['javelin_res'][i]['cont_fit_y'].tolist() )
            self.assertListEqual(c_yerrfit.tolist(), loaded_res['javelin_res'][i]['cont_fit_yerr'].tolist() )

            self.assertListEqual(xfit.tolist(), loaded_res['javelin_res'][i]['fit_x'].tolist() )
            self.assertListEqual(yfit.tolist(), loaded_res['javelin_res'][i]['fit_y'].tolist() )
            self.assertListEqual(yerrfit.tolist(), loaded_res['javelin_res'][i]['fit_yerr'].tolist() )



    def test_javelin_together_match(self):

        ###############################################################################################
        # RUN JAVELIN
        ###############################################################################################

        javelin_params = {
            'nchain': 40,
            'nburn': 20,
            'nwalker': 20,
            'rm_type': "spec",
            'together': True
        }

        lag_bounds = 'baseline'

        #Run pypetal
        res = pl.run_pipeline(self.output_dir, self.filenames, self.line_names,
                            run_javelin=True, javelin_params=javelin_params,
                            lag_bounds=lag_bounds,
                            file_fmt='ascii',
                            time_unit='d',
                            lc_unit=['Jy', 'mJy', 'mJy'])

        loaded_res = load('.tmp/')

        ###############################################################################################
        # MATCH JAVELIN
        ###############################################################################################

        jav_keys = ['tau', 'sigma', 'tophat_params',
                    'cont_fit_x', 'cont_fit_y', 'cont_fit_yerr']
        res_keys = list(res['javelin_res'].keys())
        load_keys = list(loaded_res['javelin_res'].keys())

        for key in jav_keys:
            self.assertIn(key, load_keys)
            for name in self.line_names[1:]:
                self.assertIn( name + '_fit_x', load_keys )
                self.assertIn( name + '_fit_y', load_keys )
                self.assertIn( name + '_fit_yerr', load_keys )


        for key in ['tau', 'sigma', 'tophat_params']:

            if key in ['tau','sigma']:
                self.assertListEqual(res['javelin_res'][key].tolist(), loaded_res['javelin_res'][key].tolist())
            elif key == 'tophat_params':
                for j in range(6):
                    self.assertListEqual(res['javelin_res'][key][j].tolist(), loaded_res['javelin_res'][key][j].tolist())


        c_xfit, c_yfit, c_yerrfit = np.loadtxt( '.tmp/javelin/continuum_lc_fits.dat', unpack=True, delimiter=',' )
        self.assertListEqual(c_xfit.tolist(), loaded_res['javelin_res']['cont_fit_x'].tolist() )
        self.assertListEqual(c_yfit.tolist(), loaded_res['javelin_res']['cont_fit_y'].tolist() )
        self.assertListEqual(c_yerrfit.tolist(), loaded_res['javelin_res']['cont_fit_yerr'].tolist() )


        for i in range(1, len(self.line_names)):
            xfit, yfit, yerrfit = np.loadtxt( '.tmp/javelin/' + self.line_names[i] + '_lc_fits.dat', unpack=True, delimiter=',' )
            self.assertListEqual(xfit.tolist(), loaded_res['javelin_res'][ self.line_names[i] + '_fit_x'].tolist() )
            self.assertListEqual(yfit.tolist(), loaded_res['javelin_res'][ self.line_names[i] + '_fit_y'].tolist() )
            self.assertListEqual(yerrfit.tolist(), loaded_res['javelin_res'][ self.line_names[i] + '_fit_yerr'].tolist() )



    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
