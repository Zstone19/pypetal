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
            'jitter': False
        }

        #Run pypetal
        res = pl.run_pipeline(self.output_dir, self.filenames, self.line_names,
                            run_drw_rej=True, drw_rej_params=drw_rej_params,
                            file_fmt='ascii',
                            time_unit='d',
                            lc_unit='mag')

        reject_data = drw_rej_params['reject_data']
        loaded_res = load('.tmp/')


        ###############################################################################################
        # MATCH
        ###############################################################################################

        res_keys = list(res['drw_rej_res'].keys())
        load_keys = list(loaded_res['drw_rej_res'].keys())

        for key in res_keys:
            if key != 'jitters':
                self.assertIn(key, load_keys)


        for key in ['fit_x', 'fit_y', 'fit_err', 'names']:
            self.assertIn(key, load_keys)

        self.assertNotIn( 'jitters', res_keys )


        #Names
        self.assertListEqual( loaded_res['drw_rej_res']['names'], self.line_names )

        #Reject_data
        self.assertListEqual( loaded_res['drw_rej_res']['reject_data'].tolist(), res['drw_rej_res']['reject_data'] )

        #Chains and masks
        for i in range(len(self.filenames)):

            if not reject_data[i]:
                self.assertIsNone( loaded_res['drw_rej_res']['taus'][i] )
                self.assertIsNone( loaded_res['drw_rej_res']['sigmas'][i] )

                self.assertIsNone( res['drw_rej_res']['taus'][i] )
                self.assertIsNone( res['drw_rej_res']['sigmas'][i] )
                continue

            self.assertListEqual( list(res['drw_rej_res']['taus'][i]), list(loaded_res['drw_rej_res']['taus'][i]) )
            self.assertListEqual( list(res['drw_rej_res']['sigmas'][i]), list(loaded_res['drw_rej_res']['sigmas'][i]) )
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


    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
