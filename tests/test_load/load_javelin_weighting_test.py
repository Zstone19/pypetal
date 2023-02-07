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

        javelin_params = {
            'nchain': 40,
            'nburn': 20,
            'nwalker': 20
        }

        weighting_params = {
            'k': 3,
            'gap_size': 20
        }

        lag_bounds = [-500,500]

        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                            run_javelin=True, javelin_params=javelin_params,
                            run_weighting=True, weighting_params=weighting_params,
                            lag_bounds=lag_bounds,
                            file_fmt='ascii',
                            time_unit='d',
                            lc_unit=['Jy', 'mJy', 'mJy'],
                            threads=2)

        self.filenames = filenames
        self.line_names = line_names
        self.res = res

        self.loaded_res = load('.tmp/')


    def test_match(self):

        res_keys = list(self.res['weighting_res'].keys())
        load_keys = list(self.loaded_res['weighting_res'].keys())

        self.assertIn( 'names', load_keys )
        for key in res_keys:
            self.assertIn(key, load_keys)

        self.assertIn('rmax', load_keys)
        self.assertEquals( len(self.loaded_res['weighting_res']['rmax']), 0 )
        self.assertEquals( len(self.loaded_res['weighting_res']['pyccf'].keys()), 0 )

        self.assertListEquals( self.line_names, list(self.loaded_res['weighting_res']['names']) )

        weight_keys2 = ['tophat_lag', 'bounds', 'acf', 'lags', 'weight_dist',
                        'smoothed_dist', 'ntau', 'downsampled_lag_dist', 'frac_rejected']

        for key in weight_keys2:
            self.assertIn(key, self.loaded_res['weighting_res']['javelin'].keys())

        self.assertNotIn( 'lag_dist', self.loaded_res['weighting_res']['javelin'].keys() )

        for key in ['bounds', 'acf', 'lags', 'weight_dist', 'smoothed_dist', 'ntau']:
            for i in range(len(self.line_names[1:])):
                self.assertListEqual( list(self.res['weighting_res']['javelin'][key][i]), list(self.loaded_res['weighting_res']['javelin'][key][i]) )

        for i in range(len(self.line_names[1:])):
            self.assertListEqual( list(self.res['weighting_res']['javelin']['downsampled_lag_dist'][i]), list(self.loaded_res['weighting_res']['javelin']['downsampled_lag_dist'][i]) )
            self.assertListEqual( list(self.res['weighting_res']['javelin']['tophat_lag'][i]), list(self.loaded_res['weighting_res']['javelin']['tophat_lag'][i]) )

        self.assertListEqual( list(self.res['weighting_res']['javelin']['frac_rejected']), list(self.loaded_res['weighting_res']['javelin']['frac_rejected']) )


    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
