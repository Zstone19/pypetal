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

        pyccf_params = {
            'nsim': 100,
            'interp': 1.75,
            'sigmode': .25,
            'thres': 0.79
        }

        weighting_params = {
            'k': 3,
            'gap_size': 20
        }

        lag_bounds = [-500, 500]

        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                            run_pyccf=True, pyccf_params=pyccf_params,
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
        self.assertEquals( len(self.loaded_res['weighting_res']['javelin'].keys()), 0 )

        self.assertListEqual( self.line_names[1:], list(self.loaded_res['weighting_res']['names']) )

        weight_keys1 = ['centroid', 'bounds', 'acf', 'lags', 'weight_dist',
                        'smoothed_dist', 'ntau', 'downsampled_CCCD', 'frac_rejected']

        for key in weight_keys1:
            self.assertIn(key, self.loaded_res['weighting_res']['pyccf'].keys())

        self.assertNotIn( 'CCCD', self.loaded_res['weighting_res']['pyccf'].keys() )

        for key in ['bounds', 'acf', 'lags', 'weight_dist', 'smoothed_dist', 'ntau']:
            for i in range(len(self.line_names[1:])):
                self.assertListEqual( list(self.res['weighting_res']['pyccf'][key][i]), list(self.loaded_res['weighting_res']['pyccf'][key][i]) )

        for i in range(len(self.line_names[1:])):
            self.assertListEqual( list(self.res['weighting_res']['pyccf']['downsampled_CCCD'][i]), list(self.loaded_res['weighting_res']['pyccf']['downsampled_CCCD'][i]) )
            self.assertListEqual( list(self.res['weighting_res']['pyccf']['centroid'][i]), list(self.loaded_res['weighting_res']['pyccf']['centroid'][i]) )

        self.assertListEqual( list(self.res['weighting_res']['pyccf']['frac_rejected']), list(self.loaded_res['weighting_res']['pyccf']['frac_rejected']) )



    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
