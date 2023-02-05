import glob
import os
import shutil
import unittest

import numpy as np

import pypetal.pipeline as pl
from pypetal.load import read_weighting_summary
from pypetal.weighting import prob_tau


class TestWeighting(unittest.TestCase):

    def setUp(self):

        main_dir = 'examples/dat/javelin_'
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']

        detrend_params = {
            'K': 2
        }

        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                            run_detrend=True, detrend_params=detrend_params,
                            file_fmt='ascii')

        self.filenames = filenames
        self.line_names = line_names
        self.res = res

    #Make sure the layout of the files is correct
    def test_file_layout(self):

        main_directories = glob.glob('.tmp/*/')
        subdirs = ['.tmp/continuum/', '.tmp/yelm/', '.tmp/zing/']

        #Make sure each line has a subdirectory
        self.assertIn( '.tmp/continuum/' , main_directories )
        self.assertIn( '.tmp/yelm/' , main_directories )
        self.assertIn( '.tmp/zing/' , main_directories )

        #Make sure each has a detrend plot
        self.assertIn( '.tmp/continuum/detrend.pdf', glob.glob('.tmp/continuum/*') )
        self.assertIn( '.tmp/yelm/detrend.pdf', glob.glob('.tmp/yelm/*') )
        self.assertIn( '.tmp/zing/detrend.pdf', glob.glob('.tmp/zing/*') )


        #Make sure the light curves get saved
        self.assertIn( '.tmp/light_curves/', main_directories )

        lc_files = glob.glob('.tmp/light_curves/*')
        for name in self.line_names:
            self.assertIn( '.tmp/light_curves/' + name + '.dat', lc_files )


        #Make sure the detrended light curves get saved
        self.assertIn( '.tmp/processed_lcs/', main_directories )

        lc_files = glob.glob('.tmp/processed_lcs/*')
        for i, name in enumerate(self.line_names):
            self.assertIn( '.tmp/processed_lcs/' + name + '_data.dat', lc_files )

            #Make sure the processed lcs have same length as original
            x, _, _ = np.loadtxt('.tmp/processed_lcs/' + name + '_data.dat', unpack=True, delimiter=',', usecols=[0,1,2])
            xdata, _, _ = np.loadtxt(self.filenames[i], unpack=True, usecols=[0,1,2])
            self.assertEqual( len(x), len(xdata) )

    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
