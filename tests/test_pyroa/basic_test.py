import glob
import os
import shutil
import unittest

import numpy as np
from PyROA import Fit

import pypetal.pipeline as pl
from pypetal.pyroa.utils import get_samples_chunks


class TestPyROA(unittest.TestCase):

    def setUp(self):

        main_dir = 'examples/dat/javelin_'
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']

        params = {
            'nchain': 1500,
            'nburn': 1000,
            'add_var': False,
            'delay_dist': False,
            'subtract_mean': True,
            'together': True
        }
        lag_bounds = [ [-1000, 1000], [-500, 500] ]

        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                            run_pyroa=True,
                            pyroa_params=params,
                            lag_bounds=lag_bounds,
                            file_fmt='ascii')

        self.filenames = filenames
        self.line_names = line_names
        self.res = res


    def test_all(self):

        #################################################################################################
        # RES
        #################################################################################################
        #Make sure the lengths and keys of each of the resulting arrays are correct

        self.assertIs( isinstance( self.res['pyroa_res'], Fit ), True )
        mc_length = 1500
        nparams = 2 + 3*2

        #Make sure lengths of arrays are correct
        self.assertEqual( self.res['pyroa_res'].samples.shape, (nparams, mc_length) )
        self.assertEqual( len(self.res['pyroa_res'].models), 3 )


        #################################################################################################
        # FILE LAYOUT
        #################################################################################################
        #Make sure the layout of the files is correct

        main_directories = glob.glob('.tmp/*/')
        subdirs = ['.tmp/continuum/', '.tmp/yelm/', '.tmp/zing/']

        mc_length = 1500

        #Make sure each line has a subdirectory
        self.assertIn( '.tmp/continuum/' , main_directories )
        self.assertIn( '.tmp/yelm/' , main_directories )
        self.assertIn( '.tmp/zing/' , main_directories )

        #Make sure continuum has no files
        self.assertEqual( len(glob.glob('.tmp/continuum/*')), 0 )

        #Make sure each line doesn't have a "pyroa" subdirectory
        for fdir in subdirs[1:]:
            subdir_dirs = glob.glob(fdir + '*')
            self.assertNotIn( fdir + 'pyroa' , subdir_dirs )

        #Make sure there is one main pyroa directory
        self.assertIn( '.tmp/pyroa/', main_directories )

        #Make sure "pyroa" subdirectory has PyROA files
        files = glob.glob('.tmp/pyroa/*')
        self.assertIn( '.tmp/pyroa/samples.obj' , files )
        self.assertIn( '.tmp/pyroa/samples_flat.obj' , files )
        self.assertIn( '.tmp/pyroa/Lightcurve_models.obj' , files )
        self.assertIn( '.tmp/pyroa/X_t.obj' , files )

        #Make sure "pyroa" subdirectory has PyROA plots
        files = glob.glob('.tmp/pyroa/*')
        self.assertIn( '.tmp/pyroa/samples.obj' , files )
        self.assertIn( '.tmp/pyroa/samples_flat.obj' , files )
        self.assertIn( '.tmp/pyroa/Lightcurve_models.obj' , files )
        self.assertIn( '.tmp/pyroa/X_t.obj' , files )

        #Make sure PyROA light curves get saved
        self.assertIn( '.tmp/pyroa_lcs/', main_directories )
        lc_files = glob.glob('.tmp/pyroa_lcs/*')
        for name in self.line_names:
            self.assertIn( '.tmp/pyroa_lcs/pyroa_' + name + '.dat', files )


        #Make sure the light curves get saved
        self.assertIn( '.tmp/light_curves/', main_directories )

        lc_files = glob.glob('.tmp/light_curves/*')
        for name in self.line_names:
            self.assertIn( '.tmp/light_curves/' + name + '.dat', lc_files )



        #################################################################################################
        # FUNCTIONS
        #################################################################################################

        #get_samples_chunks
        samples = self.res['pyroa_res'].samples
        samples_chunks = get_samples_chunks(samples, 1000)

        for i in range(3):
            self.assertEqual( len(samples_chunks[i]), 3 )
            for j in range(3):
                self.assertEqual( len(samples_chunks[i][j]), 500 )

        self.assertEqual( len(samples_chunks[-1]), 1 )
        self.assertEqual( len(samples_chunks[-1][0]), 500 )


        #################################################################################################
        # RES FILE MATCH
        #################################################################################################

        #Make sure the pyroa_lcs are equal to the original ones (with the mean subtracted)
        for i in range(len(self.line_names)):
            x1, y1, yerr1 = np.loadtxt( self.filenames[i], unpack=True, usecols=[0,1,2] )
            x2, y2, yerr2 = np.loadtxt( '.tmp/pyroa_lcs/pyroa_' + self.line_names[i] + '.dat', unpack=True, usecols=[0,1,2] )

            self.assertTrue( np.allclose(x1, x2) )
            self.assertTrue( np.allclose(y1 - np.mean(y1), y2) )
            self.assertTrue( np.allclose(yerr1, yerr2) )


    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
