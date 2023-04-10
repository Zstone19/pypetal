import glob
import os
import shutil
import unittest

import numpy as np

import pypetal.pipeline as pl
from pypetal.pyroa.utils import MyFit, get_samples_chunks


class TestPyROA(unittest.TestCase):

    def setUp(self):

        main_dir = 'examples/dat/javelin_'
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']

        self.delay_dist = [True, False]
        self.psi_types = ['LogGaussian', 'TruncGaussian']

        params = {
            'nchain': 100,
            'nburn': 50,
            'add_var': True,
            'delay_dist': self.delay_dist,
            'psi_types': self.psi_types,
            'subtract_mean': True,
            'div_mean': True,
            'together': False
        }
        lag_bounds = [-1000, 1000]

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

        mc_length = 100
        nburn = 50


        self.assertEqual( len(self.res['pyroa_res']), len(self.line_names)-1 )
        for i in range(len(self.line_names)-1):

            if self.delay_dist[i]:
                nparams = 3 + 5 + 1
            else:
                nparams = 3 + 4 + 1

            self.assertTrue( isinstance( self.res['pyroa_res'][i], MyFit ) )

            #Make sure lengths of arrays are correct
            self.assertEqual( self.res['pyroa_res'][i].samples.shape[0], mc_length )
            self.assertEqual( self.res['pyroa_res'][i].samples.shape[2], nparams )
            self.assertEqual( len(self.res['pyroa_res'][i].models), 2 )


        #################################################################################################
        # FILE LAYOUT
        #################################################################################################
        #Make sure the layout of the files is correct

        main_directories = glob.glob('.tmp/*/')
        subdirs = ['.tmp/continuum/', '.tmp/yelm/', '.tmp/zing/']

        #Make sure each line has a subdirectory
        self.assertIn( '.tmp/continuum/' , main_directories )
        self.assertIn( '.tmp/yelm/' , main_directories )
        self.assertIn( '.tmp/zing/' , main_directories )

        #Make sure continuum has no files
        self.assertEqual( len(glob.glob('.tmp/continuum/*')), 0 )

        #Make sure each line has a "pyroa" subdirectory
        for fdir in subdirs[1:]:
            subdir_dirs = glob.glob(fdir + '*')
            self.assertIn( fdir + 'pyroa' , subdir_dirs )

        #Make sure there is not a main pyroa directory
        self.assertNotIn( '.tmp/pyroa/', main_directories )


        #Make sure each "pyroa" subdirectory has PyROA files
        for i in range(len(self.line_names)-1):
            files = glob.glob( subdirs[i+1] + 'pyroa/*')
            self.assertIn( subdirs[i+1] + 'pyroa/samples.obj' , files )
            self.assertIn( subdirs[i+1] + 'pyroa/samples_flat.obj' , files )
            self.assertIn( subdirs[i+1] + 'pyroa/Lightcurve_models.obj' , files )
            self.assertIn( subdirs[i+1] + 'pyroa/X_t.obj' , files )

        #Make sure "pyroa" subdirectory has PyROA plots
        for i in range(len(self.line_names)-1):
            files = glob.glob( subdirs[i+1] + 'pyroa/*')
            self.assertIn( subdirs[i+1] + 'pyroa/corner_plot.pdf' , files )
            self.assertIn( subdirs[i+1] + 'pyroa/fits_plot.pdf' , files )
            self.assertIn( subdirs[i+1] + 'pyroa/histogram_plot.pdf' , files )
            self.assertIn( subdirs[i+1] + 'pyroa/trace_plot.pdf' , files )

        #Make sure PyROA light curves get saved
        self.assertIn( '.tmp/pyroa_lcs/', main_directories )

        lc_files = glob.glob('.tmp/pyroa_lcs/*')
        for i in range(len(self.line_names)):
            self.assertIn( '.tmp/pyroa_lcs/pyroa_' + self.line_names[i] + '.dat', lc_files )


        #Make sure the light curves get saved
        self.assertIn( '.tmp/light_curves/', main_directories )

        lc_files = glob.glob('.tmp/light_curves/*')
        for name in self.line_names:
            self.assertIn( '.tmp/light_curves/' + name + '.dat', lc_files )



        #################################################################################################
        # FUNCTIONS
        #################################################################################################

        #get_samples_chunks
        for i in range(len(self.line_names)-1):

            if self.delay_dist[i]:
                nparam = 5
            else:
                nparam = 4

            samples = self.res['pyroa_res'][i].samples
            nwalker = samples.shape[1]

            samples_chunks = get_samples_chunks(samples, nburn, add_var=True, delay_dist=self.delay_dist[i])
            nchain_tot = (mc_length-nburn)*nwalker

            self.assertEqual( len(samples_chunks), 3 )

            for i in range(2):
                self.assertEqual( len(samples_chunks[i]), nparam )
                for j in range(nparam):
                    self.assertEqual( len(samples_chunks[i][j]), nchain_tot )

            self.assertEqual( len(samples_chunks[-1]), 1 )
            self.assertEqual( len(samples_chunks[-1][0]), nchain_tot )


        #################################################################################################
        # RES FILE MATCH
        #################################################################################################

        #Make sure the pyroa_lcs are equal to the original ones (with the mean subtracted)
        for i in range(len(self.line_names)):

            x1, y1, yerr1 = np.loadtxt( self.filenames[i], unpack=True, usecols=[0,1,2] )
            x2, y2, yerr2 = np.loadtxt( '.tmp/pyroa_lcs/pyroa_' + self.line_names[i] + '.dat', unpack=True, usecols=[0,1,2] )

            #div_mean and subtract_mean
            yerr1 /= np.mean(y1)
            y1 /= np.mean(y1)
            y1 -= np.mean(y1)

            self.assertTrue( np.allclose(x1, x2) )
            self.assertTrue( np.allclose(y1, y2) )
            self.assertTrue( np.allclose(yerr1, yerr2) )


    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
