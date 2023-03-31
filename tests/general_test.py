import glob
import os
import shutil
import unittest

import numpy as np

import pypetal.pipeline as pl


#Test arg2 with filenames
class TestGeneral1(unittest.TestCase):

    def setUp(self):

        main_dir = 'examples/dat/javelin_'
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']

        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                              file_fmt='ascii')

        self.filenames = filenames
        self.line_names = line_names
        self.res = res


    #Make sure the layout of the files is correct
    def test_file_layout(self):
        main_directories = glob.glob('.tmp/*/')

        #Make sure each line has a subdirectory
        self.assertIn( '.tmp/continuum/' , main_directories )
        self.assertIn( '.tmp/yelm/' , main_directories )
        self.assertIn( '.tmp/zing/' , main_directories )

        #Make sure each line subdir has no files
        self.assertEqual( len(glob.glob('.tmp/continuum/*')), 0 )
        self.assertEqual( len(glob.glob('.tmp/yelm/*')), 0 )
        self.assertEqual( len(glob.glob('.tmp/zing/*')), 0 )


        #Make sure the light curves get saved
        self.assertIn( '.tmp/light_curves/', main_directories )

        lc_files = glob.glob('.tmp/light_curves/*')
        for name in self.line_names:
            self.assertIn( '.tmp/light_curves/' + name + '.dat', lc_files )



        #Make sure the light curves are correct
        for i in range(len(self.filenames)):
            x, y, yerr = np.loadtxt(self.filenames[i], unpack=True, usecols=[0,1,2])
            x_file, y_file, yerr_file = np.loadtxt('.tmp/light_curves/' + self.line_names[i] + '.dat', delimiter=',', unpack=True, usecols=[0,1,2])

            self.assertListEqual( list(x), list(x_file) )
            self.assertListEqual( list(y), list(y_file) )
            self.assertListEqual( list(yerr), list(yerr_file) )


    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')






#Test arg2 with arrays
class TestGeneral2(unittest.TestCase):

    def setUp(self):

        main_dir = 'examples/dat/javelin_'
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)


        arg2 = []
        for i in range(len(filenames)):
            x, y, yerr = np.loadtxt(filenames[i], unpack=True, usecols=[0,1,2])
            arg2.append( [x, y, yerr] )


        #Run pypetal
        res = pl.run_pipeline(output_dir, arg2, line_names,
                              file_fmt='ascii')

        self.filenames = filenames
        self.line_names = line_names
        self.res = res



    #Make sure the layout of the files is correct
    def test_file_layout(self):
        main_directories = glob.glob('.tmp/*/')

        #Make sure each line has a subdirectory
        self.assertIn( '.tmp/continuum/' , main_directories )
        self.assertIn( '.tmp/yelm/' , main_directories )
        self.assertIn( '.tmp/zing/' , main_directories )

        #Make sure each line subdir has no files
        self.assertEqual( len(glob.glob('.tmp/continuum/*')), 0 )
        self.assertEqual( len(glob.glob('.tmp/yelm/*')), 0 )
        self.assertEqual( len(glob.glob('.tmp/zing/*')), 0 )


        #Make sure the light curves get saved
        self.assertIn( '.tmp/light_curves/', main_directories )

        lc_files = glob.glob('.tmp/light_curves/*')
        for name in self.line_names:
            self.assertIn( '.tmp/light_curves/' + name + '.dat', lc_files )



        #Make sure the light curves are correct
        for i in range(len(self.filenames)):
            x, y, yerr = np.loadtxt(self.filenames[i], unpack=True, usecols=[0,1,2])
            x_file, y_file, yerr_file = np.loadtxt('.tmp/light_curves/' + self.line_names[i] + '.dat', delimiter=',', unpack=True, usecols=[0,1,2])

            self.assertListEqual( list(x), list(x_file) )
            self.assertListEqual( list(y), list(y_file) )
            self.assertListEqual( list(yerr), list(yerr_file) )


    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
