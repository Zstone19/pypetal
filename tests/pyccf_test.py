import pypetal.pipeline as pl
import numpy as np 

import glob
import os
import shutil

import unittest

class TestPyCCF(unittest.TestCase):

    def setUp(self):

        main_dir = '../examples/dat/javelin_'
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']


        if os.path.exists(output_dir) and os.path.isdir(output_dir):
            shutil.rmtree(output_dir)


        params = {
            'nsim': 100,
            'interp': 1.75,
            'sigmode': .25,
            'thres': 0.79
        }
        lag_bounds = [ [-1000, 1000], [-500, 500] ]

        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                            run_pyccf=True,
                            pyccf_params=params, 
                            lag_bounds=lag_bounds,
                            file_fmt='ascii')

        self.filenames = filenames
        self.line_names = line_names
        self.res = res

    #Make sure the lengths and keys of each of the resulting arrays are correct
    def test_res(self):
        
        expected_keys = ['CCF', 'CCF_lags', 'centroid', 'centroid_err_lo', 'centroid_err_hi',
                            'peak', 'peak_err_lo', 'peak_err_hi', 'CCCD_lags', 'CCPD_lags', 'name']

        num_dicts = len(self.line_names) - 1
        self.assertEqual( num_dicts, len(self.res['pyccf_res']) )

        mc_length = 100

        #Make sure the keys of the output are correct
        for i in range(num_dicts):
            res_keys = list(self.res['pyccf_res'][i].keys())
            for key in expected_keys:
                self.assertIn( key, res_keys )


        #Make sure the names in the dicts match the line names
        for i in range(num_dicts):
            self.assertEqual( self.res['pyccf_res'][i]['name'], self.line_names[i+1] )

        #Make sure data types are correct
        array_mask = [True, True, False, False, False, False, False, False, True, True, False]
        data_type = [np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, str]
        for i in range(num_dicts):
    
            for j, key in enumerate(expected_keys):
                if array_mask[j]:
                    self.assertIs( type(self.res['pyccf_res'][i][key][0]), data_type[j] )
                else:
                    self.assertIs( type(self.res['pyccf_res'][i][key]), data_type[j] )

        
        #Make sure lengths of arrays are correct
        for i in range(num_dicts):
            self.assertEqual( len(self.res['pyccf_res'][i]['CCCD_lags']), mc_length )
            self.assertEqual( len(self.res['pyccf_res'][i]['CCPD_lags']), mc_length )


        #Make sure that the reported centroid and peak are the medians
        for i in range(num_dicts):
            self.assertEqual( np.median(self.res['pyccf_res'][i]['CCCD_lags']), self.res['pyccf_res'][i]['centroid'] )
            self.assertEqual( np.median(self.res['pyccf_res'][i]['CCPD_lags']), self.res['pyccf_res'][i]['peak'] )


        lag_bounds = [ [-1000, 1000], [-500, 500] ]
        #Make sure lag bounds worked
        for i in range(num_dicts):

            self.assertGreaterEqual( np.min(self.res['pyccf_res'][i]['CCF_lags']), lag_bounds[i][0] )
            self.assertLessEqual( np.max(self.res['pyccf_res'][i]['CCF_lags']), lag_bounds[i][1] )

            self.assertGreaterEqual( np.min(self.res['pyccf_res'][i]['CCCD_lags']), lag_bounds[i][0] )
            self.assertLessEqual( np.max(self.res['pyccf_res'][i]['CCCD_lags']), lag_bounds[i][1] )

            self.assertGreaterEqual( np.min(self.res['pyccf_res'][i]['CCPD_lags']), lag_bounds[i][0] )
            self.assertLessEqual( np.max(self.res['pyccf_res'][i]['CCPD_lags']), lag_bounds[i][1] )




    #Make sure the layout of the files is correct
    def test_file_layout(self):
        main_directories = glob.glob('.tmp/*/')
        subdirs = ['.tmp/continuum/', '.tmp/yelm/', '.tmp/zing/']

        mc_length = 100

        #Make sure each line has a subdirectory
        self.assertIn( '.tmp/continuum/' , main_directories )
        self.assertIn( '.tmp/yelm/' , main_directories )
        self.assertIn( '.tmp/zing/' , main_directories )

        #Make sure continuum has no files
        self.assertEqual( len(glob.glob('.tmp/continuum/*')), 0 )

        #Make sure each line has a "pyccf" subdirectory
        for fdir in subdirs[1:]:
            subdir_dirs = glob.glob(fdir + '*')
            self.assertIn( fdir + 'pyccf' , subdir_dirs )

        #Make sure each "pyccf" subdirectory has a ccf and dist file
        for i, fdir in zip([1,2], subdirs[1:]):
            files = glob.glob(fdir + 'pyccf/*')
            self.assertIn( fdir + 'pyccf/' + self.line_names[i] + '_ccf.dat' , files )
            self.assertIn( fdir + 'pyccf/' + self.line_names[i] + '_ccf_dists.dat' , files )
            self.assertIn( fdir + 'pyccf/' + self.line_names[i] + '_ccf.pdf' , files )

            #Make sure the dist file has the correct number of lines
            cccd, ccpd = np.loadtxt(fdir + 'pyccf/' + self.line_names[i] + '_ccf_dists.dat', unpack=True, delimiter=',', usecols=[0,1])
            self.assertEqual( len(cccd), mc_length )


        #Make sure the light curves get saved
        self.assertIn( '.tmp/light_curves/', main_directories )

        lc_files = glob.glob('.tmp/light_curves/*')
        for name in self.line_names:
            self.assertIn( '.tmp/light_curves/' + name + '.dat', lc_files )

        


    def test_res_file_match(self):

        #Match CCF
        for i, name in enumerate(self.line_names[1:]):
            file_lags, file_ccf = np.loadtxt('.tmp/' + name + '/pyccf/' + name + '_ccf.dat', unpack=True, delimiter=',')
            res_lags, res_ccf = self.res['pyccf_res'][i]['CCF_lags'], self.res['pyccf_res'][i]['CCF']

            self.assertListEqual( list(file_lags), list(res_lags) )
            self.assertListEqual( list(file_ccf), list(res_ccf) )


        #Match dists
        for i, name in enumerate(self.line_names[1:]):
            file_cccd, file_ccpd = np.loadtxt('.tmp/' + name + '/pyccf/' + name + '_ccf_dists.dat', unpack=True, delimiter=',')
            res_cccd, res_ccpd = self.res['pyccf_res'][i]['CCCD_lags'], self.res['pyccf_res'][i]['CCPD_lags']

            self.assertListEqual( list(file_cccd), list(res_cccd) )
            self.assertListEqual( list(file_ccpd), list(res_ccpd) )


    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')

