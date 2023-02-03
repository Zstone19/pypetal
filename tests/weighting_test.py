import pypetal.pipeline as pl
from pypetal.load import read_weighting_summary
from pypetal.weighting import prob_tau
import numpy as np 

import glob
import os
import shutil

import unittest

class TestWeighting(unittest.TestCase):

    def setUp(self):

        main_dir = 'examples/dat/javelin_'            
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']


        if os.path.exists(output_dir) and os.path.isdir(output_dir):
            shutil.rmtree(output_dir)


        pyccf_params = {
            'nsim': 300,
            'interp': 1.5
        }

        javelin_params = {
            'nchain': 50,
            'nburn': 20,
            'nwalker': 20
        }

        weighting_params = {
            'k': 3,
            'gap_size': 20
        }

        lag_bounds = [ [-1000, 1000], [-500, 500] ]

        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                            run_pyccf=True, pyccf_params=pyccf_params,
                            run_javelin=True, javelin_params=javelin_params, 
                            run_weighting=True, weighting_params=weighting_params,
                            lag_bounds=lag_bounds,
                            file_fmt='ascii',
                            threads=45)

        self.filenames = filenames
        self.line_names = line_names
        self.res = res




    #Make sure the lengths and keys of each of the resulting arrays are correct
    def test_res(self):
        
        pyccf_keys = ['centroid', 'bounds', 'acf', 'lags', 'weight_dist', 
                        'smoothed_dist', 'ntau', 'downsampled_CCCD', 'frac_rejected']
        javelin_keys = ['tophat_lag', 'bounds', 'acf', 'lags', 'weight_dist', 
                        'smoothed_dist', 'ntau', 'downsampled_lag_dist', 'frac_rejected']

        #Make sure the keys of the output are correct
        for module, expected_keys in zip(['pyccf', 'javelin'], [pyccf_keys, javelin_keys]):
            res_keys = list(self.res['weighting_res'][module].keys())
            for key in expected_keys:
                self.assertIn( key, res_keys )

        self.assertIn( 'rmax', list(self.res['weighting_res'].keys()) )


        array_mask = [True, True, True, True, True, True, True, True, False]

        #Make sure data types are correct
        for i in range(len(self.filenames)-1):

            for module, expected_keys in zip(['pyccf', 'javelin'], [pyccf_keys, javelin_keys]):
                for j, key in enumerate(expected_keys):

                    if array_mask[j]:
                        self.assertIs( type(self.res['weighting_res'][module][key][i][0]), np.float64 )
                    else:
                        self.assertIs( type(self.res['weighting_res'][module][key][i]), float )

        
        #Make sure lengths of arrays are correct
        for key in pyccf_keys:
            self.assertEqual( len(self.res['weighting_res']['pyccf'][key]), len(self.filenames)-1 )

        for key in javelin_keys:
            self.assertEqual( len(self.res['weighting_res']['javelin'][key]), len(self.filenames)-1 )

        for i in range(len(self.filenames)-1):
            self.assertEqual( len(self.res['weighting_res']['pyccf']['centroid'][i]), 3 )
            self.assertEqual( len(self.res['weighting_res']['pyccf']['bounds'][i]), 3 )

            self.assertEqual( len(self.res['weighting_res']['javelin']['tophat_lag'][i]), 3 )
            self.assertEqual( len(self.res['weighting_res']['javelin']['bounds'][i]), 3 )


        #Make sure the reported lag is the med and error for the downsampled dist,
        #frac rejected is correct, ntau is correct

        xcont, _, _ = np.loadtxt(self.filenames[0], unpack=True, usecols=[0,1,2])
        for i in range(1, len(self.filenames)):
            xline, _, _ = np.loadtxt(self.filenames[i], unpack=True, usecols=[0,1,2])

            cccd = self.res['pyccf_res'][i-1]['CCCD_lags']
            lag_dist = self.res['javelin_res'][i-1]['tophat_params'][0]

            cccd_ds = self.res['weighting_res']['pyccf']['downsampled_CCCD'][i-1]
            lag_dist_ds = self.res['weighting_res']['javelin']['downsampled_lag_dist'][i-1]


            #Lag vals
            self.assertEqual( self.res['weighting_res']['pyccf']['centroid'][i-1][1], np.median(cccd_ds) )
            self.assertEqual( self.res['weighting_res']['pyccf']['centroid'][i-1][0], np.median(cccd_ds) - np.percentile(cccd_ds, 16) )
            self.assertEqual( self.res['weighting_res']['pyccf']['centroid'][i-1][2], np.percentile(cccd_ds, 84) - np.median(cccd_ds) )

            self.assertEqual( self.res['weighting_res']['javelin']['tophat_lag'][i-1][1], np.median(lag_dist_ds) )
            self.assertEqual( self.res['weighting_res']['javelin']['tophat_lag'][i-1][0], np.median(lag_dist_ds) - np.percentile(lag_dist_ds, 16) )
            self.assertEqual( self.res['weighting_res']['javelin']['tophat_lag'][i-1][2], np.percentile(lag_dist_ds, 84) - np.median(lag_dist_ds) )

            #Frac rejected
            self.assertEqual( self.res['weighting_res']['pyccf']['frac_rejected'][i-1], 1- len(cccd_ds) / len(cccd) )
            self.assertEqual( self.res['weighting_res']['javelin']['frac_rejected'][i-1], 1- len(lag_dist_ds) / len(lag_dist) )

            #Ntau
            pyccf_lags = self.res['weighting_res']['pyccf']['lags'][i-1]
            jav_lags = self.res['weighting_res']['javelin']['lags'][i-1]

            _, pyccf_ntau, pyccf_n0, _ = prob_tau(xcont, xline, lagvals=pyccf_lags, gap_size=20, k=3)
            _, jav_ntau, jav_n0, _ = prob_tau(xcont, xline, lagvals=jav_lags, gap_size=20, k=3)

            self.assertListEqual( list(self.res['weighting_res']['pyccf']['ntau'][i-1]), list(pyccf_ntau) )
            self.assertListEqual( list(self.res['weighting_res']['javelin']['ntau'][i-1]), list(jav_ntau) )



    #Make sure the layout of the files is correct
    def test_file_layout(self):
        main_directories = glob.glob('.tmp/*/')
        subdirs = ['.tmp/continuum/', '.tmp/yelm/', '.tmp/zing/']

        #Make sure each line has a subdirectory
        self.assertIn( '.tmp/continuum/' , main_directories )
        self.assertIn( '.tmp/yelm/' , main_directories )
        self.assertIn( '.tmp/zing/' , main_directories )

        #Make sure eachline has a "pyccf", "javelin", and "weights" subdir
        for fdir, name in zip(subdirs[1:], self.line_names[1:]):
            for subdir in ['pyccf', 'javelin', 'weights']:
                self.assertIn( '.tmp/' + name + '/' + subdir, glob.glob( fdir + '*' ) )
                
        #Make sure continuum has no files
        self.assertEqual( len(glob.glob('.tmp/continuum/*')), 0 )

        #Make sure plots were output
        self.assertIn( '.tmp/javelin_weights_res.pdf', glob.glob('.tmp/*') )
        self.assertIn( '.tmp/pyccf_weights_res.pdf', glob.glob('.tmp/*') )


        #Make sure each "weights" subdir has the proper files
        for i, fdir in zip([1,2], subdirs[1:]):
            files = glob.glob(fdir + 'weights/*')

            self.assertIn( fdir + 'weights/javelin_weighted_lag_dist.dat', files )
            self.assertIn( fdir + 'weights/javelin_weights.dat', files )
            self.assertIn( fdir + 'weights/pyccf_weighted_cccd.dat', files )
            self.assertIn( fdir + 'weights/pyccf_weights.dat', files )
            self.assertIn( fdir + 'weights/weight_summary.txt', files )

            self.assertIn( fdir + 'weights/' + self.line_names[i] + '_weights.pdf', files )

        


    def test_res_file_match(self):

        #Match output files
        for i, name in zip([1,2], self.line_names[1:]):
            lag_dist_ds = np.loadtxt('.tmp/' + name + '/weights/javelin_weighted_lag_dist.dat')
            cccd_ds = np.loadtxt('.tmp/' + name + '/weights/pyccf_weighted_cccd.dat')

            self.assertEqual( list(self.res['weighting_res']['javelin']['downsampled_lag_dist'][i-1]), list(lag_dist_ds) )
            self.assertEqual( list(self.res['weighting_res']['pyccf']['downsampled_CCCD'][i-1]), list(cccd_ds) )



            lags, ntau, weight_dist, acf, smooth_dist, smooth_weight_dist = np.loadtxt( '.tmp/' + name + '/weights/javelin_weights.dat', unpack=True, delimiter=',' )
            self.assertListEqual( list(self.res['weighting_res']['javelin']['lags'][i-1]), list(lags) )
            self.assertListEqual( list(self.res['weighting_res']['javelin']['ntau'][i-1]), list(ntau) )
            self.assertListEqual( list(self.res['weighting_res']['javelin']['weight_dist'][i-1]), list(weight_dist) )
            self.assertListEqual( list(self.res['weighting_res']['javelin']['acf'][i-1]), list(acf) )
            self.assertListEqual( list(self.res['weighting_res']['javelin']['smoothed_dist'][i-1]), list(smooth_weight_dist) )

            lags, ntau, weight_dist, acf, smooth_dist, smooth_weight_dist = np.loadtxt( '.tmp/' + name + '/weights/pyccf_weights.dat', unpack=True, delimiter=',' )
            self.assertListEqual( list(self.res['weighting_res']['pyccf']['lags'][i-1]), list(lags) )
            self.assertListEqual( list(self.res['weighting_res']['pyccf']['ntau'][i-1]), list(ntau) )
            self.assertListEqual( list(self.res['weighting_res']['pyccf']['weight_dist'][i-1]), list(weight_dist) )
            self.assertListEqual( list(self.res['weighting_res']['pyccf']['acf'][i-1]), list(acf) )
            self.assertListEqual( list(self.res['weighting_res']['pyccf']['smoothed_dist'][i-1]), list(smooth_weight_dist) )       



        #Match summary files
        xcont, _, _ = np.loadtxt(self.filenames[0], unpack=True, usecols=[0,1,2])
        for i, name in zip([1,2], self.line_names[1:]):
            weight_dict = read_weighting_summary( '.tmp/' + name + '/weights/weight_summary.txt' )
            xline, _, _ = np.loadtxt(self.filenames[i], unpack=True, usecols=[0,1,2])

            #k
            self.assertEqual( weight_dict['k'], 3 )

            #Get N0
            pyccf_lags = self.res['weighting_res']['pyccf']['lags'][i-1]
            jav_lags = self.res['weighting_res']['javelin']['lags'][i-1]
            _, _, pyccf_n0, _ = prob_tau(xcont, xline, lagvals=pyccf_lags, gap_size=20, k=3)
            _, _, jav_n0, _ = prob_tau(xcont, xline, lagvals=jav_lags, gap_size=20, k=3)

            self.assertEqual( weight_dict['pyccf_n0'], pyccf_n0 )
            self.assertEqual( weight_dict['javelin_n0'], jav_n0 )


            #Bounds
            pyccf_bounds = [ self.res['weighting_res']['pyccf']['bounds'][i-1][0], self.res['weighting_res']['pyccf']['bounds'][i-1][2] ]
            pyccf_peak = self.res['weighting_res']['pyccf']['bounds'][i-1][1]

            jav_bounds = [ self.res['weighting_res']['javelin']['bounds'][i-1][0], self.res['weighting_res']['javelin']['bounds'][i-1][2] ]
            jav_peak = self.res['weighting_res']['javelin']['bounds'][i-1][1]

            self.assertListEqual( list(weight_dict['pyccf_lag_bounds']), list(pyccf_bounds) )
            self.assertEqual( weight_dict['pyccf_peak'], pyccf_peak )
            self.assertListEqual( list(weight_dict['javelin_lag_bounds']), list(jav_bounds) )
            self.assertEqual( weight_dict['javelin_peak'], jav_peak )


            #Lag
            pyccf_lag_err = [ self.res['weighting_res']['pyccf']['centroid'][i-1][0], self.res['weighting_res']['pyccf']['centroid'][i-1][2] ]
            pyccf_lag = self.res['weighting_res']['pyccf']['centroid'][i-1][1]

            jav_lag_err = [ self.res['weighting_res']['javelin']['tophat_lag'][i-1][0], self.res['weighting_res']['javelin']['tophat_lag'][i-1][2] ]
            jav_lag = self.res['weighting_res']['javelin']['tophat_lag'][i-1][1]

            self.assertListEqual( list(weight_dict['pyccf_lag_uncertainty']), list(pyccf_lag_err) )
            self.assertEqual( weight_dict['pyccf_lag_value'], pyccf_lag )
            self.assertListEqual( list(weight_dict['javelin_lag_uncertainty']), list(jav_lag_err) )
            self.assertEqual( weight_dict['javelin_lag_value'], jav_lag )


            #Rmax
            self.assertEqual( weight_dict['rmax'], self.res['weighting_res']['rmax'][i-1] )


            #Frac rejected
            cccd = self.res['weighting_res']['pyccf']['CCCD'][i-1]
            lag_dist = self.res['weighting_res']['javelin']['lag_dist'][i-1]

            cccd_ds = self.res['weighting_res']['pyccf']['downsampled_CCCD'][i-1]
            lag_dist_ds = self.res['weighting_res']['javelin']['downsampled_lag_dist'][i-1]

            self.assertEqual( weight_dict['pyccf_frac_rejected'], 1- len(cccd_ds)/len(cccd) )
            self.assertEqual( weight_dict['javelin_frac_rejected'], 1- len(lag_dist_ds)/len(lag_dist) )


    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')

