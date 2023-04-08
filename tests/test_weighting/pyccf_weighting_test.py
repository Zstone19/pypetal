import glob
import os
import shutil
import unittest

import numpy as np
from astropy.table import Table

import pypetal.pipeline as pl
from pypetal.weighting.utils import prob_tau


class TestWeighting(unittest.TestCase):

    def setUp(self):

        main_dir = 'examples/dat/javelin_'
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']

        pyccf_params = {
            'nsim': 100,
            'interp': 1.5
        }

        weighting_params = {
            'k': 3,
            'gap_size': 20
        }

        lag_bounds = [ [-1000, 1000], [-500, 500] ]

        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                            run_pyccf=True, pyccf_params=pyccf_params,
                            lag_bounds=lag_bounds,
                            file_fmt='ascii')

        weighting_res, summary_dicts = pl.run_weighting(output_dir, line_names,
                                         run_pyccf=True, pyccf_params=pyccf_params,
                                         weighting_params=weighting_params,
                                         lag_bounds=lag_bounds,
                                         file_fmt='ascii')

        self.filenames = filenames
        self.line_names = line_names

        self.res1 = res
        self.res2 = weighting_res
        self.summary = summary_dicts


    def test_all(self):

        #################################################################################################
        # RES
        #################################################################################################
        #Make sure the lengths and keys of each of the resulting arrays are correct


        pyccf_keys = ['centroid', 'bounds', 'acf', 'lags', 'weight_dist',
                        'smoothed_dist', 'ntau', 'CCCD', 'downsampled_CCCD', 'frac_rejected']

        #Make sure the keys of the output are correct
        res_keys = list(self.res2['pyccf'].keys())
        for key in pyccf_keys:
            self.assertIn( key, res_keys )

        for mod in ['javelin', 'pyroa']:
            self.assertEqual( len(self.res2[mod].keys()), 0 )
            self.assertIn( 'rmax_'+mod, list(self.res2.keys()) )
            self.assertEqual( len( self.res2['rmax_'+mod] ), 0 )

        self.assertIn( 'rmax_pyccf', list(self.res2.keys()) )
        self.assertEqual( len( self.res2['rmax_pyccf'] ), len(self.filenames)-1 )

        array_mask = [True, True, True, True, True, True, True, True, True, False]

        #Make sure data types are correct
        for i in range(len(self.filenames)-1):

            for module, expected_keys in zip(['pyccf'], [pyccf_keys]):
                for j, key in enumerate(expected_keys):

                    if array_mask[j]:
                        self.assertIs( type(self.res2[module][key][i][0]), np.float64 )
                    else:
                        self.assertIs( type(self.res2[module][key][i]), float )


        #Make sure lengths of arrays are correct
        for key in pyccf_keys:
            self.assertEqual( len(self.res2['pyccf'][key]), len(self.filenames)-1 )

        for i in range(len(self.filenames)-1):
            self.assertEqual( len(self.res2['pyccf']['centroid'][i]), 3 )
            self.assertEqual( len(self.res2['pyccf']['bounds'][i]), 3 )


        #Make sure the reported lag is the med and error for the downsampled dist,
        #frac rejected is correct, ntau is correct

        xcont, _, _ = np.loadtxt(self.filenames[0], unpack=True, usecols=[0,1,2])
        for i in range(1, len(self.filenames)):
            xline, _, _ = np.loadtxt(self.filenames[i], unpack=True, usecols=[0,1,2])

            cccd = self.res1['pyccf_res'][i-1]['CCCD_lags']
            cccd_ds = self.res2['pyccf']['downsampled_CCCD'][i-1]


            #Lag vals
            self.assertEqual( self.res2['pyccf']['centroid'][i-1][1], np.median(cccd_ds) )
            self.assertEqual( self.res2['pyccf']['centroid'][i-1][0], np.median(cccd_ds) - np.percentile(cccd_ds, 16) )
            self.assertEqual( self.res2['pyccf']['centroid'][i-1][2], np.percentile(cccd_ds, 84) - np.median(cccd_ds) )

            #Frac rejected
            self.assertEqual( self.res2['pyccf']['frac_rejected'][i-1], 1- len(cccd_ds) / len(cccd) )

            #Ntau
            pyccf_lags = self.res2['pyccf']['lags'][i-1]

            _, pyccf_ntau, pyccf_n0, _ = prob_tau(xcont, xline, lagvals=pyccf_lags, gap_size=20, k=3)
            self.assertListEqual( list(self.res2['pyccf']['ntau'][i-1]), list(pyccf_ntau) )


        #################################################################################################
        # SUMMARY DICT
        #################################################################################################
        #Make sure the keys are correct and have correct length

        general_keys = ['n0', 'peak_bounds', 'peak', 'lag', 'lag_err', 'frac_rejected', 'rmax']
        isarray = [False, True, False, False, True, False, False]


        for i in range(len(self.filenames)-1):
            self.assertIn( 'k', list(self.summary[i].keys()) )

            for mod in ['pyccf', 'javelin', 'pyroa']:
                for key, arr in zip(general_keys, isarray):
                    self.assertIn( key + '_' + mod, list(self.summary[i].keys()) )

                    if arr:
                        self.assertTrue( isinstance( self.summary[i][key + '_' + mod], list ) )
                        self.assertEqual( len(self.summary[i][key+'_'+mod]), 2 )

                    if mod in ['javelin', 'pyroa']:
                        if arr:
                            self.assertTrue( np.all(np.isnan(self.summary[i][key+'_'+mod])) )
                        else:
                            self.assertTrue( np.isnan(self.summary[i][key+'_'+mod]) )


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

        #Make sure eachline has a "pyccf"and "weights" subdir
        for fdir, name in zip(subdirs[1:], self.line_names[1:]):
            for subdir in ['pyccf', 'weights']:
                self.assertIn( '.tmp/' + name + '/' + subdir, glob.glob( fdir + '*' ) )

        #Make sure continuum has no files
        self.assertEqual( len(glob.glob('.tmp/continuum/*')), 0 )

        #Make sure plots were output
        self.assertNotIn( '.tmp/javelin_weights_res.pdf', glob.glob('.tmp/*') )
        self.assertNotIn( '.tmp/pyroa_weights_res.pdf', glob.glob('.tmp/*') )
        self.assertIn( '.tmp/pyccf_weights_res.pdf', glob.glob('.tmp/*') )


        #Make sure each "weights" subdir has the proper files
        for i, fdir in zip([1,2], subdirs[1:]):
            files = glob.glob(fdir + 'weights/*')

            self.assertNotIn( fdir + 'weights/javelin_weighted_lag_dist.dat', files )
            self.assertNotIn( fdir + 'weights/javelin_weights.dat', files )
            self.assertNotIn( fdir + 'weights/pyroa_weighted_lag_dist.dat', files )
            self.assertNotIn( fdir + 'weights/pyroa_weights.dat', files )
            self.assertIn( fdir + 'weights/pyccf_weighted_cccd.dat', files )
            self.assertIn( fdir + 'weights/pyccf_weights.dat', files )
            self.assertIn( fdir + 'weights/weight_summary.fits', files )

            self.assertIn( fdir + 'weights/' + self.line_names[i] + '_weights.pdf', files )



        #################################################################################################
        # FILE RES MATCH
        #################################################################################################

        #Match output files
        for i, name in zip([1,2], self.line_names[1:]):
            cccd_ds = np.loadtxt('.tmp/' + name + '/weights/pyccf_weighted_cccd.dat')
            self.assertEqual( list(self.res2['pyccf']['downsampled_CCCD'][i-1]), list(cccd_ds) )

            lags, ntau, weight_dist, acf, smooth_dist, smooth_weight_dist = np.loadtxt( '.tmp/' + name + '/weights/pyccf_weights.dat', unpack=True, delimiter=',' )
            self.assertListEqual( list(self.res2['pyccf']['lags'][i-1]), list(lags) )
            self.assertListEqual( list(self.res2['pyccf']['ntau'][i-1]), list(ntau) )
            self.assertListEqual( list(self.res2['pyccf']['weight_dist'][i-1]), list(weight_dist) )
            self.assertListEqual( list(self.res2['pyccf']['acf'][i-1]), list(acf) )
            self.assertListEqual( list(self.res2['pyccf']['smoothed_dist'][i-1]), list(smooth_weight_dist) )



        #Match summary files
        xcont, _, _ = np.loadtxt(self.filenames[0], unpack=True, usecols=[0,1,2])
        for i, name in zip([1,2], self.line_names[1:]):
            read_dict = Table.read( '.tmp/' + name + '/weights/weight_summary.fits' )
            res_dict = self.summary[i-1]

            xline, _, _ = np.loadtxt(self.filenames[i], unpack=True, usecols=[0,1,2])

            #k
            self.assertEqual( read_dict['k'][0], 3 )
            self.assertEqual( res_dict['k'], 3 )

            #Get N0
            pyccf_lags = self.res2['pyccf']['lags'][i-1]
            _, _, pyccf_n0, _ = prob_tau(xcont, xline, lagvals=pyccf_lags, gap_size=20, k=3)

            self.assertEqual( read_dict['n0_pyccf'][0], pyccf_n0 )
            self.assertEqual( res_dict['n0_pyccf'], pyccf_n0 )


            #Bounds
            pyccf_bounds = [ self.res2['pyccf']['bounds'][i-1][0], self.res2['pyccf']['bounds'][i-1][2] ]
            pyccf_peak = self.res2['pyccf']['bounds'][i-1][1]


            self.assertTrue( np.allclose(list(read_dict['peak_bounds_pyccf'][0]), list(pyccf_bounds)) )
            self.assertTrue( np.allclose(list(res_dict['peak_bounds_pyccf']), list(pyccf_bounds)) )
            self.assertTrue( np.isclose(read_dict['peak_pyccf'][0], pyccf_peak) )
            self.assertEqual( res_dict['peak_pyccf'], pyccf_peak )


            #Lag
            pyccf_lag_err = [ self.res2['pyccf']['centroid'][i-1][0], self.res2['pyccf']['centroid'][i-1][2] ]
            pyccf_lag = self.res2['pyccf']['centroid'][i-1][1]

            self.assertTrue( np.allclose( list(read_dict['lag_err_pyccf'][0]), list(pyccf_lag_err)) )
            self.assertListEqual( list(res_dict['lag_err_pyccf']), list(pyccf_lag_err) )
            self.assertEqual( read_dict['lag_pyccf'], pyccf_lag )
            self.assertEqual( res_dict['lag_pyccf'], pyccf_lag )


            #Rmax
            self.assertTrue( np.isclose(self.res2['rmax_pyccf'][i-1], read_dict['rmax_pyccf'][0]) )
            self.assertTrue( np.isclose(read_dict['rmax_pyccf'][0], res_dict['rmax_pyccf']) )


            #Frac rejected
            cccd = self.res2['pyccf']['CCCD'][i-1]
            cccd_ds = self.res2['pyccf']['downsampled_CCCD'][i-1]

            self.assertEqual( read_dict['frac_rejected_pyccf'][0], 1- len(cccd_ds)/len(cccd) )
            self.assertEqual( res_dict['frac_rejected_pyccf'], 1- len(cccd_ds)/len(cccd) )


            for mod in ['javelin','pyroa']:
                self.assertTrue( isinstance(read_dict['n0_'+mod][0], np.ma.core.MaskedConstant) )
                self.assertTrue( isinstance(read_dict['peak_bounds_'+mod][0], np.ma.core.MaskedArray) )
                self.assertTrue( isinstance(read_dict['peak_'+mod][0], np.ma.core.MaskedConstant) )
                self.assertTrue( isinstance(read_dict['lag_err_'+mod][0], np.ma.core.MaskedArray) )
                self.assertTrue( isinstance(read_dict['lag_'+mod][0], np.ma.core.MaskedConstant) )
                self.assertTrue( isinstance(read_dict['frac_rejected_'+mod][0], np.ma.core.MaskedConstant) )
                self.assertTrue( isinstance(read_dict['rmax_'+mod][0], np.ma.core.MaskedConstant) )


    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
