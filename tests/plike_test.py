import glob
import os
import shutil
import unittest

import numpy as np
from astropy.table import Table

import pypetal.pipeline as pl


def format_float(x):
    return float("%11.3e" % x)

def format_int(x):
    return int("%5.1i" % x)


class TestPLIKE(unittest.TestCase):

    def setUp(self):

        main_dir = 'examples/dat/javelin_'
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']


        if os.path.exists(output_dir) and os.path.isdir(output_dir):
            shutil.rmtree(output_dir)


        params = {
            'nsim': 100,
            'minpts': 12,
            'prefix': 'myname',
            'run_plike': True,
            'plike_dir': 'plike_v4/'
        }

        lag_bounds = [-500, 500]

        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                            run_pyzdcf=True,
                            pyzdcf_params=params,
                            lag_bounds=lag_bounds,
                            file_fmt='ascii')

        self.lag_bounds = lag_bounds
        self.filenames = filenames
        self.line_names = line_names
        self.res = res


    def test_plike(self):

        #################################################################################################
        # RES
        #################################################################################################
        #Make sure the lengths and keys of each of the resulting arrays are correct

        df_cols = [ 'num', 'lag', 'r', '-dr', '+dr', 'likelihood' ]
        plike_keys = ['output', 'ML_lag', 'ML_lag_err_lo', 'ML_lag_err_hi']

        num_df = len(self.line_names) - 1
        self.assertEqual( num_df, len(self.res['plike_res']) )
        for i in range(num_df):
            for key in plike_keys:
                self.assertIn( key, self.res['plike_res'][i].keys() )

        #Make sure colnames are correct
        for i in range(num_df):
            self.assertListEqual( df_cols, list( self.res['plike_res'][i]['output'].colnames ) )


        #Make sure data types are correct
        for i in range(num_df):

            for col in df_cols:
                if col == 'num':
                    self.assertIs( type(self.res['plike_res'][i]['output'][col][0]), np.int64 )
                else:
                    self.assertIs( type(self.res['plike_res'][i]['output'][col][0]), np.float64 )

            for key in plike_keys[1:]:
                self.assertIn( type(self.res['plike_res'][i][key]), [float, np.float64] )


        #Make sure lags have proper values
        for i in range(num_df):
            self.assertGreaterEqual( np.min(self.res['plike_res'][i]['output']['lag']), self.lag_bounds[0] )
            self.assertLessEqual( np.max(self.res['plike_res'][i]['output']['lag']), self.lag_bounds[1] )



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

        #Make sure each line has a "pyzdcf" subdirectory
        for fdir in subdirs[1:]:
            subdir_dirs = glob.glob(fdir + '*')
            self.assertIn( fdir + 'pyzdcf' , subdir_dirs )


        #Make sure each "pyzdcf" subdirectory has a plike file
        for i, fdir in zip([1,2], subdirs[1:]):
            files = glob.glob(fdir + 'pyzdcf/*')
            self.assertIn( fdir + 'pyzdcf/' + self.line_names[i] + '_plike.out' , files )


        #Make sure the light curves get saved
        self.assertIn( '.tmp/light_curves/', main_directories )

        lc_files = glob.glob('.tmp/light_curves/*')
        for name in self.line_names:
            self.assertIn( '.tmp/light_curves/' + name + '.dat', lc_files )




        #################################################################################################
        # RES FILE MATCH
        #################################################################################################

        df_cols = [ 'num', 'lag', 'r', '-dr', '+dr', 'likelihood' ]

        for i in range(len(self.line_names) - 1):
            file_df = Table.read('.tmp/' + self.line_names[i+1] + '/pyzdcf/' + self.line_names[i+1] + '_plike.out',
                                 format='ascii',
                                 names=['num', 'lag', 'r', '-dr', '+dr', 'likelihood'])

            res_df = self.res['plike_res'][i]['output']

            # for col in df_cols:
            #     if col != '#bin':
            #         res_df[col] = list(map( format_float, res_df[col] ))
            #     else:
            #         res_df[col] = list(map( format_int, res_df[col] ))

            for col in df_cols:
                self.assertListEqual( list(res_df[col]), list(file_df[col]) )


            with open('.tmp/' + self.line_names[i+1] + '/pyzdcf/' + self.line_names[i+1] + '_plike.out', 'r') as f:
                output_str = list(f)[-3:]

                ml_lag = float( output_str[1].split()[7] )
                ml_lag_err_hi = np.abs( float( output_str[1].split()[8] )  )
                ml_lag_err_lo = np.abs( float( output_str[1].split()[9] )  )


            self.assertEqual( ml_lag, self.res['plike_res'][i]['ML_lag'] )
            self.assertEqual( ml_lag_err_hi, self.res['plike_res'][i]['ML_lag_err_hi'] )
            self.assertEqual( ml_lag_err_lo, self.res['plike_res'][i]['ML_lag_err_lo'] )



    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
