import pypetal.pipeline as pl
import numpy as np 
import pandas as pd

import glob
import os
import shutil

import unittest


def format_float(x):
    return float("%11.3e" % x)

def format_int(x):
    return int("%5.1i" % x)


class TestPyZDCF(unittest.TestCase):

    def setUp(self):

        main_dir = '../examples/dat/javelin_'
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']


        if os.path.exists(output_dir) and os.path.isdir(output_dir):
            shutil.rmtree(output_dir)


        params = {
            'nsim': 300,
            'minpts': 12,
            'prefix': 'myname'
        }

        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                            run_pyzdcf=True,
                            pyzdcf_params=params, 
                            file_fmt='ascii')

        self.filenames = filenames
        self.line_names = line_names
        self.res = res

    #Make sure the lengths and keys of each of the resulting arrays are correct
    def test_res(self):
        
        df_cols = ['tau', '-sig(tau)', '+sig(tau)', 'dcf', '-err(dcf)', '+err(dcf)', '#bin' ]

        num_df = len(self.line_names) - 1
        self.assertEqual( num_df, len(self.res['pyzdcf_res']) )

        #Make sure colnames are correct
        for i in range(num_df):
            self.assertListEqual( df_cols, list( self.res['pyzdcf_res'][i].columns ) )


        #Make sure data types are correct
        for i in range(num_df):
            for col in df_cols:
                self.assertIs( type(self.res['pyzdcf_res'][i][col].iloc[0]), np.float64 )

        #Make sure ZDCF has proper values
        for i in range(num_df):
            for j in range(len(self.res['pyzdcf_res'][i]['dcf'])):
                self.assertGreaterEqual( np.min(self.res['pyzdcf_res'][i]['dcf'].iloc[j]), -1 )
                self.assertLessEqual( np.max(self.res['pyzdcf_res'][i]['dcf'].iloc[j]), 1 )




    #Make sure the layout of the files is correct
    def test_file_layout(self):
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


        #Make sure each "pyzdcf" subdirectory has a dcf file
        for i, fdir in zip([1,2], subdirs[1:]):
            files = glob.glob(fdir + 'pyzdcf/*')
            self.assertIn( fdir + 'pyzdcf/' + self.line_names[i] + '_myname.dcf' , files )
            self.assertIn( fdir + 'pyzdcf/' + self.line_names[i] + '_zdcf.pdf' , files )


        #Make sure the light curves get saved
        self.assertIn( '.tmp/light_curves/', main_directories )

        lc_files = glob.glob('.tmp/light_curves/*')
        for name in self.line_names:
            self.assertIn( '.tmp/light_curves/' + name + '.dat', lc_files )




    def test_res_file_match(self):

        df_cols = ['tau', '-sig(tau)', '+sig(tau)', 'dcf', '-err(dcf)', '+err(dcf)', '#bin' ]

        for i in range(len(self.line_names) - 1):
            file_df = np.loadtxt('.tmp/' + self.line_names[i+1] + '/pyzdcf/' + self.line_names[i+1] + '_myname.dcf')
            file_df = pd.DataFrame(file_df, columns=df_cols)
            
            res_df = self.res['pyzdcf_res'][i].copy(deep=True)
            
            for col in df_cols:
                if col != '#bin':                
                    res_df[col] = list(map( format_float, res_df[col] ))
                else:
                    res_df[col] = list(map( format_int, res_df[col] ))


            for col in df_cols:
                self.assertListEqual( list(res_df[col]), list(file_df[col]) )
