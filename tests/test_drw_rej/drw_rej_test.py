import pypetal.pipeline as pl
import numpy as np 

import os
import glob
import shutil

import unittest


def str2bool(string_, default='raise'):

    true = ['true', 't', '1', 'y', 'yes', 'enabled', 'enable', 'on']
    false = ['false', 'f', '0', 'n', 'no', 'disabled', 'disable', 'off']
    if string_.lower() in true:
        return True
    elif string_.lower() in false or (not default):
        return False
    else:
        raise ValueError('The value \'{}\' cannot be mapped to boolean.'
                         .format(string_))



class TestDrwRej(unittest.TestCase):

    def setUp(self):

        main_dir = 'examples/dat/javelin_'            
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)


        params = {
            'nchain': 100,
            'nburn': 50,
            'nwalker': 20,
            'reject_data': [True, False, True],
            'nsig': 1.0
        }


        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                            run_drw_rej=True,
                            drw_rej_params=params, 
                            file_fmt='ascii')

        self.filenames = filenames
        self.line_names = line_names
        self.res = res



    #Make sure the lengths and keys of each of the resulting arrays are correct
    def test_res(self):
        
        expected_keys = ['masks', 'reject_data', 'taus', 'sigmas', 'jitters']
        expected_lengths = np.full( len(expected_keys), 3 )

        lc_lengths = []
        for i in range(len(self.filenames)):
            x, _, _ = np.loadtxt(self.filenames[i], unpack=True, usecols=[0,1,2])
            lc_lengths.append(len(x))

        mc_length = 100*20

        #Make sure the keys of the output are correct
        res_keys = list(self.res['drw_rej_res'].keys())
        for key in expected_keys:
            self.assertIn( key, res_keys )


        #Make sure lengths of key values are correct
        for key in expected_keys:
            self.assertEqual( len(self.res['drw_rej_res'][key]), expected_lengths[i] )


        #Make sure yelm mask is all false
        self.assertEqual( np.all(self.res['drw_rej_res']['masks'][1]), False )

        #Make sure yelm tau, sig, and jitter are None
        for key in['taus', 'sigmas', 'jitters']:
            self.assertIs( self.res['drw_rej_res'][key][1], None )

        #Make sure lengths of masks are correct
        for i in range(3):
            self.assertEqual( len(self.res['drw_rej_res']['masks'][i]), lc_lengths[i] )

        #Make sure MCMC lengths are correct
        for key in ['taus', 'sigmas', 'jitters']:
            for i in [0,2]:
                self.assertEqual( len(self.res['drw_rej_res'][key][i]), mc_length )




    #Make sure the layout of the files is correct
    def test_file_layout(self):
        main_directories = glob.glob('.tmp/*/')
        subdirs = ['.tmp/continuum/', '.tmp/zing/']

        #Make sure each line has a subdirectory
        self.assertIn( '.tmp/continuum/' , main_directories )
        self.assertIn( '.tmp/yelm/' , main_directories )
        self.assertIn( '.tmp/zing/' , main_directories )

        #Make sure yelm has no files
        self.assertEqual( len(glob.glob('.tmp/yelm/*')), 0 )

        #Make sure each line has a "drw_rej" subdirectory
        for fdir in subdirs:
            subdir_dirs = glob.glob(fdir + '*')
            self.assertIn( fdir + 'drw_rej' , subdir_dirs )

        #Make sure each "drw_rej" subdirectory has a chain,drw_fit, and mask file
        for i, fdir in zip([0,2], subdirs):
            files = glob.glob(fdir + 'drw_rej/*')
            self.assertIn( fdir + 'drw_rej/' + self.line_names[i] + '_chain.dat' , files )
            self.assertIn( fdir + 'drw_rej/' + self.line_names[i] + '_drw_fit.dat' , files )
            self.assertIn( fdir + 'drw_rej/' + self.line_names[i] + '_mask.dat' , files )
            self.assertIn( fdir + 'drw_rej/' + self.line_names[i] + '_drw_fit.pdf', files )

            #Make sure the mask file has the correct number of lines
            x, _, _ = np.loadtxt(self.filenames[i], unpack=True, usecols=[0,1,2])

            mask = np.loadtxt(fdir + 'drw_rej/' + self.line_names[i] + '_mask.dat', unpack=True, dtype=str)
            mask = list(map(str2bool, mask))

            self.assertEqual( len(mask), len(x) )



            #Make sure the chain file has the correct number of lines
            taus, sigs, jits = np.loadtxt(fdir + 'drw_rej/' + self.line_names[i] + '_chain.dat', unpack=True, delimiter=',', usecols=[0,1,2])            
            self.assertEqual( len(taus), 100*20 )


        #Make sure the light curves get saved
        self.assertIn( '.tmp/light_curves/', main_directories )

        lc_files = glob.glob('.tmp/light_curves/*')
        for name in self.line_names:
            self.assertIn( '.tmp/light_curves/' + name + '.dat', lc_files )

        #Make sure the masked light curves get saved
        self.assertIn( '.tmp/processed_lcs/', main_directories )

        lc_files = glob.glob('.tmp/processed_lcs/*')
        for i, name in enumerate(self.line_names):
            self.assertIn( '.tmp/processed_lcs/' + name + '_data.dat', lc_files )

            #Make sure the processed lcs have reasonable number of values
            x, _, _ = np.loadtxt('.tmp/processed_lcs/' + name + '_data.dat', unpack=True, delimiter=',', usecols=[0,1,2])
            xdata, _, _ = np.loadtxt(self.filenames[i], unpack=True, usecols=[0,1,2])
            self.assertLessEqual( len(x), len(xdata) )



    def test_res_file_match(self):

        #Match chains
        for i, name in zip([0,2], ['continuum', 'zing']):
            file_tau, file_sig, file_jit = np.loadtxt('.tmp/' + name + '/drw_rej/' + name + '_chain.dat', unpack=True, delimiter=',', usecols=[0,1,2])
            res_sig, res_tau, res_jit = self.res['drw_rej_res']['taus'][i], self.res['drw_rej_res']['sigmas'][i], self.res['drw_rej_res']['jitters'][i]

            self.assertListEqual(list(file_tau), list(res_tau))
            self.assertListEqual(list(file_sig), list(res_sig))
            self.assertListEqual(list(file_jit), list(res_jit))


        #Match masks
        for i, name in zip([0,2], ['continuum', 'zing']):
            file_mask = np.loadtxt('.tmp/' + name + '/drw_rej/' + name + '_mask.dat', unpack=True, dtype=str)
            file_mask = list(map(str2bool, file_mask))

            res_mask = self.res['drw_rej_res']['masks'][i]
            self.assertListEqual( list(file_mask), list(res_mask) )


    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')

