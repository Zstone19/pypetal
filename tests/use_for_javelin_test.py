import pypetal.pipeline as pl
import numpy as np 

import os
import glob
import shutil

import unittest


class TestDrwRej(unittest.TestCase):

    def setUp(self):

        main_dir = 'test_repo/pypetal/examples/dat/javelin_'
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)


        drw_rej_params = {
            'nchain': 100,
            'nburn': 50,
            'nwalker': 20,
            'reject_data': [True, False, True],
            'nsig': 1.0,
            'use_for_javelin': True
        }

        javelin_params = {
            'nchain': 20,
            'nburn': 10,
            'nwalker': 10
        }

        lag_bounds = [ [-1000, 1000], [-500, 500] ]

        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                            run_drw_rej=True, drw_rej_params=drw_rej_params, 
                            run_javelin=True, javelin_params=javelin_params,
                            lag_bounds=lag_bounds,
                            file_fmt='ascii')

        self.filenames = filenames
        self.line_names = line_names
        self.res = res


    #Make chain is only one value for sigma and tau
    def test_chain(self):
        drw_sigma = np.median(self.res['drw_rej_res']['sigmas'][0])
        drw_tau = np.median(self.res['drw_rej_res']['taus'][0])

        sigma_arr = np.full( 20*10, drw_sigma )
        tau_arr = np.full( 20*10, drw_tau )

        for i in range(len(self.filenames)-1):
            for j in range(20*10):
                self.assertAlmostEqual( self.res['javelin_res'][i]['sigma'][j], drw_sigma, places=6 )
                self.assertAlmostEqual( self.res['javelin_res'][i]['tau'][j], drw_tau, places=6 )
