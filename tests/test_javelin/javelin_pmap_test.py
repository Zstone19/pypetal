import pypetal.pipeline as pl
from pypetal.load import get_javelin_together

import javelin.lcmodel

import os
import shutil

import unittest

class TestJAVELIN(unittest.TestCase):

    def setUp(self):

        main_dir = 'examples/dat/javelin_'
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']


        if os.path.exists(output_dir) and os.path.isdir(output_dir):
            shutil.rmtree(output_dir)


        params = {
            'nchain': 40,
            'nburn': 20,
            'nwalker': 20,
            'rm_type': "phot",
            'together': True
        }

        lag_bounds = [-500, 500]

        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                            run_javelin=True, javelin_params=params,
                            lag_bounds=lag_bounds,
                            file_fmt='ascii')

        self.filenames = filenames
        self.line_names = line_names
        self.res = res

        self.mc_length = params['nchain'] * params['nwalker']



    def test_pmap(self):

        #Make sure it sets together=False
        self.assertEqual( get_javelin_together('.tmp/'), False )

        #Make sure it returns a Pmap_Model object
        for i in range(len(self.line_names[1:])):
            self.assertIs( type(self.res['javelin_res'][i]['rmap_model']), javelin.lcmodel.Pmap_Model)

        #Make sure there are 4 tophat parameters
        for i in range(len(self.line_names[1:])):
            self.assertEqual( len(self.res['javelin_res'][i]['tophat_params']), 4)
            for j in range(4):
                self.assertEqual( len(self.res['javelin_res'][i]['tophat_params'][j]), self.mc_length)


    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
