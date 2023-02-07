import os
import shutil
import unittest

import numpy as np
from scipy.interpolate import interp1d

import pypetal.pipeline as pl


class TestLoad(unittest.TestCase):

    def setUp(self):
        main_dir = 'examples/dat/javelin_'
        self.filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        self.output_dir = '.tmp/'
        self.line_names = ['continuum', 'yelm', 'zing']


    def test_drw_rejection(self):

        ###############################################################################################
        # RUN PYPETAL
        ###############################################################################################

        drw_rej_params = {
            'nchain': 100,
            'nburn': 50,
            'nwalker': 20,
            'reject_data': True,
            'nsig': 2.5
        }

        #Run pypetal
        res = pl.run_pipeline(self.output_dir, self.filenames, self.line_names,
                            run_drw_rej=True, drw_rej_params=drw_rej_params,
                            file_fmt='ascii',
                            time_unit='d',
                            lc_unit=['Jy', 'mJy', 'mJy'])

        ###############################################################################################
        # TEST REJECTION
        ###############################################################################################

        for i, name in enumerate(self.line_names):
            xfit, yfit, yerrfit = np.loadtxt( '.tmp/' + name + '/drw_rej/' + name + '_drw_fit.dat', unpack=True, delimiter=',' )
            xdat, ydat, yerrdat = np.loadtxt( '.tmp/light_curves/' + name + '.dat' , unpack=True, delimiter=',', usecols=[0,1,2] )
            res_mask = res['drw_rej_res']['masks'][i]

            interp_func = interp1d(xfit, yfit, kind='linear')

            for j in range(len(xdat)):
                fit_val = interp_func(xdat[j])
                self.assertEqual( res_mask[j], np.abs(ydat[j]-yfit[fit_ind]) > 2.5*yerrdat[j] )


    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
