import glob
import os
import shutil
import unittest

import javelin.lcmodel
import javelin.zylc
import numpy as np

import pypetal.pipeline as pl


class TestJavelinTogether(unittest.TestCase):

    def setUp(self):

        main_dir = 'examples/dat/javelin_'
        filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

        output_dir = '.tmp/'
        line_names = ['continuum', 'yelm', 'zing']


        if os.path.exists(output_dir) and os.path.isdir(output_dir):
            shutil.rmtree(output_dir)


        params = {
            'nchain': 30,
            'nburn': 20,
            'nwalker': 20,
            'together': True
        }
        lag_bounds = [ [-1000, 1000], [-500, 500] ]

        #Run pypetal
        res = pl.run_pipeline(output_dir, filenames, line_names,
                            run_javelin=True,
                            javelin_params=params,
                            lag_bounds=lag_bounds,
                            file_fmt='ascii')

        self.filenames = filenames
        self.line_names = line_names
        self.res = res




    def test_all(self):

        #################################################################################################
        # RES
        #################################################################################################
        #Make sure the lengths and keys of each of the resulting arrays are correct

        expected_keys = ['cont_hpd', 'tau', 'sigma', 'tophat_params', 'hpd',
                         'cont_model', 'rmap_model', 'cont_dat', 'tot_dat', 'bestfit_model']

        mc_length = 30*20

        #Make sure the keys of the output are correct
        res_keys = list(self.res['javelin_res'].keys())
        for key in expected_keys:
            self.assertIn( key, res_keys )


        #Make sure data types are correct
        for j, key in enumerate(expected_keys):
            if key == 'cont_model':
                self.assertIs( type(self.res['javelin_res'][key]), javelin.lcmodel.Cont_Model )
            elif key == 'rmap_model':
                self.assertIs( type(self.res['javelin_res'][key]), javelin.lcmodel.Rmap_Model )
            elif (key == 'cont_dat') or (key == 'tot_dat') or (key == 'bestfit_model'):
                self.assertIs( type(self.res['javelin_res'][key]), javelin.zylc.LightCurve )
            else:
                self.assertIs( self.res['javelin_res'][key].dtype.type, np.float64 )


        #Make sure lengths of arrays are correct
        self.assertEqual( len(self.res['javelin_res']['tau']), mc_length )
        self.assertEqual( len(self.res['javelin_res']['sigma']), mc_length )

        self.assertEqual( len(self.res['javelin_res']['tophat_params']), 6 )
        for i in range(6):
            self.assertEqual( len(self.res['javelin_res']['tophat_params'][i]), mc_length )

        self.assertEqual( len(self.res['javelin_res']['cont_hpd']), 3 )
        self.assertEqual( len(self.res['javelin_res']['cont_hpd'][0]), 2 )

        self.assertEqual( len(self.res['javelin_res']['hpd']), 3 )
        self.assertEqual( len(self.res['javelin_res']['hpd'][0]), 8 )






        #Make sure that the LC objetcs have the right light curves
        x_cont, y_cont, yerr_cont = np.loadtxt( self.filenames[0], unpack=True, usecols=[0,1,2] )

        #cont_dat
        jav_x = self.res['javelin_res']['cont_dat'].jlist[0]
        jav_y = self.res['javelin_res']['cont_dat'].mlist[0] + self.res['javelin_res'][ 'cont_dat'].blist[0]
        jav_yerr = self.res['javelin_res']['cont_dat'].elist[0]

        self.assertListEqual( list(x_cont), list(jav_x) )
        self.assertListEqual( list(yerr_cont), list(jav_yerr) )
        for j in range(len(y_cont)):
            self.assertAlmostEqual( y_cont[j], jav_y[j], places=6 )


        #Tot_dat
        jav_x = self.res['javelin_res']['tot_dat'].jlist[0]
        jav_y = self.res['javelin_res']['tot_dat'].mlist[0] + self.res['javelin_res'][ 'tot_dat'].blist[0]
        jav_yerr = self.res['javelin_res']['tot_dat'].elist[0]

        self.assertListEqual( list(x_cont), list(jav_x) )
        self.assertListEqual( list(yerr_cont), list(jav_yerr) )
        for j in range(len(y_cont)):
            self.assertAlmostEqual( y_cont[j], jav_y[j], places=6 )

        for i in range(1, len(self.filenames)):
            x, y, yerr = np.loadtxt(self.filenames[i], unpack=True, usecols=[0,1,2])

            jav_x = self.res['javelin_res']['tot_dat'].jlist[i]
            jav_y = self.res['javelin_res']['tot_dat'].mlist[i] + + self.res['javelin_res'][ 'tot_dat'].blist[i]
            jav_yerr = self.res['javelin_res']['tot_dat'].elist[i]

            self.assertListEqual( list(x), list(jav_x) )
            self.assertListEqual( list(yerr), list(jav_yerr) )

            for j in range(len(y)):
                self.assertAlmostEqual( y[j], jav_y[j], places=6 )







        #Make sure that the LCModel objects have the right light curves

        #Cont_Model
        jav_x = self.res['javelin_res'][ 'cont_model'].zydata.jlist[0]
        jav_y = self.res['javelin_res'][ 'cont_model'].zydata.mlist[0] + self.res['javelin_res'][ 'cont_model'].zydata.blist[0]
        jav_yerr = self.res['javelin_res'][ 'cont_model'].zydata.elist[0]

        self.assertListEqual( list(x_cont), list(jav_x) )
        self.assertListEqual( list(yerr_cont), list(jav_yerr) )
        for j in range(len(y_cont)):
            self.assertAlmostEqual( y_cont[j], jav_y[j], places=6 )


        #Rmap_Model
        jav_x = self.res['javelin_res'][ 'rmap_model'].zydata.jlist[0]
        jav_y = self.res['javelin_res'][ 'rmap_model'].zydata.mlist[0] + self.res['javelin_res'][ 'rmap_model'].zydata.blist[0]
        jav_yerr = self.res['javelin_res'][ 'rmap_model'].zydata.elist[0]

        self.assertListEqual( list(x_cont), list(jav_x) )
        self.assertListEqual( list(yerr_cont), list(jav_yerr) )
        for j in range(len(y_cont)):
            self.assertAlmostEqual( y_cont[j], jav_y[j], places=6 )



        for i in range(1, len(self.filenames)):
            x, y, yerr = np.loadtxt(self.filenames[i], unpack=True, usecols=[0,1,2])

            jav_x = self.res['javelin_res'][ 'rmap_model'].zydata.jlist[i]
            jav_y = self.res['javelin_res'][ 'rmap_model'].zydata.mlist[i] + self.res['javelin_res'][ 'rmap_model'].zydata.blist[i]
            jav_yerr = self.res['javelin_res'][ 'rmap_model'].zydata.elist[i]

            self.assertListEqual( list(x), list(jav_x) )
            self.assertListEqual( list(yerr), list(jav_yerr) )
            for j in range(len(y)):
                self.assertAlmostEqual( y[j], jav_y[j], places=6 )



        lag_bounds = [ [-1000, 1000], [-500, 500] ]
        #Make sure lag bounds worked
        for i in range(len(self.filenames)-1):
            self.assertGreaterEqual( np.min(self.res['javelin_res']['tophat_params'][3*i]), lag_bounds[i][0] )
            self.assertLessEqual( np.max(self.res['javelin_res']['tophat_params'][3*i]), lag_bounds[i][1] )



        #################################################################################################
        # FILE LAYOUT
        #################################################################################################
         #Make sure the layout of the files is correct

        main_directories = glob.glob('.tmp/*/')
        subdirs = ['.tmp/continuum/', '.tmp/yelm/', '.tmp/zing/']

        mc_length = 30*20
        burn_length = 20*20

        #Make sure each line has a subdirectory
        self.assertIn( '.tmp/continuum/' , main_directories )
        self.assertIn( '.tmp/yelm/' , main_directories )
        self.assertIn( '.tmp/zing/' , main_directories )
        self.assertIn( '.tmp/javelin/', main_directories )

        #Make sure line subdirs has no files
        self.assertEqual( len(glob.glob('.tmp/continuum/*')), 0 )
        self.assertEqual( len(glob.glob('.tmp/yelm/*')), 0 )
        self.assertEqual( len(glob.glob('.tmp/zing/*')), 0 )


        #Make sure "javelin" subdirectory has chain and burn files
        files = glob.glob('.tmp/javelin/*')
        for x in ['cont', 'rmap']:
            self.assertIn( '.tmp/javelin/burn_' + x + '.txt', files )
            self.assertIn( '.tmp/javelin/chain_' + x + '.txt', files )
            self.assertIn( '.tmp/javelin/logp_' + x + '.txt', files )

        for name in self.line_names:
            self.assertIn( '.tmp/javelin/' + name + '_lc_fits.dat', files )

        self.assertIn( '.tmp/javelin/cont_lcfile.dat', files )
        self.assertIn( '.tmp/javelin/tot_lcfile.dat', files )

        self.assertIn( '.tmp/javelin/javelin_bestfit.pdf', files )
        self.assertIn( '.tmp/javelin/javelin_corner.pdf', files )
        self.assertIn( '.tmp/javelin/javelin_histogram.pdf', files )




        #Make sure the lc files are correct
        xcont, ycont, yerrcont = np.loadtxt( self.filenames[0], unpack=True, usecols=[0,1,2] )

        cont_dat = javelin.zylc.get_data( '.tmp/javelin/cont_lcfile.dat' )
        jav_x = cont_dat.jlist[0]
        jav_y = cont_dat.mlist[0] + cont_dat.blist[0]
        jav_yerr = cont_dat.elist[0]

        self.assertListEqual( list(xcont), list(cont_dat.jlist[0]) )
        self.assertListEqual( list(ycont), list(cont_dat.mlist[0] + cont_dat.blist[0]) )
        for j in range(len(ycont)):
            self.assertAlmostEqual( ycont[j], jav_y[j], places=6 )




        tot_dat = javelin.zylc.get_data( '.tmp/javelin/tot_lcfile.dat' )
        jav_x = tot_dat.jlist[0]
        jav_y = tot_dat.mlist[0] + tot_dat.blist[0]
        jav_yerr = tot_dat.elist[0]

        self.assertListEqual( list(xcont), list(tot_dat.jlist[0]) )
        self.assertListEqual( list(yerrcont), list(tot_dat.elist[0]) )
        for j in range(len(jav_x)):
            self.assertAlmostEqual( ycont[j], jav_y[j], places=6 )



        for i in range(1, len(self.filenames)):
            xline, yline, yerrline = np.loadtxt( self.filenames[i], unpack=True, usecols=[0,1,2] )

            jav_x = tot_dat.jlist[i]
            jav_y = tot_dat.mlist[i] + tot_dat.blist[i]
            jav_yerr = tot_dat.elist[i]

            self.assertListEqual( list(xline), list(tot_dat.jlist[i]) )
            self.assertListEqual( list(yerrline), list(tot_dat.elist[i]) )
            for j in range(len(jav_x)):
                self.assertAlmostEqual( yline[j], jav_y[j], places=6 )



        #Make sure chain and burn files have correct length
        for x in ['cont', 'rmap']:
            for prefix in ['burn', 'chain', 'logp']:
                chain = np.loadtxt( '.tmp/javelin/' + prefix + '_' + x + '.txt' )
                self.assertEqual( len(chain), burn_length if prefix == 'burn' else mc_length )


        #Make sure the light curves get saved
        self.assertIn( '.tmp/light_curves/', main_directories )

        lc_files = glob.glob('.tmp/light_curves/*')
        for name in self.line_names:
            self.assertIn( '.tmp/light_curves/' + name + '.dat', lc_files )




        #################################################################################################
        # RES FILE MATCH
        #################################################################################################

        #Match dists
        file_sig, file_tau, t1, w1, s1, t2, w2, s2 = np.loadtxt( '.tmp/javelin/chain_rmap.txt', unpack=True )
        file_sig = np.exp(file_sig)
        file_tau = np.exp(file_tau)
        file_tophat = [t1, w1, s1, t2, w2, s2]

        self.assertListEqual( list(file_sig), list(self.res['javelin_res']['sigma']) )
        self.assertListEqual( list(file_tau), list(self.res['javelin_res']['tau']) )
        for j in range(6):
            self.assertListEqual( list(file_tophat[j]), list(self.res['javelin_res']['tophat_params'][j]) )


    def tearDown(self):
        if os.path.exists('.tmp/') and os.path.isdir('.tmp/'):
            shutil.rmtree('.tmp/')
