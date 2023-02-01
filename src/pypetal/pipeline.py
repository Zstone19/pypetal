import pypetal.modules as modules
from pypetal import weighting, defaults
from pypetal.petalio import write_data, make_directories
import pypetal.detrending as dtr

import os

import numpy as np
from astropy.table import Table

    
def run_pipeline(output_dir, arg2, 
                 line_names=None, 
                 run_drw_rej=False, drw_rej_params={},
                 run_detrend=False, detrend_params={},
                 run_pyccf=False, pyccf_params={},
                 run_pyzdcf=False, pyzdcf_params={},
                 run_javelin=False, javelin_params={},
                 run_weighting=False, weighting_params={},
                 **kwargs):
    
    """ Run the pyPetal pipeline on a list of files. Individual modules can be specified, but are not run by default.
    
    Parameters
    ----------
    
    output_dir : str
        The directory to save the output files to.

    arg2 : list of str, list of floats
        If a list of strings, the filenames of the light curves to run the pipeline on. If a list of floats, the light curves themselves.

    line_names : list of str, optional
        The names of the lines to be used in the output files. If ``None``, the lines will be named as ``line1``, ``line2``,... Default is ``None``.

    run_drw_rej : bool, optional
        Whether to run the DRW rejection module. Default is ``False``.

    drw_rej_params : dict, optional
        The parameters to pass to the DRW rejection module. Default is ``{}``.

    run_detrend : bool, optional
        Whether to run the detrending module. Default is ``False``.

    detrend_params : dict, optional
        The parameters to pass to the detrending module. Default is ``{}``.

    run_pyccf : bool, optional
        Whether to run the pyCCF module. Default is ``False``.

    pyccf_params : dict, optional
        The parameters to pass to the pyCCF module. Default is ``{}``.

    run_pyzdcf : bool, optional
        Whether to run the pyZDCF module. Default is ``False``.

    pyzdcf_params : dict, optional
        The parameters to pass to the pyZDCF module. Default is ``{}``.

    run_javelin : bool, optional
        Whether to run the Javelin module. Default is ``False``.

    javelin_params : dict, optional
        The parameters to pass to the Javelin module. Default is ``{}``.

    run_weighting : bool, optional
        Whether to run the weighting module. Default is ``False``.

    weighting_params : dict, optional
        The parameters to pass to the weighting module. Default is ``{}``.




    Returns
    -------

    output : dict
        A dictionary containing the output of the pipeline modules.

    """

    output_dir = os.path.abspath(output_dir) + r'/'

    if arg2 is None:
        raise Exception('Please provide a list of light curve filenames or the light curves themselves')

    
    if not isinstance(arg2[0], str):
        os.makedirs( output_dir + 'input_lcs/', exist_ok=True )        
        fnames = []
        
        for i in range( len(arg2) ):
            
            if i == 0:
                name = 'continuum' 
            else:
                name = 'line{}'.format(i+1)
            
            write_data( arg2[i], output_dir + 'input_lcs/' + name + '.dat' )
            fnames.append( output_dir + 'input_lcs/' + name + '.dat' )
            
        fnames = np.array(fnames)
        kwargs['file_fmt'] = 'csv'
    else:
        fnames = arg2
        
    


    if len(fnames) < 2:
        print('ERROR: Requires at least two light curves to run pipeline.')
        return {}
    
    if len(line_names) != len(fnames):
        print('ERROR: Must have the same number of line names as light curves.')
        return {}
    
    
    cont_fname = fnames[0]
    line_fnames = fnames[1:]
    
    if type(line_fnames) is str:
        line_fnames = [line_fnames]
    
    
    if (run_weighting) & (not run_pyccf) & (not run_javelin):
        print('ERROR: Either JAVELIN or pyCCF must be run before weighting can be done.')
        print('Setting run_weighting=False')
        run_weighting = False

    
    #Read in general kwargs
    general_kwargs = defaults.set_general(kwargs, fnames)

    #Get "reject_data" and "together"    
    _, _, _, _, _, _, reject_data, use_for_javelin = defaults.set_drw_rej(drw_rej_params, fnames)
    _, fixed, p_fix, _, _, _, _, _, _, _, _, _, _, together, _ = defaults.set_javelin(javelin_params, fnames)   
    
    javelin_params['fixed'] = fixed
    javelin_params['p_fix'] = p_fix

    
    
    #Name lines if unnamed
    if line_names is None:
        line_names = np.zeros( len(line_fnames), dtype=str )
        for i in range(len(line_fnames)):
            line_names.append('Line {}'.format(i+1))


    #Get directory common to all files
    input_dir = os.path.dirname( os.path.realpath(cont_fname) )
    for i in range(len(line_fnames)):
        assert input_dir == os.path.dirname( os.path.realpath(line_fnames[i]) )        
    
    
    
    #Create subdirectories for each line and module
    make_directories(output_dir, fnames, line_names, 
                     run_drw_rej, run_pyccf, run_pyzdcf, 
                     run_javelin, run_weighting,
                     reject_data, together)
        
    if general_kwargs['file_fmt'] != 'csv':
        #Make temp directory to store csv files
        os.makedirs( output_dir + 'temp/', exist_ok=True )
        input_dir = output_dir + 'temp/'

        #Read in file        
        try:
            dat = Table.read( cont_fname, format=file_fmt)
        except:
            dat = Table.read( cont_fname, format='ascii')
        
    
        #Get filename
        fname = os.path.basename( cont_fname )
        
        #Write to csv
        write_data( [ dat[ dat.colnames[0] ], dat[ dat.colnames[1] ], dat[ dat.colnames[2] ] ], input_dir + fname )        
        cont_fname = input_dir + fname        
        
        
        for i in range(len(line_fnames)):
            
            #Read in file
            try:
                dat = Table.read( line_fnames[i], format=file_fmt)
            except:
                dat = Table.read( line_fnames[i], format='ascii')
                
            #Get filename
            fname = os.path.basename( line_fnames[i] )
                
            #Write to csv
            write_data( [ dat[ dat.colnames[0] ], dat[ dat.colnames[1] ], dat[ dat.colnames[2] ] ], input_dir + fname )
            line_fnames[i] = input_dir + fname

        
    drw_rej_res = {}
    detrend_res = {}
    pyccf_res = {}
    pyzdcf_res = {}
    javelin_res = {}

    if run_drw_rej:
        
        if np.any( reject_data ):
            os.makedirs( output_dir + 'processed_lcs/', exist_ok=True )
        
        drw_rej_res = modules.drw_rej_tot( cont_fname, line_fnames, line_names, output_dir, general_kwargs, drw_rej_params ) 
                
        #Output light curve files with masks
        for i in range(len(fnames)):
            
            if i == 0:
                x, y, yerr = np.loadtxt( cont_fname, delimiter=',', unpack=True, usecols=[0,1,2] )
            else:
                x, y, yerr = np.loadtxt( line_fnames[i-1], delimiter=',', unpack=True, usecols=[0,1,2] )
            
            mask = drw_rej_res['masks'][i]
            write_data( [x,y,yerr,mask], output_dir + 'light_curves/' + line_names[i] + '.dat', '#x,y,yerr,mask' )
                            
                
        #If rejecting any data, make the new files the ones without rejected data
        if np.any( reject_data ):
            cont_fname = output_dir + 'processed_lcs/' + line_names[0] + '_data.dat'
            
            for i in range(len(line_fnames)):
                line_fnames[i] = output_dir + 'processed_lcs/' + line_names[i+1] + '_data.dat'            
            
            
        if use_for_javelin:            
            tau_cont = np.median(drw_rej_res['taus'][0])
            sig_cont = np.median(drw_rej_res['sigmas'][0])
            
            if together:                
            
                if fixed is not None:
                    javelin_params['fixed'][0] = 0
                    javelin_params['fixed'][1] = 0
                    
                    javelin_params['p_fix'][0] = np.log(sig_cont)
                    javelin_params['p_fix'][1] = np.log(tau_cont)
                else:                        
                    javelin_params['fixed'] = np.ones( 2 + 3*len(line_fnames) )
                    javelin_params['p_fix'] = np.zeros( 2 + 3*len(line_fnames) )
                
                    javelin_params['fixed'][0] = 0
                    javelin_params['fixed'][1] = 0
                    
                    javelin_params['p_fix'][0] = np.log(sig_cont)
                    javelin_params['p_fix'][1] = np.log(tau_cont)                                  
                    
            else:
                for i in range(len(line_fnames)):
                    if fixed[i] is not None:
                        javelin_params['fixed'][i][0] = 0
                        javelin_params['fixed'][i][1] = 0
                        
                        javelin_params['p_fix'][i][0] = np.log(sig_cont)
                        javelin_params['p_fix'][i][1] = np.log(tau_cont)
                        
                    else:
                        javelin_params['fixed'][i] = np.ones(5)
                        javelin_params['p_fix'][i] = np.zeros(5)

                        javelin_params['fixed'][i][0] = 0
                        javelin_params['fixed'][i][1] = 0
                        
                        javelin_params['p_fix'][i][0] = np.log(sig_cont)
                        javelin_params['p_fix'][i][1] = np.log(tau_cont)
                            
    else:
        
        #Output light curve files with masks
        for i in range(len(fnames)):
            
            if i == 0:
                x, y, yerr = np.loadtxt( cont_fname, delimiter=',', unpack=True, usecols=[0,1,2] )
            else:
                x, y, yerr = np.loadtxt( line_fnames[i-1], delimiter=',', unpack=True, usecols=[0,1,2] )
            
            write_data( [x,y,yerr,np.full_like(x, False, dtype=bool)], output_dir + 'light_curves/' + line_names[i] + '.dat' )
            
            
    if run_detrend:
        
        if not run_drw_rej:
            os.makedirs( output_dir + 'processed_lcs/', exist_ok=True )
        elif (run_drw_rej) & ( ~np.all(reject_data) ):
            os.makedirs( output_dir + 'processed_lcs/', exist_ok=True )
        
        detrend_res = dtr.detrend_tot(output_dir, cont_fname, line_fnames, line_names, general_kwargs, detrend_params)
        
        if not run_drw_rej:
            cont_fname = output_dir + 'processed_lcs/' + line_names[0] + '_data.dat'
            
            for i in range(len(line_fnames)):
                line_fnames[i] = output_dir + 'processed_lcs/' + line_names[i+1] + '_data.dat'                    

    if run_pyccf:
        pyccf_res = modules.pyccf_tot(cont_fname, line_fnames, line_names, output_dir, general_kwargs, pyccf_params)        
        
    if run_pyzdcf:
        pyzdcf_res, plike_res = modules.pyzdcf_tot(cont_fname, line_fnames, line_names, output_dir, general_kwargs, pyzdcf_params)
        
    if run_javelin:
        javelin_res = modules.javelin_tot(cont_fname, line_fnames, line_names, output_dir, general_kwargs, javelin_params)
    
    if run_weighting:
        weighting_res = weighting.run_weighting(cont_fname, line_fnames, output_dir, line_names, 
                                run_pyccf, run_javelin,
                                pyccf_res, javelin_res,
                                pyccf_params, javelin_params,
                                general_kwargs, weighting_params)


    if general_kwargs['file_fmt'] != 'csv':
        import shutil
        shutil.rmtree(input_dir)


    #Compile all results into a single dict    
    tot_res = {}
    
    if run_drw_rej:
        tot_res['drw_rej_res'] = drw_rej_res
    
    if run_detrend:
        tot_res['detrend_res'] = detrend_res
    
    if run_javelin:
        tot_res['javelin_res'] = javelin_res
    
    if run_pyccf:
        tot_res['pyccf_res'] = pyccf_res
        
    if run_pyzdcf:
        tot_res['pyzdcf_res'] = pyzdcf_res
        tot_res['plike_res'] = plike_res
        
    if run_weighting:
        tot_res['weighting_res'] = weighting_res



    if not isinstance(arg2[0], str):
        import shutil
        shutil.rmtree( output_dir + 'input_lcs/' ) 

    
    return tot_res