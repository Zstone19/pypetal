import glob
import os

import numpy as np
from astropy.table import Table

from pypetal.utils.petalio import make_directories, write_data, print_header, print_error, print_warning

try:
    import pypetal.utils.detrending as dtr
    detrend_exist = True
except ModuleNotFoundError:
    print_warning('WARNING: Could not successfully import the detrending module. Assuming run_detrend=False.')
    detrend_exist = False

try:
    from pypetal.mica2.module import mica2_tot
    mica2_exist = True
except ModuleNotFoundError:
    print_warning('WARNING: Could not successfully import the MICA2 module. Assuming run_mica2=False.')
    mica2_exist = False
    
    
if mica2_exist:
    from mpi4py import MPI



from pypetal.drw_rej.module import drw_rej_tot
from pypetal.pyccf.module import pyccf_tot
from pypetal.pyroa.module import pyroa_tot
from pypetal.pyzdcf.module import pyzdcf_tot
from pypetal.utils import defaults
from pypetal.weighting.module import run_weighting_tot


def run_pipeline(output_dir, arg2,
                 line_names=None,
                 run_drw_rej=False, drw_rej_params={},
                 run_detrend=False, detrend_params={},
                 run_pyccf=False, pyccf_params={},
                 run_pyzdcf=False, pyzdcf_params={},
                 run_pyroa=False, pyroa_params={},
                 run_mica2=False, mica2_params={},
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

    run_pyroa : bool, optional
        Whether to run the PyROA module. Default is ``False``.

    pyroa_params : dict, optional
        The parameters to pass to the PyROA module. Default is ``{}``.
        
    run_mica2 : bool, optional
        Whether to run the MICA2 module. Default is ``False``.
        
    mica2_params : dict, optional
        The parameters to pass to the MICA2 module. Default is ``{}``.



    Returns
    -------

    output : dict
        A dictionary containing the output of the pipeline modules.

    """

    if mica2_exist:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
    else:
        comm = None
        rank = 0


    if not detrend_exist:
        run_detrend = False
    if not mica2_exist:
        run_mica2 = False

    output_dir = os.path.abspath(output_dir) + r'/'

    if arg2 is None:
        print_error('ERROR: Please provide a list of light curve filenames or the light curves themselves')
        raise Exception


    if not isinstance(arg2[0], str):
        os.makedirs( output_dir + 'input_lcs/', exist_ok=True )
        fnames = []

        for i in range( len(arg2) ):

            if i == 0:
                name = 'continuum'
            else:
                name = 'line{}'.format(i+1)


            fnames.append( output_dir + 'input_lcs/' + name + '.dat' )        
            if rank == 0:
                write_data( arg2[i], output_dir + 'input_lcs/' + name + '.dat' )

        fnames = np.array(fnames)
        kwargs['file_fmt'] = 'csv'
    else:
        fnames = arg2


    cont_fname = fnames[0]
    line_fnames = fnames[1:]

    if type(line_fnames) is str:
        line_fnames = [line_fnames]


    #Read in general kwargs
    general_kwargs = defaults.set_general(kwargs, fnames)
    
    if rank == 0:
        if general_kwargs['verbose']:
            print_header('RUNNING PYPETAL')
        
    if len(fnames) < 2:
        if rank == 0:
            print_error('ERROR: Requires at least two light curves to run pipeline.')
    
        return {}

    if len(line_names) != len(fnames):
        if rank == 0:
            print_error('ERROR: Must have the same number of line names as light curves.')
    
        return {}


    #Get "reject_data"
    _, _, _, _, _, _, reject_data, _, _ = defaults.set_drw_rej(drw_rej_params, fnames)

    #Get "together_pyroa"
    _, _, _, _, _, _, _, _, together_pyroa, _, _, _ = defaults.set_pyroa(pyroa_params, len(fnames))
    
    #Get "together_mica2"
    _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, together_mica2, no_order_mica2 = defaults.set_mica2(mica2_params)


    #Name lines if unnamed
    if line_names is None:
        line_names = np.zeros( len(line_fnames), dtype=str )
        for i in range(len(line_fnames)):
            line_names.append('Line {}'.format(i+1))


    #Get directory common to all files
    input_dir = os.path.dirname( os.path.realpath(cont_fname) )
    for i in range(len(line_fnames)):
        assert input_dir == os.path.dirname( os.path.realpath(line_fnames[i]) )


    if rank == 0:
        #Create subdirectories for each line and module
        make_directories(output_dir, fnames, line_names,
                        run_drw_rej, run_pyccf, run_pyzdcf, run_pyroa, run_mica2,
                        reject_data, together_pyroa, together_mica2, no_order_mica2)


    if general_kwargs['file_fmt'] != 'csv':
        
        input_dir = output_dir + 'temp/'
        if rank == 0:
            #Make temp directory to store csv files
            os.makedirs( output_dir + 'temp/', exist_ok=True )

        #Read in file
        try:
            dat = Table.read( cont_fname, format=general_kwargs['file_fmt'])
        except:
            dat = Table.read( cont_fname, format='ascii')


        #Get filename
        fname = os.path.basename( cont_fname )

        #Write to csv
        cont_fname = input_dir + fname
        if rank == 0:
            write_data( [ dat[ dat.colnames[0] ], dat[ dat.colnames[1] ], dat[ dat.colnames[2] ] ], input_dir + fname )


        for i in range(len(line_fnames)):

            #Read in file
            try:
                dat = Table.read( line_fnames[i], format=general_kwargs['file_fmt'])
            except:
                dat = Table.read( line_fnames[i], format='ascii')

            #Get filename
            fname = os.path.basename( line_fnames[i] )

            #Write to csv
            line_fnames[i] = input_dir + fname
            if rank == 0:
                write_data( [ dat[ dat.colnames[0] ], dat[ dat.colnames[1] ], dat[ dat.colnames[2] ] ], input_dir + fname )

    if rank == 0:

        drw_rej_res = {}
        detrend_res = {}
        pyccf_res = {}
        pyzdcf_res = {}
        pyroa_res = {}
        mica2_res = {}

        if run_drw_rej:

            if np.any( reject_data ):
                os.makedirs( output_dir + 'processed_lcs/', exist_ok=True )

            drw_rej_res = drw_rej_tot( cont_fname, line_fnames, line_names, output_dir, general_kwargs, drw_rej_params )

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
            pyccf_res = pyccf_tot(cont_fname, line_fnames, line_names, output_dir, general_kwargs, pyccf_params)

        if run_pyzdcf:
            pyzdcf_res, plike_res = pyzdcf_tot(cont_fname, line_fnames, line_names, output_dir, general_kwargs, pyzdcf_params)

        if run_pyroa:
            pyroa_res = pyroa_tot(cont_fname, line_fnames, line_names, output_dir, general_kwargs, pyroa_params)


    if run_mica2:
        cont_fname = comm.bcast(cont_fname, root=0)
        line_fnames = comm.bcast(line_fnames, root=0)
        line_names = comm.bcast(line_names, root=0)
        output_dir = comm.bcast(output_dir, root=0)
        general_kwargs = comm.bcast(general_kwargs, root=0)
        mica2_params = comm.bcast(mica2_params, root=0)

        mica2_res = mica2_tot(cont_fname, line_fnames, line_names, output_dir, general_kwargs, mica2_params, comm, rank)


    if rank == 0:

        if general_kwargs['file_fmt'] != 'csv':
            import shutil
            shutil.rmtree(input_dir)


        #Compile all results into a single dict
        tot_res = {}

        if run_drw_rej:
            tot_res['drw_rej_res'] = drw_rej_res

        if run_detrend:
            tot_res['detrend_res'] = detrend_res

        if run_pyccf:
            tot_res['pyccf_res'] = pyccf_res

        if run_pyzdcf:
            tot_res['pyzdcf_res'] = pyzdcf_res
            tot_res['plike_res'] = plike_res

        if run_pyroa:
            tot_res['pyroa_res'] = pyroa_res
            
        if run_mica2:
            tot_res['mica2_res'] = mica2_res


        if not isinstance(arg2[0], str):
            import shutil
            shutil.rmtree( output_dir + 'input_lcs/' )
            
        if general_kwargs['verbose']:
            print_header('RUN FINISHED')


        return tot_res





def run_weighting(output_dir, line_names,
                 run_pyccf=False, pyccf_params={},
                 run_pyroa=False, pyroa_params={},
                 run_mica2=False, mica2_params={},
                 run_javelin=False, javelin_params={},
                 weighting_params={},
                 **kwargs):


    """Perform weighting on the results from the pyPetal pipeline.
    NOTE: This requires that you have run the pyPetal pipeline (i.e. pypetal.pipeline.run.pipeline) with
    at least run_pyccf=True.
    This will perform the weighting scheme described in Grier et al. (2019) to the lag distributions from
    each of the modules run (e.g., pyCCF, PyROA, JAVELIN, MICA2).


    Parameters
    ----------
    output_dir : str
        The directory where the output files will be saved. Must be the same output directory as used
        in pyPetal and pyPetal-jav.

    line_names : list of str
        The names of the lines to be weighted. Must be the same as the line names used in pyPetal and
        pyPetal-jav.

    run_pyccf : bool, optional
        Whether or not pyCCF was run in the pyPetal pipeline. Default is False.

    pyccf_params : dict, optional
        The parameters used in the pyPetal pipeline for pyCCF. Default is {}.

    run_pyroa : bool, optional
        Whether or not PyROA was run in the pyPetal pipeline. Default is False.

    pyroa_params : dict, optional
        The parameters used in the pyPetal pipeline for PyROA. Default is {}.
        
    run_mica2 : bool, optional
        Whether or not MICA2 was run in the pyPetal pipeline. Default is False.
        
    mica2_params : dict, optional
        The parameters used in the pyPetal pipeline for MICA2. Default is {}.

    run_javelin : bool, optional
        Whether or not the pyPetal-jav pipeline was run. Default is False.

    javelin_params : dict, optional
        The parameters used in the pyPetal-jav pipeline. Default is {}.

    weighting_params : dict, optional
        The parameters to be used for the weighting. Default is {}.



    Returns
    -------

    res : dict
        A dictionary containing the results of the weighting.

    """


    if (not run_pyccf) & (not run_javelin) & (not run_pyroa) & (not run_mica2):
        print_error('ERROR: Either JAVELIN, pyCCF, PyROA, or MICA2 must be run before weighting can be done.')
        raise Exception

    output_dir = os.path.abspath(output_dir) + r'/'

    if output_dir + 'processed_lcs/' in glob.glob( output_dir + '*/' ):
        fnames = [ output_dir + 'processed_lcs/' + x + '_data.dat' for x in line_names ]
    else:
        fnames = [ output_dir + 'light_curves/' + x + '.dat' for x in line_names ]

    #Read in general kwargs
    general_kwargs = defaults.set_general(kwargs, fnames)
    if general_kwargs['verbose']:
        print_header('RUNNING PYPETAL WEIGHTING')
    

    #Get "interp"
    interp, _, _, _, _, _ = defaults.set_pyccf(pyccf_params)

    #Get "together" for javelin
    if 'together' in javelin_params:
        together_jav = javelin_params['together']
    else:
        together_jav = False

    #Get "together" for pyroa
    _, _, _, _, _, _, _, _, together_pyroa, _, _, _ = defaults.set_pyroa( pyroa_params, len(line_names) )

    #Get "together" for mica2
    _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, together_mica2, _ = defaults.set_mica2(mica2_params)


    javelin_chain_fnames = None
    pyccf_iccf_fnames = None
    pyccf_dist_fnames = None
    pyroa_sample_fnames = None
    mica2_sample_fnames = None

    if run_pyccf:
        pyccf_iccf_fnames = []
        pyccf_dist_fnames = []

        for x in line_names[1:]:
            pyccf_iccf_fnames.append( output_dir + x + '/pyccf/' + x + '_ccf.dat' )
            pyccf_dist_fnames.append( output_dir + x + '/pyccf/' + x + '_ccf_dists.dat' )


    if run_javelin:
        if together_jav:
            javelin_chain_fnames = output_dir + 'javelin/chain_rmap.txt'
        else:
            javelin_chain_fnames = []

            for x in line_names[1:]:
                javelin_chain_fnames.append( output_dir + x + '/javelin/chain_rmap.txt' )



    if run_pyroa:
        if together_pyroa:
            pyroa_sample_fnames = output_dir + 'pyroa/samples.obj'
        else:
            pyroa_sample_fnames = []

            for x in line_names[1:]:
                pyroa_sample_fnames.append( output_dir + x + '/pyroa/samples.obj' )
                
                
    
    if run_mica2:
        if together_mica2:
            mica2_sample_fnames = []
            
            for x in line_names[1:]:
                mica2_sample_fnames.append( output_dir + x + '/mica2/{}_lag_samples.dat'.format(x) )
            
        else:
            mica2_sample_fnames = []

            for x in line_names[1:]:
                mica2_sample_fnames.append( output_dir + x + '/mica2/{}_lag_samples.dat'.format(x) )



    res = run_weighting_tot(output_dir,
                                      javelin_chain_fnames, pyccf_iccf_fnames, pyccf_dist_fnames, pyroa_sample_fnames, mica2_sample_fnames,
                                      line_names, interp, together_jav,
                                      pyroa_params=pyroa_params, general_kwargs=general_kwargs,
                                      weighting_params=weighting_params)

    if general_kwargs['verbose']:
        print_header('RUN FINISHED')


    return res
