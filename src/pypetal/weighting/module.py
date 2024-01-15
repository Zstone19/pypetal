import glob
import os
import warnings

import numpy as np

from pypetal.utils import defaults
from pypetal.utils.load import get_ordered_line_names
from pypetal.weighting.plotting import plot_weight_output, plot_weights
from pypetal.weighting.utils import (combine_weight_outputs,
                                     run_weighting_single)


def run_weighting_tot(output_dir,
                      jav_chain_fnames=None, pyccf_iccf_fnames=None, pyccf_dist_fnames=None,
                      pyroa_sample_fnames=None,
                      line_names=None, interp=2, together_jav=False,
                      pyroa_obj_inds=None, pyroa_params={},
                      general_kwargs={}, weighting_params={}, share_lag_bounds=True):


    """Run the weighting module for a number of lines, with pyPetal (and possible pyPetal-jav) already run.
    This will use the weighting scheme described in Grier et al. (2019) for all lag distributions provided.
    NOTE: If filenaes aren't given for a certain module (i.e. its input is None), it will be assumed not to have run.


    Parameters
    ----------

    output_dir : str
        The directory where the output will be saved. Must be the same output directory used for the pyPetal run.

    jav_chain_fnames : list of str, optional
        The filenames of the JAVELIN chains for each line. If None, JAVELIN will be assumed not to have run.
        Default is None.

    pyccf_iccf_fnames : list of str, optional
        The filenames of the PyCCF ICCF outputs for each line. If None, PyCCF will be assumed not to have run.
        Default is None.

    pyccf_dist_fnames : list of str, optional
        The filenames of the PyCCF lag distributions for each line. If None, PyCCF will be assumed not to have run.
        Default is None.

    pyroa_sample_fnames : list of str, optional
        The filenames of the PyROA sample for each line (i.e. "samples.obj"). If None, PyROA will be assumed not to have run.
        Default is None.

    line_names : list of str, optional
        The names of the lines to be weighted. If None, the output directory will be searched to determine the
        line_names. Default is None.

    interp : int, optional
        The interpolation factor to use for the PyCCF lag distributions. Default is 2.

    together_jav : bool, optional
        The together parameter used for the JAVELIN run. Default is False.

    pyroa_obj_inds : list of int, optional
        The indices in the samples_chunked array for each line (in the same order as line_names).
        Typically, the index for a given line will be one more than its index in line_names.
        Must be specified if pyroa_sample_fnames is not None. Default is None.

    pyroa_params : dict, optional
        The parameters used for the PyROA run. Default is {}.

    general_kwargs : dict, optional
        The general kwargs used for the pyPetal run. Default is {}.

    weighting_params : dict, optional
        The parameters to be used for the weighting. Default is {}.

    share_lag_bounds : bool, optional
        Whether to share the lag bounds for all lines when computing the weights. Default is True.


    Returns
    -------

    res_tot : list of dict
        A list of dictionaries containing the results for each line.

    summary_dicts : list of dict
        A list of dictionaries containing the summary information for each line. This information will
        be stored in the weight summary file in the output directory as well.

    """

    output_dir = os.path.abspath(output_dir) + r'/'
    pyccf_params = {'interp':interp}

    if line_names is None:
        warnings.warn('Assuming that the filenames are in the same order as the line names. Line names will be acquired in chronological order from the given directory, except the first will be the continuum', RuntimeWarning)
        line_names = get_ordered_line_names(output_dir)

    _, _, _, _, _, _, _, _, together_pyroa, _, _ = defaults.set_pyroa( pyroa_params, len(line_names) )

    #---------------------------
    #Get data fnames

    if output_dir + 'processed_lcs/' in glob.glob( output_dir + '*' ):
        line_fnames = np.array([ output_dir + 'processed_lcs/' + x + '_data.dat' for x in line_names ])
    else:
        line_fnames = np.array([ output_dir + 'light_curves/' + x + '.dat' for x in line_names ])


    general_kwargs = defaults.set_general(general_kwargs, line_fnames)
    _, _, _, _, zoom = defaults.set_weighting(weighting_params)

    #---------------------------
    #Share lag bounds?
    if share_lag_bounds:

        baselines = []

        x_cont, _, _ = np.loadtxt(line_fnames[0], unpack=True, delimiter=',', usecols=[0,1,2])
        for i in range(len(line_fnames)):
            x_line, _, _ = np.loadtxt(line_fnames[i], unpack=True, delimiter=',', usecols=[0,1,2])

            bl = np.max([ x_cont.max(), x_line.max() ]) - np.min([ x_cont.min(), x_line.min() ])
            baselines.append(bl)



        lag_bounds = []
        lag_bounds_i = [-np.max(baselines), np.max(baselines)]
        for i in range(len(line_fnames)-1):
            lag_bounds.append(lag_bounds_i)


    else:
        lag_bounds = general_kwargs['lag_bounds']


    #---------------------------
    #Account for None inputs
    run_javelin = True
    run_pyccf = True
    run_pyroa = True

    if jav_chain_fnames is None:
        warnings.warn('Assuming JAVELIN was not run.', RuntimeWarning)

    if (pyccf_iccf_fnames is None) or (pyccf_iccf_fnames is None):
        warnings.warn('Assuming PyCCF was not run.', RuntimeWarning)


    if jav_chain_fnames is None:
        jav_chain_fnames = [None] * ( len(line_names)-1)
        run_javelin = False
    elif together_jav:
        jav_chain_fnames = [jav_chain_fnames] * ( len(line_names)-1)
        run_javelin = True


    if pyccf_iccf_fnames is None:
        pyccf_iccf_fnames = [None] * ( len(line_names)-1)
        run_pyccf = False

    if pyccf_dist_fnames is None:
        pyccf_dist_fnames = [None] * ( len(line_names)-1)
        run_pyccf = False



    #Look for PyROA files if none input
    if pyroa_sample_fnames is None:
        warnings.warn('Looking for PyROA files, no files input', RuntimeWarning)


        if together_pyroa:
            fnames_i = glob.glob(output_dir + 'pyroa/*.obj')
            if len(fnames_i) == 0:
                warnings.warn('No PyROA files found. Assuming PyROA was not run.', RuntimeWarning)
                run_pyroa = False
            else:
                pyroa_sample_fnames = output_dir + 'pyroa/samples.obj'


        else:

            pyroa_sample_fnames = []
            for i in range(len(line_fnames[1:])):
                fnames_i = glob.glob(output_dir + r'/' + line_names[i+1] + 'pyroa/*.obj')

                if len(fnames_i) == 0:
                    warnings.warn('No PyROA files found for ' + line_names[i+1], RuntimeWarning)
                    pyroa_sample_fnames.append( None )

                else:
                   pyroa_sample_fnames.append( output_dir + r'/' + line_names[i+1] + 'pyroa/samples.obj' )

            if np.all([ x is None for x in pyroa_sample_fnames ]):
                warnings.warn('No PyROA files found. Assuming PyROA was not run.', RuntimeWarning)
                run_pyroa = False
                pyroa_sample_fnames = None


    #Check if PyROA files exist
    if isinstance(pyroa_sample_fnames, list):
        for i in range(len(pyroa_sample_fnames)):

            if os.path.isfile(pyroa_sample_fnames[i]):
                continue
            else:
               warnings.warn( pyroa_sample_fnames[i] + ' not found.', RuntimeWarning )
               pyroa_sample_fnames[i] = None

        if np.all([ x is None for x in pyroa_sample_fnames ]):
            warnings.warn('No PyROA files found. Assuming PyROA was not run.', RuntimeWarning)
            run_pyroa = False
            pyroa_sample_fnames = None


    elif isinstance(pyroa_sample_fnames, str):
        if not os.path.isfile(pyroa_sample_fnames):
            warnings.warn( pyroa_sample_fnames + ' not found. Assuming PyROA was not run.', RuntimeWarning )
            run_pyroa = False
            pyroa_sample_fnames = None


    if pyroa_sample_fnames is None:
        pyroa_sample_fnames = [None] * ( len(line_names)-1)
        run_pyroa=False
    elif together_pyroa:
        pyroa_sample_fnames = [pyroa_sample_fnames] * ( len(line_names)-1)
        run_pyroa = True


    if run_pyroa:
        if pyroa_obj_inds is None:
            if together_pyroa:
                pyroa_obj_inds = range(1, len(line_names))
            else:
                pyroa_obj_inds = []

                for i in range(len(line_names)-1):
                    if pyroa_sample_fnames[i] is None:
                        pyroa_obj_inds.append(None)
                    else:
                        pyroa_obj_inds.append(1)
    else:
        pyroa_obj_inds = [None] * ( len(line_names)-1)

    #---------------------------
    #Make weights directories
    for name in line_names[1:]:
        os.makedirs(output_dir + name + r'/weights/', exist_ok=True)

    #---------------------------
    #Run weighting

    summary_dicts = []
    outputs = []

    summary_fnames = []


    for i in range(len(line_fnames)-1):

        if together_jav:
            javelin_lag_col = 2 + 3*i
        else:
            javelin_lag_col = 2

        res, summary_dict = run_weighting_single(output_dir + line_names[i+1] + r'/weights/',
                                                 line_fnames[0], line_fnames[i+1],
                                                 weighting_params, lag_bounds[i],
                                                 jav_chain_fnames[i], pyccf_iccf_fnames[i], pyccf_dist_fnames[i],
                                                 pyroa_sample_fnames[i],
                                                 javelin_lag_col=javelin_lag_col,
                                                 pyroa_obj_ind=pyroa_obj_inds[i],
                                                 pyccf_params=pyccf_params, pyroa_params=pyroa_params)

        summary_dicts.append(summary_dict)
        outputs.append(res)
        summary_fnames.append(output_dir + name + r'/weights/weight_summary.fits')



    #---------------------------
    #Get total results

    res_tot = combine_weight_outputs(outputs, run_pyccf, run_javelin, run_pyroa)

    plot_mod = 'pyccf'
    if not run_pyccf:
        if run_javelin:
            plot_mod = 'javelin'
        elif run_pyroa:
            plot_mod = 'pyroa'

    for i in range(len(line_fnames)-1):
        plot_weights(output_dir, line_names[i+1], res_tot[plot_mod],
                        summary_dicts[i]['n0_'+plot_mod], summary_dicts[i]['k'],
                        general_kwargs['time_unit'], general_kwargs['plot'])

    if run_pyccf:
        plot_weight_output(output_dir, line_fnames[0], line_fnames, line_names,
                            res_tot['pyccf'], general_kwargs, 'pyccf', zoom,
                            general_kwargs['plot'])

    if run_javelin:
        plot_weight_output(output_dir, line_fnames[0], line_fnames, line_names,
                            res_tot['javelin'], general_kwargs, 'javelin', zoom,
                            general_kwargs['plot'])

    if run_pyroa:
        plot_weight_output(output_dir, line_fnames[0], line_fnames, line_names,
                            res_tot['pyroa'], general_kwargs, 'pyroa', zoom,
                            general_kwargs['plot'])


    return res_tot, summary_dicts
