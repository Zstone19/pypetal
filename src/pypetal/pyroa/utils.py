import glob
import multiprocessing as mp
import os
import pickle
import shutil
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from os import devnull

import numpy as np
import PyROA
from astropy.table import Table

##############################################################
####################### SILENCE OUTPUT #######################
##############################################################


@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)

##############################################################
######################### FIT CLASS ##########################
##############################################################

class MyFit:

    def __init__(self, output_dir):

        self.output_dir = output_dir

        with open(output_dir + 'samples.obj', 'rb') as f:
            self.samples = pickle.load(f)

        with open(output_dir + 'samples_flat.obj', 'rb') as f:
            self.samples_flat = pickle.load(f)

        with open(output_dir + 'Lightcurve_models.obj', 'rb') as f:
            self.models = pickle.load(f)

        with open(output_dir + 'X_t.obj', 'rb') as f:
            dat = pickle.load(f)

            self.t = dat[0]
            self.X = dat[1]
            self.X_errs = dat[2]


##############################################################
######################## ASSIST FUNCTIONS ####################
##############################################################

def get_samples_chunks(samples, nburn=0, add_var=False, delay_dist=False):

    """Split the samples up into individual lines. The burn-in samples will be removed and the walkers will be flattened.


    Parameters
    ----------

    samples : list of float
        The samples object output from PyROA. Can be found in the "samples.obj" file.

    nburn : int, optional
        The number of burn-in samples to discard. Default is 0.



    Returns
    -------

    samples_chunks : list of float
        The chunked sample data. If the input ``samples`` array has a shape (nchain, nwalker, nobj*nparam + extra), the output ``samples_chunked`` array will
        have a shape ( nobj+1, (nchain-nburn)*nwalker, nparam ).

    """

    chunk_size = 3
    if add_var:
        chunk_size += 1

    #Remove burn-in
    samples = samples[nburn:,:,:].copy()

    #Insert tau=0 for the continuum
    samples = np.insert( samples, 2, np.zeros( (samples.shape[0], samples.shape[1]) ), axis=2 )

    #Insert delta_i=0 for the continuum
    if delay_dist:
        chunk_size += 1
        samples = np.insert(samples, 3, np.zeros( (samples.shape[0], samples.shape[1]) ), axis=2 )

    #Flatten
    samples = samples.reshape( (-1, samples.shape[2]) ).T

    #Chunk
    nlc = (samples.shape[0] - 1)//chunk_size
    samples_chunks = []
    for i in range(nlc+1):
        samples_chunks.append( samples[i*chunk_size:(i+1)*chunk_size] )

    return samples_chunks



def save_lines(line_fnames, line_names, output_dir, objname=None, delimiter=None, subtract_mean=False, div_mean=False):

    """Save light curves as PyROA-readable files. Light curves will be saved
    as ASCII files with the format: x, y, yerr. The saved files will be named "{objname}_{line_name}.dat".



    Parameters
    ----------

    line_fnames : list of str
        The filenames of the light curves to be saved.

    line_names : list of str
        The names of the light curves to be saved.

    output_dir : str
        The directory to save the light curves to.

    objname : str, optional
        The name of the object. This will Default is 'pyroa'.

    delimiter : str, optional
        The delimiter to use when reading the light curves. Default is None.

    subtract_mean : bool, optional
        If True, the mean of the light curve will be subtracted. Default is False.

    div_mean : bool, optional
        If True, the light curve will be divided by the mean. This occurs before subtracting the mean.
        Default is False.




    Returns
    -------

    new_fnames : list of str
        The filenames of the saved light curves.

    """

    if objname is None:
        objname = 'pyroa'

    new_fnames = []

    for i in range(len(line_fnames)):
        x, y, yerr = np.loadtxt( line_fnames[i], unpack=True, usecols=[0,1,2], delimiter=delimiter )

        if div_mean:
            yerr /= np.mean(y)
            y /= np.mean(y)

        if subtract_mean:
            y -= np.mean(y)

        new_fnames.append( output_dir + objname + '_' + line_names[i] + '.dat' )
        np.savetxt( new_fnames[-1], np.array([x,y,yerr]).T )

    return new_fnames




def get_priors(fnames, laglim_in, subtract_mean=False, div_mean=False, together=False, delimiter=','):


    """Get the priors used for PyROA.


    Parameters
    ----------

    fnames : list of str
        The filenames of the light curves.

    laglim_in : list of float
        The bounds to search between for each light curve. This should be a list of lists,
        each list having an upper and lower lag bound.

    subtract_mean : bool, optional
        If True, the mean of the light curve will be subtracted. Default is False.

    div_mean : bool, optional
        If True, the light curve will be divided by the mean. This occurs before subtracting the mean.
        Default is False.

    together : bool, optional
        Whether or not to fit all light curves to the continuum in one fit. Default is False.




    Returns
    -------

    prior_arr : list of float
        The array of priors used for PyROA.

    """


    if together:

        lower_lim = []
        upper_lim = []
        for i in range(len(fnames)-1):
            lower_lim.append( laglim_in[i][0] )
            upper_lim.append( laglim_in[i][1] )

        laglim = [ np.min(lower_lim), np.max(upper_lim) ]



        std_vals = []
        min_y = []
        max_y = []
        for i in range(len(fnames)):
            _, y, yerr = np.loadtxt( fnames[i], unpack=True, usecols=[0,1,2], delimiter=delimiter )

            if div_mean:
                yerr /= np.mean(y)
                y /= np.mean(y)

            if subtract_mean:
                y -= np.mean(y)

            std_vals.append( np.std(y) )
            min_y.append( np.min(y-yerr) )
            max_y.append( np.max(y+yerr) )


        if div_mean:
            a_prior = [0., 2.]
            b_prior = [0., 2.]
            err_prior = [0., 10.]
        else:
            a_prior = [0., np.max(std_vals)]
            b_prior = [np.min(min_y), np.max(max_y)]
            err_prior = [0., 10*np.max(std_vals)]

        tau_prior = laglim
        delta_prior = [5., 50.]


        return [a_prior, b_prior, tau_prior, delta_prior, err_prior]


    else:
        prior_arr = np.zeros(( len(fnames)-1, 5, 2 ))

        _, y_cont, yerr_cont = np.loadtxt( fnames[0], unpack=True, usecols=[0,1,2], delimiter=delimiter )
        if div_mean:
            yerr_cont /= np.mean(y_cont)
            y_cont /= np.mean(y_cont)

        if subtract_mean:
            y_cont -= np.mean(y_cont)


        for i in range(len(fnames)-1):
            _, y, yerr = np.loadtxt( fnames[i+1], unpack=True, usecols=[0,1,2], delimiter=',' )

            if div_mean:
                yerr /= np.mean(y)
                y /= np.mean(y)

                if subtract_mean:
                    y -= np.mean(y)

                #a
                prior_arr[i,0,0] = 0.
                prior_arr[i,0,1] = 2.

                #b
                prior_arr[i,1,0] = 0.
                prior_arr[i,1,1] = 2.

                #err
                prior_arr[i,4,0] = 0.
                prior_arr[i,4,1] = 10.

            else:

                if subtract_mean:
                    y -= np.mean(y)

                #a - RMS of the LC
                prior_arr[i,0,0] = 0.
                prior_arr[i,0,1] = 10.*np.max( [np.std(y), np.std(y_cont)] )

                #b - mean of the LC
                prior_arr[i,1,0] = np.min( [np.min(y-yerr), np.min(y_cont-yerr_cont)] )
                prior_arr[i,1,1] = np.max( [np.max(y+yerr), np.max(y_cont+yerr_cont)] )

                #err - extra error
                prior_arr[i,4,0] = 0.
                prior_arr[i,4,1] = 10.*np.max([ np.std(y), np.std(y_cont) ])


            #tau
            prior_arr[i,2,0] = laglim_in[i][0]
            prior_arr[i,2,1] = laglim_in[i][1]

            #delta - window function width
            prior_arr[i,3,0] = 5.
            prior_arr[i,3,1] = 50.

        return prior_arr


def move_output_files(file_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    file_dir = os.path.abspath(file_dir) + r'/'
    output_dir = os.path.abspath(output_dir) + r'/'

    files = glob.glob(file_dir + '*.obj')
    for f in files:
        shutil.move(f, output_dir)

    return


def samples2table(samples, line_names, nburn=0, add_var=False, delay_dist=False):

    """Convert the PyROA "samples.obj" file into an astropy.table.Table object.


    Parameters
    ----------

    samples : list of float
        The samples from the PyROA run.

    line_names : list of str
        The names of the light curves.

    nburn : int, optional
        The number of burn-in samples to discard. Default is 0.

    add_var : bool, optional
        Whether the PyROA run included an extra error term. Default is False.

    delay_dist : bool, optional
        Whether the PyROA run included a delay distribution. Default is False.



    Returns
    -------

    table_dict : dict of astropy.table.Table
        A dictionary of astropy.table.Table objects, one for each light curve. There will also be
        a Table for the windowing function width, under the key 'delta'.

    """


    #Works for together=True
    samples_chunks = get_samples_chunks(samples, nburn)
    table_dict = {}

    for i, name in enumerate(line_names[1:]):
        samp_i = samples_chunks[i]
        tab_i = Table()

        tab_i.add_column( 'A', samp_i[:,0] )
        tab_i.add_column( 'B', samp_i[:,1] )
        tab_i.add_column( 'tau', samp_i[:,2] )

        if add_var:
            tab_i.add_column( 'sigma', samp_i[:,3] )

        if delay_dist:
            tab_i.add_column( 'delta', samples_chunks[-1][:,4] )

        table_dict[name] = tab_i


    delta_tab = Table()
    delta_tab.add_column( 'delta', samples_chunks[-1][:,0] )
    table_dict['delta'] = delta_tab

    return table_dict


def get_all_samples(samples, nburn=0):

    #Remove burn-in
    samples = samples[nburn:,:,:].copy()

    #Flatten
    samples = samples.reshape( (-1, samples.shape[2]) ).T

    return samples


##############################################################
##############################################################
##############################################################

def run_pyroa(fnames, lc_dir, line_dir, line_names,
              nburn=10000, nchain=15000, lag_bounds=None,
              init_tau=None, init_delta=10, sig_level=100,
              together=True, subtract_mean=True, div_mean=False,
              add_var=False, delay_dist=False, psi_types='Gaussian',
              objname=None, verbose=True):


    """Run PyROA for a number of input light curves.
    NOTE: This will assume a log-Gaussian for the delay distribution.


    Parameters
    ----------

    fnames : list of str
        A list of paths to the light curves to be used in PyROA. The first light curve will be assumed to be the continuum. Must be in CSV format.

    lc_dir : str
        The output directory to put the light curves for PyROA to use.

    line_dir : str, list of str
        The output directory to store the data products output from PyROA.
        If together=False, this should be a list of paths for each line. If together=True, this should be a single path.

    line_names : list of str
        A list of the line names corresponding to the light curves.

    nburn : int, optional
        The number of burn-in samples to discard. Default is 10000.

    nchain : int, optional
        The number of samples to get per walker. Default is 15000.

    lag_bounds : list of str, optional
        The range of lags to consider for PyROA, one set of bounds for each line (excluding the continuum). If ``None``, the baseline of
        each light curve will be used as the lag bounds. Default is ``None``.

    init_tau : list of float
        A list of inital values for the time delay (tau) for each light curve. If ``None``, will set the inital value to 10.0 for each line.
        Default is ``None``.

    sig_level : float
        The number to use for sigma clipping in PyROA. Default is 100.

    together : bool, optional
        Whether to fit all light curves together or not. If ``together=False``, the ``line_dirs`` argument must be set. Default is ``True``.

    subtract_mean : bool, optional
        Whether to subtract the mean from all light curves before using PyROA or not. Will occur after ``div_mean`` if set to ``True``. Default is ``True``.

    div_mean : bool, optional
        Whether to divide each light curve by its mean before using PyROA. Will occur before ``subtract_mean`` if set to ``True``. Default is ``False``.

    add_var : bool or list of bool, optional
        Whether or not to add additional uncertainty in the data, same as the PyROA argument. If ``together=False``, multiple values may be given for each line.
        If only one value is given, it will be assumed for all lines. Default is ``True``.

    delay_dist : bool or list of bool, optional
        Same as the ``delay_dist`` argument for PyROA. If ``together=False``, multiple values may be given for each line.
        If only one value is given, it will be assumed for all lines. Default is ``True``.

    objname : str, optional
        The name of the object, will be used for plot and for the saved PyROA light curve data. If ``None``, will be set to "pyroa".

    line_names : list of str
        A list of directories to place the output PyROA data for each of the lines (excluding the continuum). Must be set if ``together=False``, and will only be used in such a case.
        Default is ``None``.


    Returns
    -------

    fit : PyROA.Fit or list of pyROA.Fit
        The PyROA.Fit object output from PyROA. If ``together=False``, this will be an array of the ``Fit`` objects.


    """


    if objname is None:
        objname = 'pyroa'

    if init_tau is None:
        init_tau = np.full( len(fnames)-1, 10. )

    if isinstance(psi_types, str):
        psi_types = [psi_types] * ( len(fnames)-1 )

    #If no lag bounds given, use the baseline of the light curves
    if lag_bounds is None:
        lag_bounds = []

        x_cont, _, _ = np.loadtxt( fnames[0], unpack=True, usecols=[0,1,2], delimiter=',' )
        for i in range(1, len(fnames)):
            x, _, _ = np.loadtxt( fnames[i], unpack=True, usecols=[0,1,2], delimiter=',' )

            max_x = np.max([ np.max(x), np.max(x_cont) ])
            min_x = np.min([ np.min(x), np.min(x_cont) ])
            bl = max_x - min_x

            lag_bounds.append([-bl, bl])


    if (not together) & (len(fnames)==2) & isinstance(line_dir, str):
        line_dir = [line_dir]


    if isinstance(line_dir, str):
        os.makedirs(line_dir, exist_ok=True)
    elif isinstance(line_dir, list):
        for ldir in line_dir:
            os.makedirs(ldir, exist_ok=True)
    else:
        raise ValueError('line_dir must be a string or list of strings')


    os.makedirs(lc_dir, exist_ok=True)

    prior_arr = get_priors(fnames, lag_bounds, subtract_mean=subtract_mean, div_mean=div_mean, together=together, delimiter=',')
    _ = save_lines(fnames, line_names, lc_dir, objname=objname, subtract_mean=subtract_mean, div_mean=div_mean, delimiter=',')

    cwd = os.getcwd()

    if not together:

        assert isinstance(line_dir, list), 'Must provide multiple line_dir if together=False'

        if isinstance(add_var, bool):
            add_var = np.full( len(fnames)-1, add_var )

        if isinstance(delay_dist, bool):
            delay_dist = np.full( len(fnames)-1, delay_dist )


        fit_arr = []

        for i in range(len(fnames)-1):

            filters = [line_names[0], line_names[i+1]]

            args = (lc_dir, objname, filters, prior_arr[i,:,:],)
            kwargs = {'add_var':add_var[i], 'init_tau':[init_tau[i]], 'init_delta':init_delta, 'sig_level':sig_level,
                      'delay_dist':delay_dist[i], 'psi_types':[psi_types[i]], 'Nsamples':nchain, 'Nburnin':nburn}

            if verbose:
                proc = mp.get_context('fork').Process(target=PyROA.Fit, args=args, kwargs=kwargs)

                proc.start()
                while proc.is_alive():
                    proc.is_alive()
                proc.terminate()

            else:
                with suppress_stdout_stderr():
                    proc = mp.get_context('fork').Process(target=PyROA.Fit, args=args, kwargs=kwargs)

                    proc.start()
                    while proc.is_alive():
                        proc.is_alive()
                    proc.terminate()


            move_output_files(cwd, line_dir[i])
            fit = MyFit(line_dir[i])

            fit_arr.append(fit)

        return fit_arr

    else:

        assert isinstance(line_dir, str), 'Must provide one line_dir if together=True'


        args = (lc_dir, objname, line_names, prior_arr,)
        kwargs = {'add_var':add_var, 'init_tau':init_tau, 'init_delta':init_delta, 'sig_level':sig_level,
                  'delay_dist':delay_dist, 'psi_types':psi_types, 'Nsamples':nchain, 'Nburnin':nburn}

        if verbose:
            proc = mp.get_context('fork').Process(target=PyROA.Fit, args=args, kwargs=kwargs)

            proc.start()
            while proc.is_alive():
                proc.is_alive()
            proc.terminate()

        else:
            with suppress_stdout_stderr():
                proc = mp.get_context('fork').Process(target=PyROA.Fit, args=args, kwargs=kwargs)

                proc.start()
                while proc.is_alive():
                    proc.is_alive()
                proc.terminate()

        move_output_files(cwd, line_dir)
        fit = MyFit(line_dir)

        return fit
