import numpy as np
from PyROA.PyROA import Gaussian, Uniform, LogGaussian, TruncGaussian, InverseGaussian
import scipy
from astropy.table import Table


import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import ListedColormap
import palettable
import corner


import os
import glob
import shutil




import matplotlib as mpl

mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.size'] = 3

mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.size'] = 3

mpl.rcParams["figure.autolayout"] = False

mpl.rcParams['savefig.dpi'] = 200
mpl.rcParams['savefig.format'] = 'pdf'




###########
# NOTES:

# - add_var adds a sigma parameter for each LC (except cont)
# - delay_dist adds a delta_i parameter for each LC (except cont)

# - samples_flat.obj thins the samples by a factor 15
# - samples_flat.obj uses the emcee "get_chain" function, with flat=True and discard=Nburnin
# - samples.obj uses the emcee "get_chain" function (without discarding)
# - Lightcurve_models.obj contains the fits to each light curve [ [x1, y1, err1], [x2, y2, 2err2], ... ]

# - delay_dist allows for a distribution to blur the lags between LCs.
#   If delay_dist=False, this function is a dirac delta function, which is convolved with the continuum to get the line. 
#   Different distributions can be convolved with the continuum to get the line if delay_dist=True, which is specified by psi_types.
#   This can be: 'Gaussian' (default), 'Uniform', 'LogGaussian', 'TruncGaussian', 'InverseGaussian' 

# - add_var adds an additional source of error for each LC (including the continuum) added in quadrature to the orifinal error.
#   This is a parameter fit for in the MCMC.


############################################################################################################
############################################  UTILS  #######################################################
############################################################################################################

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
            lower_lim.append( laglim_in[0] )
            upper_lim.append( laglim_in[1] )
            
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


############################################################################################################
#########################################  PLOTTING  #######################################################
############################################################################################################

def plot_histograms(samples, line_names, nburn=0, add_var=False, delay_dist=False, nbin=50, fname=None, show=False):


    """Plot histograms of the posterior distributions of the parameters from PyROA.
    
    Parameters
    ----------
    
    samples : list of float
        The samples object output from PyROA. Can be found in the "samples.obj" file.
    
    line_names: list of str
        The names of the lines corresponding to the samples.
        
    nburn : int, optional
        The number of burn-in samples to discard. Default is 0.
        
    add_var : bool, optional
        The add_var parameter from PyROA. Default is False.
        
    delay_dist : bool, optional
        The delay_dist parameter from PyROA. Default is False.
        
    nbin : int, optional
        The number of bins to use in the histograms. Default is 50.
        
    fname : str, optional
        The name of the file to save the plot to. If None, the plot will not be saved. Default is None.
        
    show : bool, optional
        Whether to show the plot. Default is False.
        
        
    Returns
    -------
    
    """


    samples_chunks = get_samples_chunks(samples, nburn, add_var, delay_dist)

    names = ['A', 'B', r'\tau']
    if delay_dist:
        names.append(r'\Delta')
    if add_var:
        names.append(r'\sigma')


    nlc = len(samples_chunks) - 1
    nrow = nlc + 1
    ncol = 3

    if add_var:
        ncol += 1
    
    if delay_dist:
        ncol += 1


    fig, ax = plt.subplots(nrow, ncol, figsize=(5*ncol, 5*nrow) )

    for n in range(nlc):
        ax[n, 0].set_ylabel('N', fontsize=20)

        for i in range(len(samples_chunks[0])):
            ax[n,i].hist(samples_chunks[n][i], nbin)

            ax[n,i].tick_params('y', labelsize=13)
            ax[n,i].tick_params('x', labelsize=15)
            ax[n,i].tick_params('both', which='major', length=7)
            ax[n,i].tick_params('both', which='minor', length=3)

            ax[n,i].set_xlabel( r'$' + names[i] + '_{' + line_names[n] + r'}$', fontsize=20)


    #For delta
    for i in range(ncol):
        if i == 0:
            ax[-1, i].set_ylabel('N', fontsize=20)
            ax[-1,i].set_xlabel( r'$\Delta$', fontsize=20 )
            ax[-1,i].hist(samples_chunks[-1][-1], nbin)

            ax[-1,i].tick_params('y', labelsize=13)
            ax[-1,i].tick_params('x', labelsize=15)
            ax[-1,i].tick_params('both', which='major', length=7)
            ax[-1,i].tick_params('both', which='minor', length=3)

        else:
            ax[-1, i].axis('off')


    plt.subplots_adjust(hspace=.3)
    
    if fname is not None:
        plt.savefig('pyroa_histogram.pdf', bbox_inches='tight')
    
    if show:
        plt.show()
        
    plt.cla()
    plt.clf()
    plt.close()
    
    return




def pyroa_trace_plot(samples, line_names, add_var=False, delay_dist=False, nburn=0, fname=None, show=False, alpha=.1):

    """Plot traces of the walkers from the MCMC sampling in PyROA.
    
    
    Parameters
    ----------
    
    samples : list of float
        The samples object output from PyROA. Can be found in the "samples.obj" file.
        
    line_names: list of str
        The names of the lines corresponding to the samples.
        
    add_var : bool, optional
        The add_var parameter from PyROA. Default is False.
        
    delay_dist : bool, optional
        The delay_dist parameter from PyROA. Default is False.
        
    nburn : int, optional
        The number of burn-in samples to discard. Default is 0.
        
    fname : str, optional
        The name of the file to save the plot to. If None, the plot will not be saved. Default is None.
        
    show : bool, optional
        Whether to show the plot. Default is False.
        
    alpha : float, optional
        The alpha value for the walker traces. Default is 0.1.
    
    
    Returns
    -------
    
    """


    row_labels = [r'$A_i$', r'$B_i$', r'$\tau_i$']
    if delay_dist:
        row_labels.append( r'$\Delta_i$' )
    if add_var:
        row_labels.append( r'$\sigma_i$' )

    nrow = 3
    if add_var:
        nrow += 1
    if delay_dist:
        nrow += 1

    ncol = (samples.shape[2]-3)//nrow + 2
    nparam = samples.shape[2] - 1

    fig, ax = plt.subplots(nrow, ncol, figsize=(5*ncol, nrow), sharex=True)


    #Labels
    for i in range(nrow):
        ax[i,0].set_ylabel(row_labels[i])

    #Do continuum
    for j in range(samples.shape[1]):
        ax[0,0].plot( samples[:,j,0], c='k', alpha=alpha )
        ax[1,0].plot( samples[:,j,1], c='k', alpha=alpha )    
        ax[2,0].plot( np.zeros( samples.shape[0] ), c='k', alpha=alpha )
        if add_var & delay_dist:
            ax[3,0].plot( np.zeros( samples.shape[0] ), c='k', alpha=alpha )
            ax[4,0].plot( samples[:,j,2], c='k', alpha=alpha )
        elif add_var:
            ax[3,0].plot( samples[:,j,2], c='k', alpha=alpha )
        elif delay_dist:
            ax[3,0].plot( np.zeros( samples.shape[0] ), c='k', alpha=alpha )

    if nburn > 0:
        for i in range(nrow):
            ax[i,0].axvline(nburn, ls='--', color='r')        

    ax[0,0].set_title(line_names[0])


    #Do lines
    if add_var: 
        min_i = 3
    else:
        min_i = 2
    
    for i in range(min_i, nparam):
        
        if add_var:
            m = (i-min_i)%nrow
            n = (i-min_i)//nrow + 1 
        else:
            m = (i-min_i)%nrow
            n = (i-min_i)//nrow + 1

        
        for j in range(samples.shape[1]):
            ax[m,n].plot( samples[:,j,i], c='k', alpha=alpha )
            
        if nburn > 0:
            ax[m,n].axvline(nburn, ls='--', color='r')
            
        if m == 0:
            ax[m,n].set_title(line_names[n])
            
            
        
    #Delta
    for j in range(samples.shape[1]):
        ax[0,-1].plot( samples[:,j,-1], c='k', alpha=alpha )
        
    if nburn > 0:
        ax[0,-1].axvline(nburn, ls='--', color='r')
        
    ax[0,-1].set_title(r'$\Delta$')

    for i in range(1, nrow):
        ax[i,-1].axis('off')


    plt.subplots_adjust(wspace=.2)

    if fname is not None:
        plt.savefig(fname, dpi=200, bbox_inches='tight')

    if show:
        plt.show()
        
    plt.cla()
    plt.clf()
    plt.close()
    
    return







def peaktomean(mu, a, rms):
    alpha = (a-mu)/rms
    phi = (1.0/np.sqrt(2.0*np.pi))*np.exp(-0.5*alpha**2)
    Z = 1 - 0.5*(1+scipy.special.erf(alpha/np.sqrt(2.0)))
    return mu + rms*phi/Z  



def plot_fits(fnames, line_names, samples, models, 
              nburn=0, add_var=False, delay_dist=False, psi_types='Gaussian', 
              nbin=50, delimiter=',', 
              time_unit='d', lc_unit='mJy', 
              output_fname=None, show=False):
    
    if isinstance(psi_types, str):
        psi_types = [psi_types]*( len(fnames)-1 )
        
    if isinstance(lc_unit, str):
        lc_unit = [lc_unit]*len(fnames)
        
    #Get colors
    cmap = ListedColormap( palettable.cartocolors.qualitative.Vivid_10.mpl_colors )
    colors = cmap.colors
        
    #Get samples
    samples_chunks = get_samples_chunks(samples, nburn, add_var, delay_dist)
    
    
    nrow = len(fnames)
    fig = plt.figure( figsize=(12, 2*nrow) )
    
    gs = fig.add_gridspec( nrow, 2, hspace=0, wspace=0, width_ratios=[4,1] )
    ax = gs.subplots(sharex='col')

    for i in range(nrow):
        x, y, yerr = np.loadtxt( fnames[i], usecols=[0,1,2], unpack=True, delimiter=delimiter )
        x_sim, y_sim, yerr_sim = models[i]
        tau_samp = samples_chunks[i][2]
        
        if delay_dist:
            tau_rms_samp = samples_chunks[i][3]
            
        if add_var:
            if delay_dist:
                sig_med = np.median(samples_chunks[i][4])
            else:
                sig_med = np.median(samples_chunks[i][3])
                
            yerr = np.sqrt( yerr**2 + sig_med**2 )
        
        
        _, _, bars = ax[i,0].errorbar(x, y, yerr=yerr, fmt='.k', mec='k', mfc=colors[9], zorder=0) 
        [bar.set_alpha(.3) for bar in bars]
        
        ax[i,0].plot(x_sim, y_sim, color=colors[i%9], zorder=-1)
        ax[i,0].fill_between(x_sim, y_sim-yerr_sim, y_sim+yerr_sim, color=colors[i%9], alpha=.3, zorder=-2)
 
        
        if i == 0:
            ax[i,1].axis('off')
        
        if i != 0:
            
            if delay_dist:
                
                gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=ax[i,1], height_ratios=[1, 3], hspace=0, wspace=0 )
                ax1 = fig.add_subplot(gs00[0], sharex=ax[i,1])
                ax2 = fig.add_subplot(gs00[1], sharex=ax[i,1]) 
                
                tau = np.median(tau_samp)
                tau_rms = np.median(tau_rms_samp)
                
                smpls = tau_samp
                if psi_types[i-1] == 'TruncGaussian':
                    smpls = peaktomean( tau_samp, np.zeros_like(tau_samp), tau_rms_samp )
                
                frq, edges = np.histogram(smpls, bins=nbin)
                ax2.bar( edges[:-1], frq, width=np.diff(edges), align='edge', color=colors[i%9], edgecolor=colors[i%9])               
                                
                
                
                length = 10*tau_rms
                taus = np.arange( tau - 5*tau_rms, tau + 5*tau_rms, length/500 )
                if (psi_types[i-1] == 'LogGaussian') or (psi_types[i-1] == 'InverseGaussian'):
                    taus = np.linspace( -tau_rms, np.median(samples_chunks[-2][2]) + 5*tau_rms, 600 )
                    
                hi_vals = []
                lo_vals = []
                med_vals = []
                for k in range(len(taus)):
                    if psi_types[i-1] == 'Gaussian':
                        psi = Gaussian( tau_samp, tau_rms_samp, taus[k], conv=False )
                    elif psi_types[i-1] == 'Uniform':
                        psi = Uniform( tau_samp, tau_rms_samp, taus[k], conv=False )
                    elif psi_types[i-1] == 'TruncGaussian':
                        psi = TruncGaussian( tau_samp, tau_rms_samp, taus[k], np.zeros_like(tau_samp), conv=False )
                    elif psi_types[i-1] == 'LogGaussian': 
                        psi = LogGaussian( tau_samp, tau_rms_samp, taus[k], np.zeros_like(tau_samp), conv=False )
                    elif psi_types[i-1] == 'InverseGaussian':
                        psi = InverseGaussian( tau_samp, tau_rms_samp, taus[k], np.zeros_like(tau_samp), conv=False )
                        
                    hi_vals.append( np.percentile(psi, 84) )
                    lo_vals.append( np.percentile(psi, 16) )
                    med_vals.append( np.median(psi) )
                

                # norm = np.max(hi_vals)/np.max(frq)
                
                ax1.plot(taus, med_vals, color='k', zorder=1)
                ax1.fill_between(taus, lo_vals, hi_vals, color='k', alpha=.2, edgecolor='none', zorder=0) 
                ax1.set_ylim(0, np.max(med_vals)*1.5)

                ax2.axvline( np.median(tau_samp), ls='--', color='k', lw=.75 )
                
                for a in [ax1, ax2]:
                    a.set_yticklabels([])
                    a.set_yticks([])
                    
                    a.label_outer()
                    a.tick_params('both', which='major', length=6)
                    a.tick_params('both', which='minor', length=3)
                
            else:
                ax[i,1].hist(tau_samp, bins=nbin, color=colors[i%9] )
                ax2.axvline( np.median(tau_samp), ls='--', color='k', lw=.75 )
        
        
            ax[i,1].set_yticklabels([])
            ax[i,1].set_yticks([])
         
         
    tau_min = []
    tau_max = []
    for i in range(1, len(fnames)):
        tau_min.append( np.min(samples_chunks[i][2]) )
        tau_max.append( np.max(samples_chunks[i][2]) )
        
    length = np.max(tau_max) - np.min(tau_min)
    ax[-1,1].set_xlim( np.min(tau_min) - .25*length, np.max(tau_max) + .25*length ) 
         
         
         
    ax[-1,0].set_xlabel('Time [{}]'.format(time_unit))
    ax[-1,1].set_xlabel(r'$\tau$ [{}]'.format(time_unit))
    
    for a in ax.flatten():
        a.label_outer()
        a.tick_params('both', which='major', length=6)
        a.tick_params('both', which='minor', length=3)
    
    for i in range(nrow):
        if (lc_unit[i] == '') or (lc_unit[i] == 'Arbitrary Unit'):
            ytxt = 'Flux'
        elif lc_unit == 'mag':
            ytxt = 'Magnitude'
            ax[i,0].invert_yaxis()
        else:
            ytxt = 'Flux [{}]'.format(lc_unit[i])
       
       
        ax[i,0].text( -.1, .5, ytxt, transform=ax[i,0].transAxes, rotation='vertical', ha='left', va='center' )        
        ax[i,0].text( .025, .85, line_names[i], transform=ax[i,0].transAxes, ha='left', va='top' )
        
        if i != 0:
            tau_med = np.median( samples_chunks[i][2] )
            tau_hi = np.percentile( samples_chunks[i][2], 84)
            tau_lo = np.percentile( samples_chunks[i][2], 16)
            txt = r'$\tau = ' + '{:.2f}'.format(tau_med) + '^{+' + '{:.2f}'.format(tau_hi-tau_med) + '}_{-' + '{:.2f}'.format(tau_med-tau_lo) + '}$'
            ax[i,0].text( .025, .75, txt, transform=ax[i,0].transAxes, ha='left', va='top' ) 
    
           
    if output_fname is not None:
        plt.savefig(output_fname, dpi=200, bbox_inches='tight')
              
    if show:     
       plt.show()
       
    plt.cla()
    plt.clf()
    plt.close()
    
    return





def get_all_samples(samples, nburn=0):
    
    #Remove burn-in
    samples = samples[nburn:,:,:].copy()

    #Flatten
    samples = samples.reshape( (-1, samples.shape[2]) ).T
        
    return samples
    
    
def pyroa_corner_plot(samples, line_names, nburn=0, add_var=False, delay_dist=False, split=True, fname=None, show=False):

    if split:
        
        if fname is not None:
            if isinstance(fname, str) or len(fname) != len(line_names):
                raise Exception('If split=True, a file name must be given for every line')
        
        samples_chunks = get_samples_chunks(samples, nburn, add_var, delay_dist)
        
        for i in range(len(line_names)):
            samps_i = []
            
            if i == 0:
                for j in [0,1]:
                    samps_i.append( samples_chunks[i][j] ) 
                if add_var:
                    if delay_dist:
                        samps_i.append( samples_chunks[i][4] )
                    else:
                        samps_i.append( samples_chunks[i][3] )
                    
                samps_i.append( samples_chunks[-1][0] )                
                samps_i = np.array(samps_i)

                labels = []
                cols = ['A', 'B']
                if add_var:
                    cols.append( r'\sigma' )
                    
                for col in cols:
                    labels.append( r'$' + col + '_{' + line_names[0] + '}$' )     
                labels.append( r'$\Delta$' )
                    
                fig = corner.corner( samps_i.T, labels=labels, show_titles=True, plot_contours=True, smooth=1, title_fmt='.3f', label_kwargs={'fontsize':18} );
                
                if fname is not None:
                    fig.savefig(fname[i], dpi=200, bbox_inches='tight')      
                    
                if show:
                    plt.show()
                    
                plt.cla()
                plt.clf()
                plt.close()   

                
            else:
                for j in range(len(samples_chunks[i])):
                    samps_i.append( samples_chunks[i][j] )
                samps_i.append( samples_chunks[-1][0] )
                samps_i = np.array(samps_i)
                    
                labels = []
                cols = ['A', 'B', r'\tau']
                if delay_dist:
                    cols.append( r'\Delta' )
                if add_var:
                    cols.append( r'\sigma' )
                    
                for col in cols:
                    labels.append( r'$' + col + '_{' + line_names[i] + '}$' )     
                labels.append( r'$\Delta$' )
                
                
                fig = corner.corner( samps_i.T, labels=labels, show_titles=True, plot_contours=True, smooth=1, title_fmt='.3f', label_kwargs={'fontsize':18} );
                
                if fname is not None:
                    fig.savefig(fname[i], dpi=200, bbox_inches='tight')      
                    
                if show:
                    plt.show()
                    
                plt.cla()
                plt.clf()
                plt.close()   
        

    else:
        #Remove burn-in
        samples = samples[nburn:,:,:].copy()
        #Flatten
        samples = samples.reshape( (-1, samples.shape[2]) )
        
        
        labels = []
        for col in ['A', 'B']:
            labels.append( r'$' + col + '_{' + line_names[0] + '}$' )
        if add_var:
            labels.append( r'$\sigma_{' + line_names[0] + '}$' )
            
        cols = ['A', 'B', r'\tau']
        if delay_dist:
            cols.append( r'\Delta' )
        if add_var:
            cols.append( r'\sigma' )
            
        for name in line_names[1:]:
            for col in cols:
                labels.append( r'$' + col + '_{' + name + '}$' )
                
        labels.append( r'$\Delta$' )
        

        fig = corner.corner( samples, labels=labels, show_titles=True, plot_contours=True, smooth=1, title_fmt='.3f', label_kwargs={'fontsize':18} );
        
        if fname is not None:
            fig.savefig(fname, dpi=200, bbox_inches='tight')

        if show:
            plt.show()
        
        plt.cla()
        plt.clf()
        plt.close()
    
    return
