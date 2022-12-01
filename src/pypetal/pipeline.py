import pypetal.modules as modules
from pypetal.formatting import write_data

import os

import numpy as np
from astropy.table import Table

    
def run_pipeline(output_dir, arg2, 
                 line_names=None, 
                 run_drw_rej=False, drw_rej_params={},
                 run_pyccf=False, pyccf_params={},
                 run_pyzdcf=False, pyzdcf_params={},
                 run_javelin=False, javelin_params={},
                 **kwargs):
    
    """ Run the pyPetal pipeline on a list of files. Individual modules can be specified, but are not run by default.
    
    
    Default parameters
    -------------------    
        
    output_dir : str
        Path to the directory containing the output.
            
    arg2 : list
        Either a list of paths to the light curve files, or an array of the light curves.
        
    .. note:: The first light curve will be assumed to be the continuum light curve.
    .. note:: pyPetal will use the first 3 columns of the light curve files, and assume they represent the time, values, and uncertainty in the light curves.
    .. warning:: By default, set to ``None``. If left to ``None``, pyPetal will raise an error.    
        
    line_names : str, list of str, optional
        Name(s) of the line(s). If ``None``, the lines will be named in chronological order (i.e. Line 1, Line 2, etc).
        Default is ``None``.
    
    
    
    DRW rejection parameters
    -------------------------
    
    nsig : float
        The algorithm will flag all data points greater than nsig*sigma away from the mean of the DRW fit to the light curve. 
        Default is 1.

    jitter : bool
        If ``True``, will include a noise term in the model to account for extra noise. Default is ``True``.

    nwalker : int
        Number of walkers for the MCMC. Default is 32.

    nburn : int
        Number of burn-in steps for the MCMC. Default is 300.
        
    nchain : int
        Number of steps for the MCMC. Default is 1000.

    reject_data : bool, list of bool
        If ``True``, will use light curves without the rejected data points for the rest of the pipeline. These
        light curves will be saved in the given output directory as csv files labeled '{line name}_data.dat'. This
        parameter can be input as a list, to specify for each light curve. If one value is given and multiple light
        curves are input, the same value will be used for all light curves. Default is ``True``.
        
    clip : bool, list of bool
        If ``True``, the light curves will be clipped so that the time between two observations is not less than $10^{-8}$ days.
        This can also be input as a list of bool for each light curve. If one value is given and multiple light
        curves are input, the same value will be used for all light curves. Default is ``True``.
        
    use_for_javelin : bool
        If ``True``, the output DRW parameters will be used as fixed DRW parameters for JAVELIN. The DRW parameters for the continuum
        will be used as these fixed parameters in JAVELIN. Default is ``False``.
             
    
    PyCCF parameters    
    -----------------
    
    lag_bounds : (2,) list
        Lower and upper bounds to search for the lag. If None, will be set to [-200, 200].
        Default is ``None``. 
    
    interp : float
        The interval with which pyCCF will interpolate the ligh curves to form the ICCF. This value must be 
        shorter than the average cadence of the ligh curves. Setting this value too low can introduce noise.
        If set to ``None``, will be set to half of the average cadence of the light curves. The default is ``None``. 
        
    nsim : int
        The number of Monte Carlo simulations to run. The default is 1000.
        
    mcmode: int
        The type of resampling to do for the Monte Carlo Simulations. 0 performs both FR and RSS, 1 performs FR, and 2 performs RSS.
        The default is 0.
        
    sigmode : float
        The threshold for considering a measurement in the ICCF significant when computing peaks and centroids. Must be within the 
        interval (0,1). All peaks and centroids with correlation coefficient r_max <= sigmode will be considered as "failed". 
        If set to 0, will exclude all peaks based on a p-value significance test (see pyCCF documentation). The default is 0.2. 
        
    thres : float
        The lower limit of correlation coefficient used when calculating the centroid of the ICCF. Must be within 
        the interval (0,1). The default is 0.8.
    
        
    PyZDCF parameters
    -----------------
    
    nsim : int
        Number of Monte Carlo simulations to run. Default is 500.
        
    minpts : int
        The minimum number of points to use in each bin when computing the ZDCF. 
        Must be larger than 11. If set to 0, it will be set to 11. Default is 0
        
    uniform_sampling : bool
        If ``True``, the light curves will be assumed to be uniformly sampled.
        Default is ``False``.
        
    omit_zero_lags : bool
        If ``True``, will omit the points with zero lags when computing the ZDCF.
        Default is ``True``.
        
    sparse : bool, str
        Determines whether to use a sparse matrix implementation for reduced RAM usage.
        This feature is suitable for longer light curves (> 3000 data points). If True, will
        use sparse matrix implementation. If set to 'auto', will use sparse matrix implementation
        if there are more than 3000 data points per light curve. Default is 'auto'.
        
    prefix : str
        Prefix to the output ZDCF file. Default is 'zdcf'.

    run_plike : bool
        Whether to run the PLIKE algorithm on the ZDCF to get a maximum likelihood time lag.
        Default is False. Note: the directory containing the PLIKE executable must be input
        with the ``plike_dir`` argument, and the range of lags to search must be input with the
        ``lag_bounds`` argument.
        
    plike_dir : str
        Path to the PLIKE executable. Default is ``None``.

    lag_bounds : (2,) list
        The lower and upper bounds of lags to search for the PLIKE algorithm. If ``None``, the lower and 
        upper bounds will be the (negative, positive) baseline of the light curves. Default is ``None``.
            
        
    JAVELIN parameters
    ------------------

    rm_type : str
        The type of analysis to use when running JAVELIN. Can either be set to 'spec' for spectroscopic RM, or 
        'phot' for photometric RM. Default is 'spec'.

    lagtobaseline : float
        A log  prior is used to logarithmically penalizes lag values larger than 
        x*baseline, where x is the value of this parameter. Default is 0.3.
    
    laglimit : (2,) list, str
        The upper and lower bounds of the region of lags to search. If 'baseline', the it will be 
        set to [-baseline, basleline]. Default is 'baseline'.

    fixed : list
        A list to determine what parameters to fix/vary when fitting the light curves. This should be an array
        with a length equal to the number of parameters in the model (i.e. 2 + 3*(number of light curves) ). The fitted 
        parameters will be the two DRW parameters ( log(sigma), log(tau) ) and three tophat parameters for each
        non-continuum light curve (lag, width, scale). Setting to 0 will fix the parameter and setting to 1 will allow it to vary.
        If ``None``, all parameters will be allowed to vary. The fixed parameters must match the fixed value in the 
        array input to the ``p_fix`` argument. default is ``None``. 
    
    p_fix : list
        A list of the fixed parameters, corresponding to the elements of the ``fixed`` array.
        If ``None``, all parameters will be allowed to vary. Default is ``None``.
    
    subtract_mean : bool
        If ``True``, will subtract the mean of all light curves before analysis. Default is ``True``.
            
    nwalkers : int
        The number of walkers to use in the MCMC. Default is 100.
    
    nchain : int
        The number of steps to take in the MCMC. Default is 100.
        
    nburn : int 
        The number of steps to take in the burn-in phase of the MCMC. Default is 100.
    
    output_chains : bool
        If ``True``, will output the MCMC chains to a file. Default is ``True``.
        
    output_burn : bool
        If ``True``, will output the MCMC burn-in chains to a file. Default is ``True``.
        
    output_logp : bool
        If ``True``, will output the MCMC log probability to a file. Default is ``True``.            
    
    .. warning:: JAVELIN cannot utilize multiple bands for photometric RM, so ``together" must be set to ``False" if ``rm_type='phot'``.

    
    Additional kwargs
    -----------------
    
    verbose : bool
        Whether to print out progress. Default is ``False``.
        
    use_weights: bool
        Whether to use weights (from Grier et al. (2017)) in calculating the time lags. Default is ``False``. 
        
    file_fmt : str
        The format of the input files. All files must be 'csv' to be used. If the files are not 'csv',
        the files will be converted to 'csv' and used for subsequent analysis. Default is 'csv'.
        
    time_unit : str
        The unit of the time values. Default is 'd'.
        
    lc_unit : str
        The unit of the light curve values and uncertainties. Default is 'mag'.
            
            
    Returns
    -------
        
    res : dict
        A dict of dicts for the output of each module.
            
    Outputs for each module:
    
    DRW rejection: 'drw_rej_res'
    
        * masks : list of bool
            A list of the masks for each light curve. ``True`` indicates that the light curve value was rejected.

        * reject_data : bool 
            ``True`` if the rejected points were not used for the rest of the pipeline.

        * taus : list of floats
            List of DRW tau values from the MCMC chains.

        * sigmas : list of floats
            List of DRW sigma values from the MCMC chains.

        * jitters : list of floats
            List of jitter values from the MCMC chains. If ``jitter=False``, this is ``None``.



    pyCCF: 'pyccf_res'
    
        Each line will have a different dict with the following outputs
        
        * CCF : list of floats
            The cross-correlation function.
            
        * CCF_lags : list of floats
            The lags corresponding to the CCF.
        
        * centroid : float
            The median of the CCCD.
        
        * centroid_err_lo : float
            The lower error on the centroid.
        
        * centroid_err_hi : float
            The upper error on the centroid.
        
        * peak : float
            The median of the CCPD.
        
        * peak_err_hi : float
            The upper error on the peak.
        
        * peak_err_lo : float
            The lower error on the peak.
        
        * CCCD_lags : list of floats
            The lags corresponding to the CCCD.
        
        * CCPD_lags : list of floats
            The lags corresponding to the CCPD.
        
        * name : str
            The name of the light curve.
            
            
        If ``use_weights=True``, the following outputs will be added:
        
        * weighted_CCCD : list of floats
            The weighted CCCD.
            
        * weighted_CCPD : list of floats
            The weighted CCPD.
            
    .. note:: The weighted distributions are ouptut as the samples that would be required to create the weighted histogram. In simpler terms, the original sample is downsampled in each bin to match that of the weighted distribution.



    pyZDCF: 'pyzdcf_res'
    
        Each line will have a different pandas DataFrame with the following columns:
        
        * tau : list of floats
            The lags.
            
        * -sig(tau) : list of floats
            The lower error on the lags.
            
        * +sig(tau) : list of floats
            The upper error on the lags.
            
        * dcf : list of floats
            The ZDCF.
            
        * -err(dcf) : list of floats
            The lower error on the ZDCF.
            
        * +err(dcf) : list of floats
            The upper error on the ZDCF.
            
            
            
    PLIKE: 'plike_res'
    
        If ``run_plike=True``, each light curve will have a dict with the following outputs:
        
        * output : astropy.table.Table 
            The output of the PLIKE algorithm. This contains the lags, likelihoods, and ZDCF (and its errors).
        
        * ML_lag : float
            The maximum likelihood lag.
        
        * ML_lag_err_lo : float
            The lower error on the maximum likelihood lag.
        
        * ML_lag_err_hi : float
            The upper error on the maximum likelihood lag.
            
            
            
    JAVELIN: 'javelin_res'
    
        * cont_hpd: list of floats
            The HPD (highest posterior density) interval for the continuum DRW parameters. If both DRW parameters are fixed, this 
            will be ``None``.
            
        * tau : list of floats
            The DRW tau values from the MCMC chains.
        
        * sigma : list of floats
            The DRW sigma values from the MCMC chains.
        
        * tophat_params : list of floats
            The tophat parameters from the MCMC chains. The order is [lag, width, scale] for each light curve.    
        
        * hpd : list of floats
            The HPD (highest posterior density) interval for the combined model.
        
        * cont_model : javelin.lcmodel.Cont_Model
            The continuum model. If both DRW parameters are fixed, this will be ``None``.
        
        * rmap_model : javelin.lcmodel.Rmap_Model, javelin.lcmodel.Pmap_Model
            The rmap model.
        
        * cont_dat : javelin.zylc.LightCurve
            The continuum light curve.
        
        * tot_dat : javelin.zylc.LightCurve
            The ``LightCurve`` object storing all light curves.
        
        * bestfit_model : javelin.lcmodel.Rmap_Model
            The combined model after using the best-fit DRW and tophat parameters to fit each light curve.
                
    """

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
        file_fmt = 'csv'
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
    
    
    #Read in general kwargs
    general_kwargs = kwargs
    default_kwargs = {
        'verbose': False,
        'time_unit': 'd',
        'lc_unit': 'Arbitrary Units',
        'file_fmt': 'csv'
    }
    general_kwargs = { **default_kwargs, **general_kwargs }
    
    verbose = general_kwargs['verbose']
    time_unit = general_kwargs['time_unit']
    lc_unit = general_kwargs['lc_unit']
    file_fmt = general_kwargs['file_fmt']
    
    
    
    #Read in DRW rejection kwargs
    default_kwargs = {
        'jitter' : True,
        'nsig' : 1,
        'nwalkers': 100,
        'nburn': 300,
        'nchain': 1000,
        'clip': np.full( len(fnames), True),
        'reject_data': np.full( len(fnames), True),
        'use_for_javelin': False
    }
    drw_kwargs = { **default_kwargs, **drw_rej_params }
    
    if isinstance(drw_kwargs['reject_data'], bool):
        drw_kwargs['reject_data'] = np.full( len(fnames), drw_kwargs['reject_data'] )
    
    
    
    
    
    #Name lines if unnamed
    if line_names is None:
        line_names = np.zeros( len(line_fnames), dtype=str )
        for i in range(len(line_fnames)):
            line_names.append('Line {}'.format(i+1))


    #Get directory common to all files
    input_dir = os.path.dirname( os.path.realpath(cont_fname) )
    for i in range(len(line_fnames)):
        assert input_dir == os.path.dirname( os.path.realpath(line_fnames[i]) )        
    
    
    if 'together' in javelin_params:
        javelin_together = javelin_params['together']
    else:
        javelin_together = False
        
        
    if 'rm_type' in javelin_params:
        rm_type = javelin_params['rm_type']
    else:
        rm_type = 'spec'
        
    if (rm_type == 'phot') & (javelin_together == True):
        print('ERROR: JAVELIN cannot do phtotometric RM with more than two lines.')
        print('Setting together=False')
        javelin_together = False
        javelin_params['together'] = False
        
    
    
    #Create subdirectories for each line and javelin
    for i in range(len(fnames)):
        os.makedirs( output_dir + line_names[i], exist_ok=True )
        
        if run_drw_rej:
            os.makedirs( output_dir + line_names[i] + '/drw_rej', exist_ok=True )
        
    for i in range(len(line_fnames)):
        if run_pyccf:
            os.makedirs( output_dir + line_names[i+1] + '/pyccf', exist_ok=True )
        if run_pyzdcf:
            os.makedirs( output_dir + line_names[i+1] + '/pyzdcf', exist_ok=True )
        
        if run_javelin:
            if javelin_together:
                os.makedirs( output_dir + 'javelin/', exist_ok=True )
            else:
                os.makedirs( output_dir + line_names[i+1] + '/javelin', exist_ok=True )
                
    #Make subdirectories for light curves
    os.makedirs( output_dir + 'light_curves/', exist_ok=True )
        
        
    if file_fmt != 'csv':
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
    pyccf_res = {}
    pyzdcf_res = {}
    javelin_res = {}

    if run_drw_rej:
        
        if np.any( drw_kwargs['reject_data'] ):
            os.makedirs( output_dir + 'rejected_lcs/', exist_ok=True )
        
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
        if np.any( drw_kwargs['reject_data'] ):
            cont_fname = output_dir + 'rejected_lcs/' + line_names[0] + '_data.dat'
            
            for i in range(len(line_fnames)):
                line_fnames[i] = output_dir + 'rejected_lcs/' + line_names[i] + '_data.dat'            

        
        if 'use_for_javelin' in drw_rej_params:
            
            if drw_rej_params['use_for_javelin']:
                
                if javelin_together:                
                    tau_cont = drw_rej_res['taus'][0]
                    sig_cont = drw_rej_res['sigmas'][0]
                
                
                    if 'fixed' in javelin_params:
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
                    taus = []
                    sigmas = []

                    tau_cont = drw_rej_res['taus'][0]
                    sig_cont = drw_rej_res['sigmas'][0]
                    
                        
                    if 'fixed' in javelin_params:
                        for i in range(len(line_fnames)):
                            javelin_params['fixed'][i][0] = 0
                            javelin_params['fixed'][i][1] = 0
                            
                            javelin_params['p_fix'][i][0] = np.log(tau_cont)
                            javelin_params['p_fix'][i][1] = np.log(sig_cont)
                            
                    else:
                        javelin_params['fixed'] = np.ones( ( len(line_fnames), 5 ) )
                        javelin_params['p_fix'] = np.zeros( ( len(line_fnames), 5 ) )
                        
                        for i in range( len(line_fnames) ):
                            javelin_params['fixed'][i][0] = 0
                            javelin_params['fixed'][i][1] = 0
                            
                            javelin_params['p_fix'][i][0] = np.log(sigmas_cont)
                            javelin_params['p_fix'][i][1] = np.log(taus_cont)
                            
    else:
        
        #Output light curve files with masks
        for i in range(len(fnames)):
            
            if i == 0:
                x, y, yerr = np.loadtxt( cont_fname, delimiter=',', unpack=True, usecols=[0,1,2] )
            else:
                x, y, yerr = np.loadtxt( line_fnames[i-1], delimiter=',', unpack=True, usecols=[0,1,2] )
            
            write_data( [x,y,yerr,np.full_like(x, False, dtype=bool)], output_dir + 'light_curves/' + line_names[i] + '.dat' )
            


    if run_pyccf:
        pyccf_res = modules.pyccf_tot(cont_fname, line_fnames, line_names, output_dir, general_kwargs, pyccf_params)        
        
    if run_pyzdcf:
        pyzdcf_res, plike_res = modules.pyzdcf_tot(cont_fname, line_fnames, line_names, output_dir, general_kwargs, pyzdcf_params)
        
    if run_javelin:
        javelin_res = modules.javelin_tot(cont_fname, line_fnames, line_names, output_dir, general_kwargs, javelin_params)
    

    if file_fmt != 'csv':
        import shutil
        shutil.rmtree(input_dir)


    #Compile all results into a single dict    
    tot_res = {}
    
    if run_drw_rej:
        tot_res['drw_rej_res'] = drw_rej_res
    
    if run_javelin:
        tot_res['javelin_res'] = javelin_res
    
    if run_pyccf:
        tot_res['pyccf_res'] = pyccf_res
        
    if run_pyzdcf:
        tot_res['pyzdcf_res'] = pyzdcf_res
        tot_res['plike_res'] = plike_res


    if not isinstance(arg2[0], str):
        import shutil
        shutil.rmtree( output_dir + 'input_lcs/' ) 

    
    return tot_res