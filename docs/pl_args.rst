.. role:: python(code)
   :language: highlight
   :class:

pyPetal Arguments
==================

Required General Arguments
---------------------------

``output_dir``
    The directory used for all output.
    
    Type: str



``arg2``
    Either the list of filenames to all light curve fils, or an 
    array of the light curves themselves. If given as a list of 
    filenames, all files must be in the same directory. The first 
    line will be considered the continuum light curve.

    Type: list of str, list of float  



Optional General Arguments
----------------------------

``line_names``
    A list of the names of all lines input in ``arg2``. If :python:`None`, the 
    lines will be named in chronological order (i.e. "Line1", "Line2", etc...)

    Type: :python:`None`, list of :python:`str`

    Default: :python:`None`



``file_fmt``: The format of the light curve files input in ``arg2``. All light curve files are 
                required to be CSV in the analysis, so if :python:`file_fmt != "csv"`, it will be saved 
                in the ``light_curves/`` directory in CSV format. Currently, "csv" and "ascii" are 
                recognised. All other formats will need to be recognised by the ``astropy.table`` 
                module.
                Type: :python:`str``
                Default: "csv"

* ``verbose``: Whether or not to display text progress of the pipeline.
               Type: :python:`bool`
               Default: :python:`False`

* ``plot``: Whether or not to display plots showing the progress of the pipeline.
            Type: bool
            Default: False

* ``time_unit``: The unit to use for figures for the time axis.
                 Type: str
                 Default: "d"

* ``lc_unit``: The unit used for figures for the light curve axis. Can be a list of units or a single unit. 
               If a single unit is given, it will be assumed for all lines. pyPetal will recognize "mag" as
               as magnitude and invert the axis of all plots. All other units will be assumed to be flux units.
               Type: str, list of str
               Default: ""

* ``lag_bounds``: The range of lags to use for all pyPetal modules when searching for a lag. If ``None`` or "baseline" are 
                  input for a given line, the baseline (both positive and negative) will be used as the lag bounds. If only one
                  set of bounds is given, it will be assumed for all lines.
                  Type: list of float, ``None``, "baseline"



Module: DRW Rejection (``run_drw_rej``)
---------------------------------------

* ``nsig``: The number of $\sigma$ from the mean DRW fit to reject data points.
            Type: float
            Default: 3.0

* ``jitter``: Whether to incluse a noise ("jitter") term in the DRW fitting process.
              Type: bool
              Default: True

* ``nchain``: The number of chains for Monte Carlo sampling.
              Type: int
              Default: 10000

* ``nburn``: The number of burn-in Monte Carlo samples.
             Type: int
             Default: 3000

* ``nwalker``: The number of walkers for Monte Carlo sampling.
               Type: int
               Default: 32

* ``clip``: Celerite will use a prior for the characteristic DRW timescale $\tau_{\rm DRW}$, 
            spanning the minimum cadence to the baseline of the input light curve. If ``clip=True`` 
            for a given light curve, instead of using the minimum difference between times given for
            the light curve, it will clip these differences for values below $10^{-8}$. If one value 
            is given, it will be assumed for all light curves.
            Type: bool, list of bool 
            Default: True  

* ``reject_data``: If ``reject_data=True`` for a given light curve, it will be fit and its values will be 
                   rejected based on the value of ``nsig``. If ``reject_data=False`` for a given light curve,
                   it will not be fit to a DRW. If one value is given, it will be assumed for all light curves.
                   Type: bool, list of bool
                   Default: ``True`` for the continuum, ``False`` for all lines

* ``use_for_javelin``: If ``True``, the resulting DRW parameters ($\sigma_{\rm DRW}$, $\tau_{\rm DRW}$), will
                       be used as input to the JAVELIN module of pyPetal. The DRW parameters in each fit will be
                       fixed to the results obtained in this module.
                       Type: bool
                       Default: False



Module: Detrending (``run_detrend``)
------------------------------------

* ``K``: The number of Gaussians to use in the LinMix model.
         Type: int
         Default: 2

* ``nchain``: The number of chains to use for the Monte Carlo simulations.
              Type: int
              Default: 4

* ``miniter``: The minimum number of iterations for the Monte Carlo simulations.
               Type: int
               Default: 5000

* ``maxiter``: The maximum number of iterations for the Monte Carlo simulations.
               Type: int
               Default: 10000



Module: pyCCF (``run_pyccf``)
-----------------------------

* ``nsim``: The number of Monte Carlo simulations to run.
            Type: int
            Default: 3000

* ``interp``: The time interval with which pyCCF will interpolate the ligh curves to form the ICCF. This value must be 
              shorter than the average cadence of the ligh curves. Setting this value too low can introduce noise. If 
              set to ``None``, ``interp`` will be set to half of the average cadence of the light curves. 
              Type: float, ``None``
              Default: 2.0

* ``mcmode``: The type of resampling to perform for the Monte Carlo simulations. 0 performs both flux randomization (FR) 
              and random subset selection (RSS). 1 performs only FR. 2 performs only RSS.
              Type: int (0, 1, or 2)
              Default: 0

* ``sigmode``: The threshold for considering a measurement in the ICCF significant when computing peaks and centroids. 
               Must be within the interval (0,1). All peaks and centroids with correlation coefficient $r_{\rm max} <=$ ``sigmode`` 
               will be considered as “failed”. If set to 0, will exclude all peaks based on a p-value significance 
               test (see pyCCF documentation). 
               Type: float 
               Default: 0.2

* ``thres``: The lower limit of correlation coefficient used when calculating the centroid of the ICCF. 
             Must be within the interval (0,1). 
             Type: float
             Default: 0.8


Module: pyZDCF (``run_pyzdcf``)
-------------------------------

* ``nsim``: The number of Monte Carlo simulations to run.
            Type: int
            Default: 1000

* ``minpts``: The minimum number of points to use in each bin when computing the ZDCF. Must be larger than 11. If set 
              to 0, it will be set to 11. 
              Type: int
              Default: 0

* ``uniform_sampling``: Whether or not the light curves are uniformly sampled.
                        Type: bool
                        Default: False

* ``omit_zero_lags``: Whether or not to omit the points with zero lags when computing the ZDCF.
                      Type: bool
                      Default: True

* ``sparse``: Determines whether to use a sparse matrix implementation for reduced RAM usage. This feature is suitable 
              for longer light curves (> 3000 data points). If True, will use sparse matrix implementation. If set to "auto", 
              will use sparse matrix implementation if there are more than 3000 data points per light curve. 
              Type: bool, str
              Default: "auto"

* ``prefix``: Prefix to the output ZDCF file. 
              Type: str
              Default: "zdcf"

* ``run_plike``: Whether or not to run the PLIKE algorithm on the ZDCF to get a maximum likelihood time lag.
                 NOTE: If ``run_plike=True``, the ``plike_dir`` argument must also be specified.
                 Type: bool
                 Default: False

* ``plike_dir``: The path to the PLIKE executable.
                 Type: str, ``None``
                 Default: ``None``



Module: JAVELIN (``run_javelin``)
---------------------------------

* ``subtract_mean``: Whether or not to subtract the mean from all light curves before analysis.
                     Type: bool
                     Default: True

* ``nchain``: The number of chains to use in the MCMC.
              Type: int
              Default: 100

* ``nburn``: The number of burn-in steps to use in the MCMC.
             Type: int
             Default: 100

* ``nwalkers``: The number of walkers to use in the MCMC.
                Type: int
                Default: 100

* ``rm_type``: The type of reverberation mapping (RM) analysis to use when running JAVELIN. Can either be set 
               to "spec" for spectroscopic RM, or "phot" for photometric RM. 
               Type: str
               Default: "spec"

* ``together``: Whether or not to fit all lines to the same model. If ``together=False`` all lines will be fit
                to the continuum separately.
                Type: bool
                Default: False

* ``lagtobaseline``: A log prior is used to logarithmically penalizes lag values larger than x*baseline, where 
                     x is the value of this parameter. 
                     Type: float 
                     Default: 0.3

* ``fixed``: A list to determine what parameters to fix/vary when fitting the light curves. This should be an 
             array with a length equal to the number of parameters in the model (i.e. 2 + 3*(number of light curves) ). 
             The fitted parameters will be the two DRW parameters ( log($sigma_{\rm DRW}$), log($tau_{\rm DRW}$) ) and 
             three tophat parameters for each non-continuum light curve (lag, width, scale). Setting to 0 will fix the 
             parameter and setting to 1 will allow it to vary. If None, all parameters will be allowed to vary. The fixed 
             parameters must match the fixed value in the array input to the ``p_fix`` argument. If ``together=False``, this 
             can be input as a list of inputs, one for each line. If only one input is given, it will be assumed for each line.
             Type: ``None``, list of int
             Default: ``None``

* ``p_fix``: A list of the fixed parameters, corresponding to the elements of the fixed array. If None, all parameters will 
             be allowed to vary. Similar to ``fixed``, if ``together=False`` this can be input as a list of inputs for each line.
             If only one input is given, it will be assumed for all lines.
             Type: ``None``, list of float
             Default: ``None``

* ``output_chains``: Whether or not to output the MCMC chains to a file.
                     Type: bool
                     Default: True

* ``output_burn``: Whether or not to output the MCMC burn-in chains to a file.
                   Type: bool
                   Default: True

* ``output_logp``: Whether or not to output the MCMC log probability to a file.
                   Type: bool
                   Default: True

* ``nbin``: The number of bins to use for the output histogram plots.
            Type: int
            Default: 100



Module: Weighting (``run_weighting``)
-------------------------------------

* ``gap_size``: The minimum gap size to use to detect gaps in the continuum light curve when obtaining $N(\tau)$.
                Type: float
                Default: 20.0

* ``k``: The exponent used when calculating $P(\tau)$.
         Type: float
         Default: 2.0

* ``width``: The width of the Gaussian used to smooth the weighted distribution to find the primary peak.
             Type: float
             Default: 20.0

* ``zoom``: Whether or not to zoom in on the peak with an inset in the output plot.
            Type: bool
            Default: True

