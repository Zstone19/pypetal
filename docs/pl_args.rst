.. role:: python(code)
   :language: python
   :class: highlight



pyPetal Arguments
==================

Required General Arguments
---------------------------

.. list-table::
    :widths: 20 60 20
    :header-rows: 1

    * - Argument
      - Description
      - Type
    * - ``output_dir``
      - The directory used for all output.
      - :python:`str`
    * - ``arg2``
      - Either the list of filenames to all light curve files, or an array of the light curves themselves. If given as a list of filenames, all files must be in the same directory. The first line will be considered the continuum light curve.
      - list of :python:`str`, list of :python:`float`


.. note:: pyPetal will use the first 3 columns of the light curve files, and assume they represent the time, values, and uncertainty in the light curves.

.. warning:: By default, :python:`arg2=None`. pyPetal will raise an error if :python:`arg2=None`.



Optional General Arguments
----------------------------

.. list-table::
    :widths: 16 50 16 16
    :header-rows: 1

    * - Argument
      - Description
      - Type
      - Default
    * - ``line_names``
      - A list of the names of all lines input in ``arg2``. If :python:`None`, the lines will be named in chronological order (i.e. "Line1", "Line2", etc...)
      - :python:`None`, list of :python:`str`
      - :python:`None`
    * - ``file_fmt``
      - The format of the light curve files input in ``arg2``. All light curve files are required to be CSV in the analysis, so if :python:`file_fmt != "csv"`, it will be saved in the ``light_curves/`` directory in CSV format. Currently, "csv" and "ascii" are recognised. All other formats will need to be recognised by the ``astropy.table`` module.
      - :python:`str`
      - :python:`"csv"`
    * - ``verbose``
      - Whether or not to display text progress of the pipeline.
      - :python:`bool`
      - :python:`False`
    * - ``plot``
      - Whether or not to display plots showing the progress of the pipeline.
      - :python:`bool`
      - :python:`False`
    * - ``time_unit``
      - The unit to use for figures for the time axis.
      - :python:`str`
      - :python:`"d"`
    * - ``lc_unit``
      - The unit used for figures for the light curve axis. Can be a list of units or a single unit. If a single unit is given, it will be assumed for all lines. pyPetal will recognize "mag" as as magnitude and invert the axis of all plots. All other units will be assumed to be flux units.
      - :python:`str`, list of :python:`str`
      - :python:`""`
    * - ``lag_bounds``
      - The range of lags to use for all pyPetal modules when searching for a lag. If :python:`None` or "baseline" are input for a given line, the baseline (both positive and negative) will be used as the lag bounds. If only one set of bounds is given, it will be assumed for all lines.
      - list of :python:`None`, :python:`float`, :python:`"baseline"`
      - :python:`None`
    * - ``threads``
      - The number of threads to use for multiprocessing. This will be applied to all modules selected.
      - :python:`int`
      - 1



Module: DRW Rejection (``run_drw_rej``)
---------------------------------------

.. list-table::
    :widths: 16 50 16 16
    :header-rows: 1

    * - Argument
      - Description
      - Type
      - Default
    * - ``nsig``
      - The number of :math:`\sigma` from the mean DRW fit to reject data points.
      - :python:`float`
      - 3.0
    * - ``jitter``
      - Whether to incluse a noise ("jitter") term :math:`\sigma_n` in the DRW fitting process.
      - :python:`bool`
      - :python:`True`
    * - ``nchain``
      - The number of chains for Monte Carlo sampling.
      - :python:`int`
      - 10000
    * - ``nburn``
      - The number of burn-in Monte Carlo samples.
      - :python:`int`
      - 3000
    * - ``nwalker``
      - The number of walkers for Monte Carlo sampling.
      - :python:`int`
      - 32
    * - ``clip``
      - ``Celerite`` will use a prior for the characteristic DRW timescale :math:`\tau_{\rm DRW}`, spanning the minimum cadence to the baseline of the input light curve. If :python:`clip=True` for a given light curve, instead of using the minimum difference between times given for the light curve, it will clip these differences for values below :math:`10^{-8}`. If one value is given, it will be assumed for all light curves.
      - :python:`bool`, list of :python:`bool`
      - :python:`True`
    * - ``reject_data``
      - If :python:`reject_data=True` for a given light curve, it will be fit and its values will be rejected based on the value of ``nsig``. If :python:`reject_data=False` for a given light curve, it will not be fit to a DRW. If one value is given, it will be assumed for all light curves.
      - :python:`bool`, list of :python:`bool`
      - :python:`True` for the continuum, :python:`False` for all lines
    * - ``use_for_javelin``
      - If :python:`True`, the resulting DRW parameters :math:`(\sigma_{\rm DRW}, $\tau_{\rm DRW})`, will used as input to the JAVELIN module of pyPetal. The DRW parameters in each fit will be fixed to the results obtained in this module.
      - :python:`bool`
      - :python:`False`



Module: Detrending (``run_detrend``)
------------------------------------

.. list-table::
    :widths: 16 50 16 16
    :header-rows: 1

    * - Argument
      - Description
      - Type
      - Default
    * - ``K``
      - The number of Gaussians to use in the ``LinMix`` model.
      - :python:`int`
      - 2
    * - ``nchain``
      - The number of chains for Monte Carlo sampling.
      - :python:`int`
      - 4
    * - ``miniter``
      - The minimum number of iterations for the Monte Carlo simulations.
      - :python:`int`
      - 5000
    * - ``maxiter``
      - The maximum number of iterations for the Monte Carlo simulations.
      - :python:`int`
      - 10000


Module: pyCCF (``run_pyccf``)
-----------------------------

.. list-table::
    :widths: 16 50 16 16
    :header-rows: 1

    * - Argument
      - Description
      - Type
      - Default
    * - ``nsim``
      - The number of Monte Carlo simulations to run.
      - :python:`int`
      - 3000
    * - ``interp``
      - The time interval with which pyCCF will interpolate the ligh curves to form the ICCF. This value must be shorter than the average cadence of the ligh curves. Setting this value too low can introduce noise. If set to :python:`None`, ``interp`` will be set to half of the average cadence of the light curves.
      - :python:`float`, :python:`None`
      - 2.0
    * - ``mcmode``
      - The type of resampling to perform for the Monte Carlo simulations. 0 performs both flux randomization (FR) and random subset selection (RSS). 1 performs only FR. 2 performs only RSS.
      - :python:`int`
      - 0
    * - ``sigmode``
      - The threshold for considering a measurement in the ICCF significant when computing peaks and centroids. Must be within the interval (0,1). All peaks and centroids with correlation coefficient :math:`r_{\rm max} \leq` ``sigmode`` will be considered as “failed”. If set to 0, will exclude all peaks based on a p-value significance test (see pyCCF documentation).
      - :python:`float`
      - 0.2
    * - ``thres``
      - The lower limit of correlation coefficient used when calculating the centroid of the ICCF. Must be within the interval (0,1).
      - :python:`float`
      - 0.8


Module: pyZDCF (``run_pyzdcf``)
-------------------------------

.. list-table::
    :widths: 16 50 16 16
    :header-rows: 1

    * - Argument
      - Description
      - Type
      - Default
    * - ``nsim``
      - The number of Monte Carlo simulations to run.
      - :python:`int`
      - 1000
    * - ``minpts``
      - The minimum number of points to use in each bin when computing the ZDCF. Must be larger than 11. If set to 0, it will be set to 11.
      - :python:`int`
      - 0
    * - ``uniform_sampling``
      - Whether or not the light curves are uniformly sampled.
      - :python:`bool`
      - :python:`False`
    * - ``omit_zero_lags``
      - Whether or not to omit the points with zero lags when computing the ZDCF.
      - :python:`bool`
      - :python:`True`
    * - ``sparse``
      - Determines whether to use a sparse matrix implementation for reduced RAM usage. This feature is suitable for longer light curves (> 3000 data points). If True, will use sparse matrix implementation. If set to "auto", will use sparse matrix implementation if there are more than 3000 data points per light curve.
      - :python:`bool`, :python:`str`
      - :python:`"auto"`
    * - ``prefix``
      - Prefix to the output ZDCF file.
      - :python:`str`
      - :python:`"zdcf"`
    * - ``run_plike``
      - Whether or not to run the PLIKE algorithm on the ZDCF to get a maximum likelihood time lag. **NOTE**: If :python:`run_plike=True`, the ``plike_dir`` argument must also be specified.
      - :python:`bool`
      - :python:`False`
    * - ``plike_dir``
      - The path to the PLIKE executable.
      - :python:`str`, :python:`None`
      - :python:`None`



Module: pyROA (``run_pyroa``)
-------------------------------

.. list-table::
    :widths: 16 50 16 16
    :header-rows: 1

    * - Argument
      - Description
      - Type
      - Default
    * - ``nchain``
      - The number of chains for Monte Carlo sampling.
      - :python:`int`
      - 20000
    * - ``nburn``
      - The number of burn-in steps to remove from the Monte Carlo samples.
      - :python:`int`
      - 15000
    * - ``together``
      - Whether or not to fit the time lags of all light curves together.
      - :python:`bool`
      - :python:`True`
    * - ``init_tau``
      - The initial guess for the time lag. If one value is given, it will be used for all lines. If :python:`None`, it will be set to 10. for each line.
      - :python:`float`, list of :python:`float`, :python:`None`
      - :python:`None`
    * - ``subtract_mean``
      - Whether or not to subtract the mean from all light curves before analysis.
      - :python:`bool`
      - :python:`True`
    * - ``div_mean``
      - Whether or not to divide the light curves by their mean before analysis. This will occur before the mean is subtracted if :python:`subtract_mean=True`.
      - :python:`bool`
      - :python:`False`
    * - ``add_var``
      - Whether or not to add additional uncertainty in the data, same as the PyROA argument. If :python:`together=False`, multiple values may be given for each line. If only one value is given, it will be assumed for all lines.
      - :python:`bool`, list of :python:`bool`
      - :python:`True`
    * - ``delay_dist``
      - Same as the ``delay_dist`` argument for PyROA. If :python:`together=False`, multiple values may be given for each line. If only one value is given, it will be assumed for all lines.
      - :python:`bool`, list of :python:`bool`
      - :python:`True`
    * - ``psi_types``
      - Same as the ``psi_types`` argument for PyROA - the form of the delay distribution. This will only affect the output if :python:`delay_dist=True`. If :python:`together=False`, multiple values may be given for each line. If only one value is given, it will be assumed for all lines. If :python:`None`, a Gaussian distribution will be assumed.
      - :python:`str`, list of :python:`str`, :python:`None`
      - :python:`None`
    * - ``objname``
      - The name of the object to use for PyROA analysis. This will apply to the output file names and figures. If :python:`None`, this will be set to "pyroa".
      - :python:`str`, :python:`None`
      - :python:`None`



Module: JAVELIN (``pypetal_jav.run_pipeline``)
----------------------------------------------

.. list-table::
    :widths: 16 50 16 16
    :header-rows: 1

    * - Argument
      - Description
      - Type
      - Default
    * - ``subtract_mean``
      - Whether or not to subtract the mean from all light curves before analysis.
      - :python:`bool`
      - :python:`True`
    * - ``nchain``
      - The number of chains to use in the MCMC.
      - :python:`int`
      - 100
    * - ``nburn``
      - The number of burn-in steps to use in the MCMC.
      - :python:`int`
      - 100
    * - ``nwalkers``
      - The number of walkers to use in the MCMC.
      - :python:`int`
      - 100
    * - ``rm_type``
      - The type of reverberation mapping (RM) analysis to use when running JAVELIN. Can either be set to "spec" for spectroscopic RM, or "phot" for photometric RM.
      - :python:`str`
      - :python:`"spec"`
    * - ``together``
      - Whether or not to fit all lines to the same model. If :python:`together=False` all lines will be fit to the continuum separately.
      - :python:`bool`
      - :python:`False`
    * - ``lagtobaseline``
      - A log prior is used to logarithmically penalizes lag values larger than `x`*baseline, where `x` is the value of this parameter.
      - :python:`float`
      - 0.3
    * - ``fixed``
      - A list to determine what parameters to fix/vary when fitting the light curves. This should be an array with a length equal to the number of parameters in the model. The fitted parameters will be the two DRW parameters :math:`( \log(\sigma_{\rm DRW}), \log(\tau_{\rm DRW}) )` and (3 or 4) tophat parameters for each non-continuum light curve. Setting to 0 will fix the parameter and setting to 1 will allow it to vary. If None, all parameters will be allowed to vary. The fixed parameters must match the fixed value in the array input to the ``p_fix`` argument. If :python:`together=False`, this can be input as a list of inputs, one for each line. If only one input is given, it will be assumed for each line.
      - :python:`None`, list of :python:`int`
      - :python:`None`
    * - ``p_fix``
      - A list of the fixed parameters, corresponding to the elements of the fixed array. If :python:`None`, all parameters will be allowed to vary. Similar to ``fixed``, if :python:`together=False` this can be input as a list of inputs for each line. If only one input is given, it will be assumed for all lines.
      - :python:`None`, list of :python:`float`
      - :python:`None`
    * - ``output_chains``
      - Whether or not to output the MCMC chains to a file.
      - :python:`bool`
      - :python:`True`
    * - ``output_burn``
      - Whether or not to output the MCMC burn-in chains to a file.
      - :python:`bool`
      - :python:`True`
    * - ``output_logp``
      - Whether or not to output the MCMC log probability to a file.
      - :python:`bool`
      - :python:`True`
    * - ``nbin``
      - The number of bins to use for the output histogram plots.
      - :python:`int`
      - 100


Determining the number of parameters in the JAVELIN model:

.. list-table::
    :widths: 20 20 25 35
    :header-rows: 1

    * - :python:`rm_type`
      - :python:`together`
      - Number of Parameters
      - Parameter Names
    * - :python:`"spec"`
      - :python:`True`
      - :math:`2 + 3 \cdot ({\rm number of light curves})`
      - :math:`\log(\sigma_{\rm DRW})`, :math:`\log(\tau_{\rm DRW})`, :math:`t_1`, :math:`w_1`, :math:`s_1`, :math:`t_2`, ...
    * - :python:`"spec"`
      - :python:`False`
      - 5 per line
      - :math:`\log(\sigma_{\rm DRW})`, :math:`\log(\tau_{\rm DRW})`, :math:`t`, :math:`w`, :math:`s`
    * - :python:`"phot"`
      - :python:`True`
      - 6 per line
      - :math:`\log(\sigma_{\rm DRW})`, :math:`\log(\tau_{\rm DRW})`, :math:`t`, :math:`w`, :math:`s`, :math:`\alpha`


.. note:: If :python:`use_for_javelin=True` in the DRW Rejection module, and ``fixed/p_fix`` are set in the JAVELIN module, the DRW fitting results will be used instead of the input fixed parameter values.

.. note:: If :python:`rm_type="phot"`, only one light curve can be modeled to a given continuum. Therefore, pyPetal will set :python:`together=False`.



Module: Weighting (``pypetal.run_weighting``)
---------------------------------------------

.. list-table::
    :widths: 16 50 16 16
    :header-rows: 1

    * - Argument
      - Description
      - Type
      - Default
    * - ``gap_size``
      - The minimum gap size to use to detect gaps in the continuum light curve when obtaining :math:`N(\tau)`.
      - :python:`float`
      - 20.0
    * - ``k``
      - The exponent used when calculating :math:`P(\tau)`.
      - :python:`float`
      - 2.0
    * - ``width``
      - The width of the Gaussian used to smooth the weighted distribution to find the primary peak.
      - :python:`float`
      - 20.0
    * - ``rel_height``
      - The relative height (0-1) to use for the peak-finding algorithm.
      - :python:`float`
      - 0.99
    * - ``zoom``
      - Whether or not to zoom in on the peak with an inset in the output plot.
      - :python:`bool`
      - :python:`True`
