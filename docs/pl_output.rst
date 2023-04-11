pyPetal Output
===============

Each module run in pyPetal has its own output, most of which are dictionaries with a variety of keys. These keys lead to the different output data from each module. Here we provide an in-depth description of all output and how to access them.
Firstly, the output from the pipeline itself (i.e. from ``pyPetal.pipeline.run_pipeline``) will be a dictionary, containing the output dictionaries for each module run in pyPetal (except the detrending module) - DRW-based outlier rejection, pyCCF, pyZDCF, PyROA.
The key referencing each output will be labeled ``(name)_res``, where ``(name)`` is the name of the module, in addition to the output from PLIKE (if it is run). Therefore, the possible keys for the pyPetal output are:

* ``drw_rej_res``
* ``pyccf_res``
* ``pyzdcf_res``
* ``plike_res``

In addition, ``pypetal_jav.pipeline.run_pipeline`` will output a dictionary of results. The weighting module ``pypetal.pipeline.run_weighting`` will also output a dictionary.


Module: DRW Rejection
---------------------

The output dictionary from this module differs from the other modules. Each value (except ``reject_data``) in the dictionary is a list of outputs, one for each line in which :python:`reject_data=True`. These will be in the same order as the input light curves.

.. list-table::
    :widths: 20 60 20
    :header-rows: 1

    * - Key
      - Description
      - Type
    * - ``masks``
      - The DRW rejection mask, where :python:`True` indicates a given data point is rejected. These will be the same length of the input light curve.
      - list of :python:`bool`
    * - ``reject_data``
      - A copy of the input ``reject_data`` argument.
      - list of :python:`bool`
    * - ``taus``
      - A list of the MCMC samples for :math:`\tau_{\rm DRW}`.
      - list of :python:`float`
    * - ``sigmas``
      - A list of the MCMC samples for :math:`\sigma_{\rm DRW}`.
      - list of :python:`float`
    * - ``jitters``
      - A list of the MCMC samples for the jitter term :math:`\sigma_n`. If the argument :python:`jitter=False`, this will not be present.
      - list of :python:`float`


Module: PyCCF
-------------

The output pyCCF module will be a list of dictionaries, one for each line, with the line specified in each dictionary. These will be in the same order as the input light curves.

.. list-table::
    :widths: 20 60 20
    :header-rows: 1

    * - Key
      - Description
      - Type
    * - ``CCF``
      - The output cross-correlation function.
      - list of :python:`float`
    * - ``CCF_lags``
      - The lags corresponding to the CCF.
      - list of :python:`float`
    * - ``centroid``
      - The median of the CCCD.
      - :python:`float`
    * - ``centroid_err_lo``
      - The lower error on the centroid.
      - :python:`float`
    * - ``centroid_err_hi``
      - The upper error on the centroid.
      - :python:`float`
    * - ``peak``
      - The median of the CCPD.
      - :python:`float`
    * - ``peak_err_lo``
      - The lower error on the peak.
      - :python:`float`
    * - ``peak_err_hi``
      - The upper error on the peak.
      - :python:`float`
    * - ``CCCD_lags``
      - The lags corresponding to the CCCD.
      - list of :python:`float`
    * - ``CCPD_lags``
      - The lags corresponding to the CCPD.
      - list of :python:`float`
    * - ``name``
      - The name of the line.
      - :python:`str`



Module: pyZDCF
--------------

The pyZDCF can have one or two outputs, depending on the value of ``run_plike``. If ``run_plike=True``, there will be a pyZDCF output and a PLIKE output in the output dictionary.

The pyZDCF output will be a list of ``pandas.DataFrame`` objects, which are output from pyZDCF itself. These will be in the same order as the input light curves. These ``DataFrame`` objects have the following columns:

.. list-table::
    :widths: 20 60 20
    :header-rows: 1

    * - Column
      - Description
      - Type
    * - ``tau``
      - The time lag.
      - :python:`float`
    * - ``-sig(tau)``
      - The lower error on the time lag.
      - :python:`float`
    * - ``+sig(tau)``
      - The upper error on the time lag.
      - :python:`float`
    * - ``dcf``
      - The ZDCF value at that lag.
      - :python:`float`
    * - ``-err(dcf)``
      - The lower error on the ZDCF value.
      - :python:`float`
    * - ``+err(dcf)``
      - The upper error on the ZDCF value.
      - :python:`float`
    * - ``#bin``
      - The number of points in the given :math:`tau` bin.
      - :python:`int`


The PLIKE output will be a list of dictionaries, one for each line. Each dictionary will contain an ``astropy.table.Table`` object under the ``output``, which contain the output from PLIKE, read from the output file. Each table will have the following columns:

.. list-table::
    :widths: 20 60 20
    :header-rows: 1

    * - Column
      - Description
      - Type
    * - ``lag``
      - The time lag.
      - :python:`float`
    * - ``r``
      - The ZDCF value at that lag.
      - :python:`float`
    * - ``-dr``
      - The lower error on the ZDCF.
      - :python:`float`
    * - ``+dr``
      - The upper error on the ZDCF.
      - :python:`float`
    * - ``likelihood``
      - The likelihood value at that lag.
      - :python:`float`


Each dictionary will have the following keys:

.. list-table::
    :widths: 20 60 20
    :header-rows: 1

    * - Key
      - Description
      - Type
    * - ``output``
      - The output from PLIKE.
      - ``astropy.table.Table``
    * - ``ML_lag``
      - The maximum likelihood lag.
      - :python:`float`
    * - ``ML_lag_err_lo``
      - The lower error on the maximum likelihood lag.
      - :python:`float`
    * - ``ML_lag_err_hi``
      - The upper error on the maximum likelihood lag.
      - :python:`float`



Module: PyROA
-------------

Normally, the PyROA code outputs a ``PyROA.PyROA.Fit`` object which contain its results. However, due to threading issues, this cannot be done in pyPetal. To circumvent this, pyPetal has a similar class called ``MyFit``, which contains some of the attributes of the original ``Fit`` class (see below).
In general, the PyROA module's output will be dependent on the input value of the ``together`` argument. If :python:`together=False`, there will be a list of ``MyFit`` objects, one for each line. If :python:`together=True`, there will only be one ``MyFit`` object.

The ``MyFit`` object will have the following attributes:

.. list-table::
    :widths: 20 60 20
    :header-rows: 1

    * - Attribute
      - Description
      - Type
    * - ``samples``
      - The MCMC samples.
      - list of :python:`float`
    * - ``samples_flat``
      - The flattened MCMC samples.
      - list of :python:`float`
    * - ``models``
      - The ROA model fit (times, values, and errors) for each light curve.
      - list of :python:`float`
    * - ``t``
      - The time for the driving light curve model fit.
      - list of :python:`float`
    * - ``X``
      - The driving light curve model fit.
      - list of :python:`float`
    * - ``X_errs``
      - The error in the driving light curve model fit.
      - list of :python:`float`



Module: JAVELIN
---------------

The JAVELIN module's output will have a different structure depending on the value of ``together``. If :python:`together=False`, there will be a list of dictionaries for each line, in the order of the light curves given. If :python:`together=True`, there will only be one output dictionary. However, in both cases, the keys will be the same.

The output dictionary(ies) will have the following keys:

.. list-table::
    :widths: 20 60 20
    :header-rows: 1

    * - Key
      - Description
      - Type
    * - ``cont_hpd``
      - The highest posterior density (HPD) interval for the initial continuum fit. If both DRW parameters are fixed, this will be None. The first column corresponds to :math:`\sigma_{\rm DRW}`, and the second corresponds to :math:`\tau_{\rm DRW}`.
      - list of :python:`float`, :python:`None`
    * - ``tau``
      - The list of MCMC samples for :math:`\tau_{\rm DRW}`.
      - list of :python:`float`
    * - ``sigma``
      - The list of MCMC samples for :math:`\sigma_{\rm DRW}`.
      - list of :python:`float`
    * - ``tophat_params``
      - The list of MCMC samples for the tophat parameters. These tophat parameters will be ordered in the same way as the input light curves.
      - list of :python:`float`
    * - ``hpd``
      - The HPD interval for the combined fit. The first column corresponds to :math:`\sigma_{\rm DRW}`, the second corresponds to :math:`\tau_{\rm DRW}`, and the rest are the tophat parameters, in the same order as described in ``tophat_params``.
      - list of :python:`float`
    * - ``cont_model``
      - The output ``javelin.lcmodel`` object for the initial continuum fit.
      - ``javelin.lcmodel.Cont_Model``, :python:`None`
    * - ``rmap_model``
      - The output ``javelin.lcmodel`` object for the final fit.
      - ``javelin.lcmodel.Rmap_Model``, ``javelin.lcmodel.Pmap_Model``
    * - ``cont_dat``
      - The continuum light curve in a ``javelin.zylc.LightCurve`` object.
      - ``javelin.zylc.LightCurve``
    * - ``tot_dat``
      - All light curves (continuum +lines) in a ``javelin.zylc.LightCurve`` object.
      - ``javelin.zylc.LightCurve``
    * - ``bestfit_model``
      - The ``javelin.zylc.LightCurve`` object for the JAVELIN fit to the light curves.
      - ``javelin.zylc.LightCurve``


.. note:: If both of the DRW parameters (i.e. the first two parameters) are fixed, the continuum will not be fit to get an estimate on :math:`\sigma_{\rm DRW}` and :math:`\tau_{\rm DRW}`. In this case, the ``cont_hpd`` and ``cont_model`` keys will be :python:`None`.

.. note:: If :python:`rm_type="spec"`, then the ``rmap_model`` key will be a ``javelin.lcmodel.Rmap_Model`` object. If :python:`rm_type="phot"`, then the ``rmap_model`` key will be a ``javelin.lcmodel.Pmap_Model`` object.



Module: Weighting
-----------------

The weighting module output dictionary will contain three dictionaries within it, for the three modules (with keys ``pyccf``, ``javelin``, and ``pyroa``). Each of these two dictionaries will have similar data, representing the results from the weighting. If any of these modules aren't run, then the value corresponding to its key will be :python:`None`.
Similar to the DRW Rejection module, the values for the keys will be lists of results, one for each line, in the order of the input liht curves.

Most of the keys for each of the three dictionaries will be the same:

.. list-table::
    :widths: 20 60 20
    :header-rows: 1

    * - Key
      - Description
      - Type
    * - ``bounds``
      - The bounds and lag value of the primary peak, given as [lower bound, peak, upper bound].
      - list of :python:`float`
    * - ``acf``
      - The ACF of the continuum light curve.
      - list of :python:`float`
    * - ``lags``
      - The lags that the weighting distributions are computed on.
      - list of :python:`float`
    * - ``weight_dist``
      - The weight distribution :math:`w(\tau)`
      - list of :python:`float`
    * - ``smooth_dist``
      - The smoothed :math:`w(\tau)`.
      - list of :python:`float`
    * - ``ntau``
      - The number of overlapping points at a given lag :math:`N(\tau)`.
      - list of :python:`float`
    * - ``frac_rejected``
      - The fraction of the original CCCD rejected when downsampling.
      - list of :python:`float`




The three modules differ with these keys:

.. list-table::
    :widths: 15 15 55 15
    :header-rows: 1

    * - Module
      - Key
      - Description
      - Type
    * - pyCCF
      - ``centroid``
      - The median of the downsampled CCCD and its uncertainties, given as [lower error, value, upper error].
      - list of :python:`float`
    * - pyCCF
      - ``downsampled_CCCD``
      - The downsampled CCCD.
      - list of :python:`float`
    * - pyCCF
      - ``rmax_pyCCF``
      - The maximum value of the downsampled CCCD.
      - list of :python:`float`
    * - JAVELIN
      - ``tophat_lag``
      - The median of the downsampled distribution for the JAVELIN tophat peak.
      - list of :python:`float`
    * - JAVELIN
      - ``downsampled_lag_dist``
      - The downsampled distribution for the JAVELIN tophat peak.
      - list of :python:`float`
    * - JAVELIN
      - ``rmax_javelin``
      - The maximum value of the CCCD within :math:`\pm 1\sigma` of the median tophat lag.
      - list of :python:`float`
    * - PyROA
      - ``time_delay``
      - The median of the downsampled distribution for the PyROA time delay.
      - list of :python:`float`
    * - PyROA
      - ``downsampled_lag_dist``
      - The downsampled distribution for the PyROA time delay.
      - list of :python:`float`
    * - PyROA
      - ``rmax_pyroa``
      - The maximum value of the CCCD within :math:`\pm 1\sigma` of the median PyROA time delay.
      - list of :python:`float`


.. note:: If pyCCF isn't run, all ``rmax_(module)`` values will be NaN.

In addition, the weighting module will output an array of dictionaries representing the data in the weighting summary files (``weighting_summary.fits``) for each line. The keys and values for the entries in these dictionaries can be found in the pyPetal output files API.
