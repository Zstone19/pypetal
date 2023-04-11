pyPetal Output Files and Directories
=====================================

Depending on the modules the user chooses to run, and the parameters chosen for each module, the output diagnostic information from pyPetal, and the structure of the output directiry will be different.
In general, running all of the modules with three lines (named "cont", "line1", and "line2") will produce the following directory structure:

.. code-block:: none

    output_directory/
    ├── cont/
    │   ├── drw_rej/
    │   └── detrend.pdf
    ├── line1/
    │   ├── drw_rej/
    │   ├── pyccf/
    │   ├── pyzdcf/
    │   ├── javelin/
    │   ├── weights/
    │   └── detrend.pdf
    ├── line2/
    │   ├── drw_rej/
    │   ├── pyccf/
    │   ├── pyzdcf/
    │   ├── javelin/
    │   ├── weights/
    │   └── detrend.pdf
    ├── processed_lcs/
    ├── pyroa/
    ├── pyroa_lcs/
    ├── light_curves/
    ├── pyccf_weights_res.pdf
    └── javelin_weights_res.pdf



Each line will have its own subdirectory labeled with the same names given in the ``line_names`` argument for ``pyPetal.pipeline.run_pipeline``. Each of these line subdirectories will have multiple subdirectories for each module, depending on
which modules are run.

The ``processed_lcs`` subdirectory contains the original light curves for all input lines after processing steps (i.e. DRW-based outlier rejection and detrending). This subdirectory will not exist if neither DRW-based outlier rejection nor detrending are run.

The ``light_curves`` subdirectory contains the original light curves for all input lines before processing steps. In addition to the light curves, these files will contain the masks produced in the DRW-based outlier rejection module if it is run.


In addition, each of these subdirectories contain files with the results of modules from pyPetal. We describe the files output from each of these modules below.



Light Curves
------------

Regardless of the modules run, pyPetal will produce a ``light_curves`` subdirectory within the main ``output_directory``. This contains the original light curves obtained through the ``arg2`` argument of ``pypetal.pipeline.run_pipeline``.
These light curves will be named ``{line_name}.dat``, where ``line_name`` is the name of the given light curve input to the pyPetal pipeline. These will be formatted as CSV files, with the first three columns
representing the times, values, and uncertainties of the light curve.

If the DRW-based outlier rejection module is run, these light curve files will contain a fourth column with the DRW-rejection mask. This mask will consist of booleans, where :python:`True` means a point was rejeted.



Module: DRW-based Outlier Rejection
-----------------------------------

The DRW Rejection module file output is unique in that it depends on the user's ``reject_data`` argument. For all lines with :python:`reject_data=True`, their subdirectory will contain a ``drw_rej/`` subdirectory.
This subdirectory will contain the following files:

.. list-table::
    :widths: 30 30 10 30
    :header-rows: 1

    * - Filename
      - Description
      - Format
      - Columns
    * - ``{line_name}_chain.dat``
      - The MCMC chains for the DRW parameters from the fit.
      - CSV
      - :math:`\tau_{\rm DRW}, \sigma_{\rm DRW}, \sigma_n`
    * - ``{line_name}_drw_fit.dat``
      - The DRW fit to the light curve.
      - CSV
      - time, value, uncertainty
    * - ``{line_name}_mask.dat``
      - The DRW-based outlier rejection mask.
      - CSV
      -
    * - ``{line_name}_drw_fit.pdf``
      - A figure describing the DRW fit to the light curve (see the DRW rejection tutorial).
      - PDF
      -

.. note:: If :python:`jitter=False` for the module, there will only be two columns in ``{line_name}_chain.dat``, and the jitter term :math:`sigma_n` will not be included.


In addition, pypetal will save the light curve excluding the rejected points to the ``processed_lcs`` subdirectory under the name ``{line_name}_data.dat``.




Module: Detrending
------------------

There is only file output from the detrending module, which will appear in each line's subdirectory. This will be a plot showing the linear fit to the original light curve before subtraction, which will be named ``detrend.pdf``.

In addition, the detrended light curve will be saved to the ``processed_lcs`` subdirectory under the name ``{line_name}_detrended.dat``.


.. warning:: The detrending module takes place after the DRW rejection module. Therefore, the detrended and rejected results will overwrite the purely rejected results in the ``processed_lcs/`` directory under the same filename.




Module: PyCCF
-------------

Each line subdirectory (excluding the continuum) will contain a subdirectory ``pyccf/`` for all results from the pyCCF module. This subdirectory will contain the following files:

.. list-table::
    :widths: 30 30 10 30
    :header-rows: 1

    * - Filename
      - Description
      - Format
      - Columns
    * - ``{line_name}_ccf_dists.dat``
      - The CCCD and CCPD.
      - CSV
      - CCCD, CCPD
    * - ``{line_name}_ccf.dat``
      - The CCF.
      - CSV
      - Time lags, CCF
    * - ``{line_name}_ccf.pdf``
      - A figure showing the CCF and output pyCCF distributions (see the pyCCF tutorial).
      - PDF
      -



Module: pyZDCF
--------------

Each line subdirectory (excluding the continuum) will contain a subdirectory ``pyzdcf/`` for all results from the pyZDCF module. This subdirectory will contain the following files:

.. list-table::
    :widths: 30 30 10 30
    :header-rows: 1

    * - Filename
      - Description
      - Format
      - Columns
    * - ``{line_name}_{prefix}.dcf``
      - The ZDCF file from pyZDCF.
      - ASCII
      - tau, -sig(tau), +sig(tau), dcf, -err(dcf), +err(dcf), #bin
    * - ``{line_name}_zdcf.pdf``
      - A figure showing the ZDCF (see the pyZDCF tutorial).
      - PDF
      -




Module: PLIKE
-------------

If PLIKE is run under the pyZDCF module, its results will be stored in the ``pyzdcf/`` directory for a given line. It will add the following additional files:

.. list-table::
    :widths: 30 30 10 30
    :header-rows: 1

    * - Filename
      - Description
      - Format
      - Columns
    * - ``{line_name}_plike.out``
      - The PLIKE results.
      - ASCII
      - num, lag, -dr, +dr, r, likelihood



Module: PyROA
--------------

Unlike the previous modules, the layout of the output directory and the structure of the files depend on the ``together`` parameter.

If :python:`together=True`, the output directory for all lines will be ``output_directory/pyroa/``. If :python:`together=False`, each line will have it's PyROA results in its own subdirectory, labeled ``pyroa/``.

In addition, PyROA necessitates a directory for all light curves with names and contents in a specific format. This will be the ``output_directory/pyroa_lcs/`` directory.

Each PyROA directory (whether ``together`` is :python:`True` or :python:`False`) will have the following files:

.. list-table::
    :widths: 30 30 10 30
    :header-rows: 1

    * - Filename
      - Description
      - Format
      - Columns
    * - ``samples.obj``
      - The PyROA MCMC samples.
      - pickle
      - see below
    * - ``samples_flat.obj``
      - The PyROA MCMC samples, flattened.
      - pickle
      - see below
    * - ``Lightcurve_models.obj``
      - The models for the light curves (including the continuum).
      - pickle
      - There will be one model for each light curve, and each model with have the time, value, and error for the modeled light curve.
    * - ``X_t.obj``
      - The drving continuum light curve model.
      - pickle
      - time, value, error
    * - ``trace_plot.pdf``
      - A figure showing the MCMC trace plots for each parameter, and the cutoff for the specified burn-in.
      - PDF
      -
    * - ``histogram_plot.pdf``
      - A figure showing the MCMC posterior histograms for each parameter (excluding burn-in).
      - PDF
      -
    * - ``corner_plot.pdf``
      - A figure showing the MCMC corner plot for all parameters (excluding burn-in).
      - PDF
      -
    * - ``fits_plot.pdf``
      - A figure analogous to the PyROA fit plots, showing the light curve fits to the data, the time lag distributions, and the ``delay_dist`` distributions (if :python:`delay_dist=True`).
      - PDF
      -



If :python:`together=True`, the columns of the ``samples`` files will be:

.. list-table::
    :widths: 30 30 10 30
    :header-rows: 1

    * - ``add_var``
      - ``delay_dist``
      - Columns
    * - :python:`False`
      - :python:`False`
      - :math:`A_0, B_0, \tau_0, A_1, B_1, \tau_1, ..., \Delta`
    * - :python:`True`
      - :python:`False`
      - :math:`A_0, B_0, \tau_0, \sigma_0, A_1, B_1, \tau_1, \sigma_0, ..., \Delta`
    * - :python:`False`
      - :python:`True`
      - :math:`A_0, B_0, \tau_0, A_1, B_1, \tau_1, \Delta_1, A_2, B_2, \tau_2, \Delta_2, ..., \Delta`
    * - :python:`True`
      - :python:`True`
      - :math:`A_0, B_0, \tau_0, \sigma_0, A_1, B_1, \tau_1, \Delta_1, \sigma_1, A_2, B_2, \tau_2, \Delta_2, \sigma_2, ..., \Delta`

If :python:`together=False`, the columns will be the same as for :python:`together=True`, except the file for each line will only contain samples for the continuum, and that line.



Module: JAVELIN
---------------

Unlike the other modules, the layout of the output directory and the structure of the files depends on multiple parameters, in particular ``together``, ``rm_type``, and ``fixed/p_fix``.

If :python:`together=True`, the output directory for all lines will be ``output_directory/javelin/``. If :python:`together=False`, each line will have it's JAVELIN results in its own subdirectory, labeled ``javelin/``.

If :python:`together=True`, the output directory will contain the following files:

.. list-table::
    :widths: 30 30 10 30
    :header-rows: 1

    * - Filename
      - Description
      - Format
      - Columns
    * - ``burn_cont.txt``
      - The burn-in samples for the initial continuum fit.
      - ASCII
      - :math:`\log(\sigma_{\rm DRW}), \log(\tau_{\rm DRW})`
    * - ``burn_rmap.txt``
      - The burn-in sampled for the total JAVELIN fit.
      - ASCII
      - :math:`\log(\sigma_{\rm DRW}), \log(\tau_{\rm DRW})`, tophat parameters for each line
    * - ``chain_cont.txt``
      - The MCMC chains for the initial continuum fit.
      - ASCII
      - :math:`\log(\sigma_{\rm DRW}), \log(\tau_{\rm DRW})`
    * - ``chain_rmap.txt``
      - The MCMC chains for the total JAVELIN fit.
      - ASCII
      - :math:`\log(\sigma_{\rm DRW}), \log(\tau_{\rm DRW})`, tophat parameters for each line
    * - ``logp_cont.txt``
      - The log-probability for the initial continuum fit.
      - ASCII
      -
    * - ``logp_rmap.txt``
      - The log-probability for the total JAVELIN fit.
      - ASCII
      -
    * - ``cont_lcfile.dat``
      - The continuum light curve in JAVELIN format.
      - ASCII
      -
    * - ``tot_lcfile.dat``
      - All light curves in JAVELIN format.
      - ASCII
      -
    * - ``{line_name}_lc_fits.dat``
      - The best-fit light curves for each line. There will be one file for each line.
      - CSV
      - time, value, uncertainty
    * - ``javelin_histogram.pdf``
      - A figure showing the histograms of the MCMC chains for each parameter.
      - PDF
      -
    * - ``javelin_bestfit.pdf``
      - A figure showing the best-fit light curves for each line.
      - PDF
      -
    * - ``javelin_corner.pdf``
      - A corner plot for all JAVELIN parameters.
      - PDF
      -




If :python:`together=False`, the output directory for each line will contain the following files:

.. list-table::
    :widths: 30 30 10 30
    :header-rows: 1

    * - Filename
      - Description
      - Format
      - Columns
    * - ``burn_cont.txt``
      - The burn-in samples for the initial continuum fit.
      - ASCII
      - :math:`\log(\sigma_{\rm DRW}), \log(\tau_{\rm DRW})`
    * - ``burn_rmap.txt``
      - The burn-in sampled for the total JAVELIN fit.
      - ASCII
      - :math:`\log(\sigma_{\rm DRW}), \log(\tau_{\rm DRW})`, tophat parameters for the line
    * - ``chain_cont.txt``
      - The MCMC chains for the initial continuum fit.
      - ASCII
      - :math:`\log(\sigma_{\rm DRW}), \log(\tau_{\rm DRW})`
    * - ``chain_rmap.txt``
      - The MCMC chains for the total JAVELIN fit.
      - ASCII
      - :math:`\log(\sigma_{\rm DRW}), \log(\tau_{\rm DRW})`, tophat parameters for the line
    * - ``logp_cont.txt``
      - The log-probability for the initial continuum fit.
      - ASCII
      -
    * - ``logp_rmap.txt``
      - The log-probability for the total JAVELIN fit.
      - ASCII
      -
    * - ``cont_lcfile.dat``
      - The continuum light curve in JAVELIN format.
      - ASCII
      -
    * - ``tot_lcfile.dat``
      - All light curves in JAVELIN format.
      - ASCII
      -
    * - ``{line_name}_lc_fits.dat``
      - The best-fit light curves for the line.
      - CSV
      - time, value, uncertainty
    * - ``javelin_histogram.pdf``
      - A figure showing the histograms of the MCMC chains for each parameter.
      - PDF
      -
    * - ``javelin_bestfit.pdf``
      - A figure showing the best-fit light curves for each line.
      - PDF
      -
    * - ``javelin_corner.pdf``
      - A corner plot for all JAVELIN parameters.
      - PDF
      -



.. note:: If both DRW parameters (i.e. the first two) are fixed, then there will not be a ``burn_cont.txt`` or ``chain_cont.txt`` file.

.. note:: If any parameters are fixed, there will not be a ``javelin_corner.pdf`` file.

The number of tophat parameters in the ``burn`` and ``chain`` files depends on the ``rm_type`` argument. If :python:`rm_type="spec"`, there will be 3 tophat parameters for each line (t, w, s).
If :python:`rm_type="phot"`, there will be 2 tophat parameters for each line (t, w, s, :math:`\alpha`).

If :python:`together=True`, the tophat parameters will be grouped by line in order. For example, if :python:`rm_type="spec"`, the columns of the ``chain`` and ``burn`` files will be
:math:`\log(\sigma_{\rm DRW}), \log(\tau_{\rm DRW}), t_1, w_1, s_1, t_2, w_2, s_2, ...`.




Module: Weighting
-----------------

The output of the weighting module depends on if the pyCCF and JAVELIN modules are run. All results will either be stored in the ``weights/`` subdirectory for each line or the main ``output_directory/``.

If the pyCCF module is run, the ``weights/`` subdirectory will contain the following files:

.. list-table::
    :widths: 30 30 10 30
    :header-rows: 1

    * - Filename
      - Description
      - Format
      - Columns
    * - ``pyccf_weights.dat``
      - The distributions needed to weight the CCCD for the line.
      - CSV
      - lags :math:`\tau` , :math:`N(\tau)`, :math:`w(\tau)`, ACF, smoothed CCCD, smoothed weighted CCCD
    * - ``pyccf_weighted_cccd.dat``
      - The downsampled CCCD after weighting and finding the primary peak.
      - CSV
      -

If the JAVELIN module is run, the ``weights/`` subdirectory will contain the following files:

.. list-table::
    :widths: 30 30 10 30
    :header-rows: 1

    * - Filename
      - Description
      - Format
      - Columns
    * - ``javelin_weights.dat``
      - The distributions needed to weight the JAVELIN lag distribution :math:`t` for the line.
      - CSV
      - lags :math:`\tau` , :math:`N(\tau)`, :math:`w(\tau)`, ACF, smoothed :math:`t`, smoothed weighted :math:`t`
    * - ``javelin_weighted_lag_dist.dat``
      - The downsampled :math:`t` after weighting and finding the primary peak.
      - CSV
      -



If the PyROA module is run, the ``weights/`` subdirectory will contain the following files:

.. list-table::
    :widths: 30 30 10 30
    :header-rows: 1

    * - Filename
      - Description
      - Format
      - Columns
    * - ``pyroa_weights.dat``
      - The distributions needed to weight the PyROA lag distribution :math:`t` for the line.
      - CSV
      - lags :math:`\tau` , :math:`N(\tau)`, :math:`w(\tau)`, ACF, smoothed :math:`t`, smoothed weighted :math:`t`
    * - ``pyroa_weighted_lag_dist.dat``
      - The downsampled :math:`t` after weighting and finding the primary peak.
      - CSV
      -


In addition, the weighting module will always output the following files in the ``weights/`` subdirectory:

.. list-table::
    :widths: 30 30 10 30
    :header-rows: 1

    * - Filename
      - Description
      - Format
      - Columns
    * - ``{line_name}_weights.pdf``
      - A figure showing the distributions needed to weight the CCCD, JAVELIN lag distribution, or PyROA lag distribution.
      - PDF
      -
    * - ``weight_summary.fits``
      - A FITS table containing the results of the weighting and auxiliary information from the weighting.
      - FITS
      - See below



The ``weight_summary.fits`` file contains the following information for each module (pyCCF, JAVELIN, and/or PyROA):

.. list-table::
    :widths: 30 30 30
    :header-rows: 1

    * - Name
      - Description
      - Type
    * - ``k``
      - The exponent used to calculate :math:`P(\tau)`
      - :python:`float`
    * - ``n0_(module)``
      - The value of :math:`N(0)`. Given for both the CCCD and :math:`t`.
      - :python:`float`
    * - ``peak_bounds_(module)``
      - The bounds of the primary peak of the weighted distribution. Given as [lower bound, peak, upper bound] for both the CCCD and :math:`t`.
      - list of :python:`float`
    * - ``peak_(module)``
      - The peak of the primary peak. Given for both the CCCD and :math:`t`.
      - :python:`float`
    * - ``lag_(module)``
      - The median of the downsampled lag distribution. Given for both the CCCD and :math:`t`.
      - :python:`float`
    * - ``lag_err_(module)``
      - The uncertainty on the lag. Given as [lower error, upper error] for both the CCCD and :math:`t`.
      - list of :python:`float`
    * - ``frac_rejected_(module)``
      - The fraction of the original distribution that was rejected to obtain the downsampled distribution. Given for both the CCCD and :math:`t`
      - :python:`float`
    * - ``rmax_(module)``
      - The maximum value of the CCCD within the region covered by the downsampled JAVELIN lag distribution.
      - :python:`float`

where ``(module)`` is either ``pyccf``, ``javelin``, or ``pyroa``.

.. note:: If a module is not run, its values in ``weight_summary.txt`` for that module will be ``NaN``.

.. note:: If pyCCF is not run, ``rmax_(module)`` will be ``NaN``.


In addition, the following files will be placed in the main ``output_directory/``:

.. list-table::
    :widths: 30 30 10 30
    :header-rows: 1

    * - Filename
      - Description
      - Format
      - Columns
    * - ``pyccf_weights_res.pdf``
      - A figure showing the output of the weighting process for the CCCD.
      - PDF
      -
    * - ``javelin_weights_res.pdf``
      - A figure showing the output of the weighting process for the JAVELIN lag distribution.
      - PDF
      -
    * - ``pyroa_weights_res.pdf``
      - A figure showing the output of the weighting process for the PyROA lag distribution.
      - PDF
      -
