.. pypetal documentation master file, created by
   sphinx-quickstart on Wed Nov 16 00:29:48 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=================================================
pyPETaL: A Pipeline for Estimating AGN Time Lags
=================================================
**pyPetal** is a pipeline made in python to obtain time lags from multi-band AGN time-series data. Normally, AGN geometric and kinematic analyses (e.g. reverberation mapping) utilize a variety of different tools to obtain time lags between two (or more) light curves. This package combines four popular algorithms for estimating time lags (pyCCF, pyZDCF, PyROA, and JAVELIN), and uses the popular Damped Random Walk algorithm to model input light curves for outlier rejection.

Currently, pyPetal has combined the functionality of pyCCF, pyZDCF, PyROA, MICA2m and JAVELIN to produce cross-correlation functions, discrete correlation functions, and mean time lags. This is only made to run on Linux-based operating systems, though this may be improved in the future.


.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   installation
   notebooks/getting_started


.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   notebooks/drw_rej
   notebooks/detrending
   notebooks/pyccf
   notebooks/pyzdcf
   notebooks/plike
   pyroa_toc
   mica2_toc
   javelin_toc
   notebooks/weighting
   notebooks/all_together
   notebooks/load
   notebooks/from_file
   notebooks/from_bash
   notebooks/bibliography


.. toctree::
   :maxdepth: 1
   :caption: API

   pl_args
   pl_output
   pl_files
   python/drw_rej
   python/detrending
   python/pyccf
   python/pyzdcf
   python/pyroa
   python/javelin
   python/weighting
   python/load
   python/plot




Citing pyPetal
--------------

To cite the pyPetal code itself, use the ASCL reference: `</https://ascl.net/2401.004>`__

To cite the paper pyPetal was used in: `<https://ui.adsabs.harvard.edu/abs/2023arXiv230501014S>`__




References
----------

pyPetal makes use of a multitude of packages:

* **pyCCF**: `<https://arxiv.org/abs/astro-ph/9802103>`__, `<http://ascl.net/code/v/1868>`__
* **pyZDCF**: `<https://doi.org/10.5281/zenodo.7253034>`__,
* **PLIKE**: `<http://arxiv.org/abs/1302.1508>`__
* **PyROA**: `<https://ui.adsabs.harvard.edu/abs/2021MNRAS.508.5449D>`__
* **MICA2**: `<https://ui.adsabs.harvard.edu/abs/2016ApJ...831..206L>`__
* **JAVELIN**: `<http://adsabs.harvard.edu/abs/2013ApJ...765..106Z>`__, `<http://adsabs.harvard.edu/abs/2011ApJ...735...80Z>`__, `<http://adsabs.harvard.edu/abs/2016ApJ...819..122Z>`__
