.. petl documentation master file, created by
   sphinx-quickstart on Wed Nov 16 00:29:48 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===============================================
PETL: A Pipeline for Estimating AGN Time Lags
===============================================
**PETL** is a pipeline made in python to obtain time lags from multi-band AGN time-series data. Normally, AGN geometric and kinematic analyses (e.g. reverberation mapping) utilize a variety of different tools to obtain time lags between two (or more) light curves. This package combines three popular algorithms for estimating time lags (pyCCF, pyZDCF, and JAVELIN), and uses the popular Damped Random Walk algorithm to model input light curves for outlier rejection.  

Currently, PETL has combined the functionality of pyCCF, pyZDCF, and JAVELIN to produce cross-correlation functions, discrete correlation functions, and mean time lags. This is only made to run on Linux-based operating systems, though this may be improved in the future.


.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   installation


.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   notebooks/drw_rej
   notebooks/pyccf
   notebooks/pyzdcf
   notebooks/javelin
   notebooks/all_together


.. toctree::
   :maxdepth: 1
   :caption: API

   python/drw_rej
   python/pyccf
   python/pyzdcf
   python/javelin
   python/pipeline
   python/plot



References
----------

PETL makes use of a multitude of packages:

* **pyCCF**: `<https://arxiv.org/abs/astro-ph/9802103>`__, `<http://ascl.net/code/v/1868>`__
* **pyZDCF**: `<https://doi.org/10.5281/zenodo.7253034>`__, 
* **PLIKE**: `<http://arxiv.org/abs/1302.1508>`__
* **JAVELIN**: `<http://adsabs.harvard.edu/abs/2013ApJ...765..106Z>`__, `<http://adsabs.harvard.edu/abs/2011ApJ...735...80Z>`__, `<http://adsabs.harvard.edu/abs/2016ApJ...819..122Z>`__
