.. petl documentation master file, created by
   sphinx-quickstart on Wed Nov 16 00:29:48 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PETL: A Pipeline for Estimating AGN Time Lags
===============================================
**PETL** is a pipeline made in python to obtain time lags from multi-band AGN time-series data. Normally, AGN geometric and kinematic analyses (e.g. reverberation mapping) utilize a variety of different tools to obtain time lags between two (or more) light curves. This package combines three popular algorithms for estimating time lags (pyCCF, pyZDCF, and JAVELIN), and uses the popular Damped Random Walk algorithm to model input light curves for outlier rejection.  
Currently, PETL has combined the functionality of pyCCF, pyZDCF, and JAVELIN to produce cross-correlation functions, discrete correlation functions, and mean time lags. This is only made to run on Linux-based operating systems, though this may be improved in the future.

Installation
------------
PETL is available on PyPI, and can be installed with pip:
```
pip install petl
```

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
   :caption: API:

   python/drw_fit
   python/utils
   python/modules
   python/pipeline
   python/plot
