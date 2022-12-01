pyPETaL: A Pipeline for Estimating AGN Time Lags
=================================================

[![Documentation Status](https://readthedocs.org/projects/petl1/badge/?version=latest)](https://petl1.readthedocs.io/en/latest/?badge=latest)

Installation
------------
External requirements (not installed by ``pip`` or ``setup.py``):
```
    A Fortran compiler (>F90)
```

pyPetal is available on PyPI and can be installed with pip: **(NOT IMPLEMENTED YET)**
```
    pip install pypetal
```

Or, if you want to install the latest development version:
```
    git clone https://github.com/Zstone19/pypetal.git
    cd pypetal
    pip install .
```

PLIKE is an optional algorithm that is used in pyPetal. There is a script available in the main directory to install and compile PLIKE (assuming that `gfortran` is installed). To install PLIKE, run the following command:
```
    sh build_plike.sh
```