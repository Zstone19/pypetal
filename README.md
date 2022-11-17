PETL: A Pipeline for Estimating Time Lags
==========================================

Installation
------------
A number of different packages and softwares are required. Most can be installed through `pip` or `python setup.py install`. Only a few are not available through `pip` and must be installed manually.
```
    A Fortran compiler (>F90)
```

PETL is available on PyPI and can be installed with pip: **(NOT IMPLEMENTED YET)**
```
    pip install petl
```

Or, if you want to install the latest development version:
```
    git clone https://github.com/Zstone19/petl.git
    cd petl
    python setup.py install
```

PLIKE is an optional algorithm that is used in PETL. There is a scipt available in PETL to install and compile PLIKE (assuming that `gfortran` is installed). To install PLIKE, run the following command:
```
    sh build_plike.sh
```