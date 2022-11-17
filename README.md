PETL: A Pipeline for Estimating Time Lags
==========================================

Installation
------------
A number of different packages and softwares are required. Most can be installed through `pip` or `python setup.py install`. Only a few are not available through `pip` and must be installed manually.
```
    A Fortran compiler (>F90)
    make
    wget
```

PETL is available on PyPI and can be installed with pip:
```
    pip install petl
```

Or, if you want to install the latest development version:
```
    git clone https://github.com/Zstone19/petl.git
    cd petl
    pip install .
```

PETL will be installed using the default parameters:
```
    python setup.py install -u True -p False -f None
```

The different arguments are:
```
    -u: if JAVELIN is installed locally
    -f: Fortran compiler to use for JAVELIN
    -p: if PLIKE is installed
```

**NOTE**: PLIKE assumes that the Fortran compiler used is `gfortran`. If you have another Fortran compiler, PLIKE may need to be installed/compiled manually from the source. However, PLIKE is an optional part of PETL in general, and may not need to be installed generally.