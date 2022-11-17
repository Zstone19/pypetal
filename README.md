PETL: A Pipeline for Estimating Time Lags
==========================================

Installation
------------
A number of different packages and softwares are required. Most can be installed through `pip` or `python setup.py install`. Only a few are not available through `pip` and must be installed manually.
```
    A Fortran compiler (>F90)
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

If an error occurs during installation, building JAVELIN, pyCCF, and PLIKE can be done manually with the `build_dep.sh` script. For example
```
    sh build_dep.sh -u false -f gnu95 -p true
```
where `u`, `f`, and `p` have the same meanings as before.