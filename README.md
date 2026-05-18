# pyPETaL: A Pipeline for Estimating AGN Time Lags

[![Documentation Status](https://readthedocs.org/projects/pypetal/badge/?version=latest)](https://pypetal.readthedocs.io/en/latest/?badge=latest)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://github.com/pre-commit/pre-commit)
[![Workflow Status](https://img.shields.io/github/actions/workflow/status/Zstone19/pypetal/python-package.yml)](https://img.shields.io/github/actions/workflow/status/Zstone19/pypetal/python-package.yml)
[![License](https://img.shields.io/github/license/Zstone19/pypetal)](https://img.shields.io/github/license/Zstone19/pypetal)
[![codecov](https://codecov.io/gh/Zstone19/pypetal/branch/main/graph/badge.svg?token=00O40N9H05)](https://codecov.io/gh/Zstone19/pypetal)

[![pypi-ver](https://img.shields.io/pypi/v/pypetal)](https://img.shields.io/pypi/v/pypetal)
[![python-ver](https://img.shields.io/pypi/pyversions/pypetal)](https://img.shields.io/pypi/pyversions/pypetal)
[![pypi-downloads](https://static.pepy.tech/badge/pypetal)](https://pepy.tech/project/pypetal)

<a href="https://ascl.net/2401.004"><img src="https://img.shields.io/badge/ascl-2401.004-blue.svg?colorB=262255" alt="ascl:2401.004" /></a>


pyPETaL is a time-series data analysis pipeline for AGN reverberation mapping (RM) data. It combines multiple different popular softwares using for AGN RM analysis, including PyCCF, PyZDCF, JAVELIN, and PyROA. This pipeline also implements outlier rejection using Damped Random Walk Gaussian process fitting, and detrending through the LinMix algorithm. pyPetal implements a weighting scheme (for all modules) in order to mitigate aliasing in peaks of time lag distributions between light curves.

pyPetal is very flexible, with almost every argument for each module allowing user input. pyPetal is designed to work with any combination of modules being run, allowing it to scale from the simplest to the most complex of projects.



## Installation

### pyPetal

pyPetal is available on PyPI and can be installed with pip:
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



### Detrending

pyPetal offers the option to detrend the input light curves via the [LinMix](https://github.com/jmeyers314/linmix.git) algorithm. This package is not offered in the base version of pyPetal, but can be installed with pip:
```
    pip install "linmix @ git+https://github.com/jmeyers314/linmix.git"
```

Or with ``poetry``:
```
    poetry install --with extra
```


### MICA2

pyPetal also offers [MICA2](https://github.com/LiyrAstroph/MICA2) as an optional module. This package is more complex to install than the others - to find out more, read the [MICA2 installation guide](https://mica2.readthedocs.io/en/latest/getting_started.html).

``poetry`` can help to install some easily-installable python dependencies for MICA2:
```
    poetry install --with extra
```

The functionality and inputs of the MICA2 module are identical to the original software, so the best way to learn how this module functions is to read the MICA2 documentation.

__NOTE__: Seeing as this is an optional module, pyPetal will still compile and run if MICA2's dependencies are not installed. Each pyPetal run will assume ``run_mica2=False``.



### pyPetal and JAVELIN

The JAVELIN software used in pyPetal runs on Python 2, though the bulk of pyPetal software relies on Python >=3.8. To circumvent this issue, pyPetal has a JAVELIN "module" (``pypetal-jav``) that can be installed as a separate package and used in conjunction with pyPetal, in a Python 2 environment.


External requirements (not installed by ``pip`` or ``setup.py``):
```
    A Fortran compiler (>F90)
```


pyPetal-jav is available on PyPI and can be installed with pip:
```
    pip install --no-deps pypetal-jav
    pip install pypetal-jav
```


Or, if you want to install the latest development version:
```
    git clone https://github.com/Zstone19/pypetal-jav.git
    cd pypetal-jav
    pip install .
```


__NOTE:__ The user may need to install ``NumPy`` before installing pyPetal-jav through ``pip`` or ``setup.py``. This is because ``JAVELIN`` requires ``NumPy`` in order to be installed. This can be done with:
```
    pip install numpy
```



## Citing pyPetal

To cite the pyPetal code itself, use the [ASCL reference](https://ascl.net/2401.004):
```
    Stone Z., 2024, version xxxx, Astrophysics Source Code Library, record ascl:2401.004
```

Cite the paper pyPetal was used in: [Shen et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230501014S/abstract)
