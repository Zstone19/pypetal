=============
Installation
=============

External requirements (not installed by ``pip`` or ``setup.py``):

* A Fortran compiler (>F90)


pyPetal is available on PyPI and can be installed with pip:
::
    pip install pypetal


Or, if you want to install the latest development version:
::
    git clone https://github.com/Zstone19/pypetal.git
    cd pypetal
    pip install .


PLIKE is an optional algorithm that is used in pyPetal. There is a script available in the main directory to install and compile PLIKE. To install PLIKE, run the following command:
::
    sh build_plike.sh

.. note:: PLIKE is not required to run pyPetal, it is an optional part of the pyZDCF module to get an estimate of the time lag. If installing through ``build_plike.sh``, pyPetal assumes that ``gfortran`` is the Fortran compiler used.
