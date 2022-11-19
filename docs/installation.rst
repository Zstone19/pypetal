=============
Installation
=============

A number of different packages and softwares are required. Most can be installed through ``pip`` or ``python setup.py install``. Only a few are not available through ``pip`` and must be installed manually.
    #. A Fortran compiler (>F90)


PETL is available on PyPI and can be installed with pip:
:: 
    pip install petl


Or, if you want to install the latest development version:
::
    git clone https://github.com/Zstone19/petl.git
    cd petl
    pip install .


PLIKE is an optional algorithm that is used in PETL. There is a script available in the main directory to install and compile PLIKE. To install PLIKE, run the following command:
::
    sh build_plike.sh

.. note:: PLIKE is not required to run PETL, it is an optional part of the pyZDCF module to get an estimate of the time lag. If installing through ``build_plike.sh``, PETL assumes that ``gfortran`` is the Fortran compiler used.
