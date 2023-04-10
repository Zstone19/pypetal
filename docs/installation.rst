=============
Installation
=============


pyPetal
-------

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





pyPetal and JAVELIN
--------------------

The JAVELIN software used in pyPetal runs on Python 2, though the bulk of pyPetal software relies on Python >=3.8. To circumvent this issue, pyPetal has a JAVELIN "module" (``pypetal-jav``) that can be installed as a separate package and used in conjunction with pyPetal, in a Python 2 environment.


External requirements (not installed by ``pip`` or ``setup.py``):

* A Fortran compiler (>F90)


pyPetal-jav is available on PyPI and can be installed with pip:
::
    pip install pypetal-jav


Or, if you want to install the latest development version:
::
    git clone https://github.com/Zstone19/pypetal-jav.git
    cd pypetal-jav
    pip install .


NOTE: The user may need to install ``NumPy`` before installing ``pypetal-jav`` through ``pip`` or ``setup.py``. This is because ``JAVELIN`` requires ``NumPy`` in order to be installed. This can be done with:
::
    pip install numpy
