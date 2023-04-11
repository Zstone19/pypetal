=========
JAVELIN
=========

JAVELIN is a well-known tool used to calculate the time lag between two AGN light curves.
However, JAVELIN usage and understanding is more complex than the other modules. We
provide in-depth tutorials on how to use the JAVELIN module for pyPetal here.

Installation
------------
Using JAVELIN requires a separate package than the rest of pyPetal, as it runs on Python 2. To circumvent this, we have created ``pypetal-jav``.
To install via ``pip``:
::
    pip install pypetal-jav

.. note:: ``pypetal-jav`` requires a Python 2 environment.
.. note:: ``pypetal-jav`` may require that ``NumPy`` be installed beforehand. To do this with ``pip``:
    ::
        pip install numpy


.. toctree::
    :maxdepth: 1
    :caption: JAVELIN Tutorials

    notebooks/javelin/basic
    notebooks/javelin/together
    notebooks/javelin/fixed
    notebooks/javelin/type
    notebooks/javelin/advanced
