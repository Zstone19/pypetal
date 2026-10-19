=======
MICA2
=======

MICA2 is a more complex tool used to correlate two light curves in a non-parametric fashion to find the transfer function between the two,
as a sum of either Gaussians or tophat functions. MICA2 is much more flexible in its functionality, though this comes with much more complex
customization from the user in terms of parameters chosen and types of analysis possible. These tutorials provide instructions
on in-depth usage of MICA2 and understanding of its arguments and outputs.

While these tutorials do a sufficient job of introducing the user to the basic functionality of MICA2, the user is encouraged to read the
MICA2 `documentation <https://mica2.readthedocs.io/en/latest/index.html>`_ and/or `repository <https://github.com/LiyrAstroph/MICA2>`_ for
more in-depth information on the software. In addition, the user can read the associated `paper <https://ui.adsabs.harvard.edu/abs/2016ApJ...831..206L/abstract>`_
for any additional details.

.. warning::
    MICA2 is a software that can normally be run on multiple cores by using ``mpiexec``, such as with::

        mpiexec -n 4 python mica2_script.py


    pyPetal extends this functionality, so that the MICA2 module can utilize ``mpiexec``. However, if this is done, all other module will
    only be run on one core (as they are optimized for python multi-threading instead of MPI). If the user wishes to only use the MICA2 module,
    it may be best to run using ``mpiexec``. If the user wishes to use the other modules, it may be best to run the MICA2 module separately
    from the other modules.


.. toctree::
    :maxdepth: 1
    :caption: MICA2 Tutorials

    notebooks/mica2/basic
    notebooks/mica2/numcomp
    notebooks/mica2/noorder
    notebooks/mica2/together
