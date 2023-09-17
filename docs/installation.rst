Installation
------------

CCBlade is available as a `WISDEM <https://github.com/WISDEM/WISDEM>`_ module and WISDEM is both pip-installable and conda-installable. For building CCBlade from source as a standalone library, first make sure that you have the necessary prerequisites installed.  After cloning the repository, do:

.. code-block:: bash

   $ pip install CCBlade

The pre-requisites for running CCBlade include: `numpy`, `scipy`, `openmdao`

The pre-requisites for building CCBlade from source include: `meson`, `ninja`, `gfortran`

To check if installation was successful run the unit tests for the NREL 5-MW model

.. code-block:: bash

   $ python test/test_ccblade.py

Additional tests for the gradients are available at:

.. code-block:: bash

   $ python test/test_gradients.py

An "OK" signifies that all the tests passed.

.. only:: latex

    To access an HTML version of this documentation that contains further details and links to the source code, open docs/index.html.

A current copy of the documentation for this code is also available online at http://nrel-wisdem.github.io/CCBlade


.. note::

    The CCBlade installation also installs the `AirfoilPrep` module. Although it is not strictly necessary to use AirfoilPrep.py with CCBlade, its inclusion is convenient when working with AeroDyn input files or doing any aerodynamic preprocessing of airfoil data.

