A blade element momentum method for analyzing wind turbine aerodynamic performance that is robust (guaranteed convergence), fast (superlinear convergence rate), and smooth (continuously differentiable).  Analytic gradients are also (optionally) provided for the distributed loads, thrust, torque, and power with respect to design variables of interest.

Author: [S. Andrew Ning](mailto:andrew.ning@nrel.gov)

## User Information

If you came to this page directly without going through the NWTC Information Portal, **we would appreciate if you could [report your user information](http://wind.nrel.gov/designcodes/simulators/ccblade/downloaders/CCBlade_github_redirect.html) before cloning the repository**.  We use this information in order to allocate resources for supporting our software, and to notify users of critical updates.

## Prerequisites

Fortran compiler, NumPy, SciPy, zope.interface

## Installation

Install CCBlade with the following command.

    $ python setup.py install

Note that the installation also includes AirfoilPrep.py.  Though not strictly necessary to use with CCBlade, it is convenient when working with AeroDyn input files or doing any aerodynamic preprocessing of airfoil data.

## Run Unit Tests

To check if installation was successful, run the unit tests

    $ python test/test_ccblade.py
    $ python test/test_gradients.py

## Detailed Documentation

Open a local copy of the documentation at `docs/_build/html/index.html`.  Or access the online version at <http://wisdem.github.io/CCBlade/>


