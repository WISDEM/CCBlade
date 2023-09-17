# CCBlade

A blade element momentum method for analyzing wind turbine aerodynamic performance that is robust (guaranteed convergence), fast (superlinear convergence rate), and smooth (continuously differentiable).  Analytic gradients are also (optionally) provided for the distributed loads, thrust, torque, and power with respect to design variables of interest.

Author: [NREL WISDEM Team](mailto:systems.engineering@nrel.gov) 

## Documentation

See local documentation in the `docs`-directory or access the online version at <http://wisdem.github.io/CCBlade/>

## Prerequisites

CCBlade execution requires: `numpy`, `scipy`, `openmdao`
CCBlade installation requires: `meson`, `ninja`, `gfortran`

## Installation

CCBlade is available as a [WISDEM](https://github.com/WISDEM/WISDEM) module and WISDEM is both pip-installable and conda-installable. For building CCBlade from source as a standalone library, first make sure that you have the necessary prerequisites installed.  After cloning the repository, do:

    $ pip install CCBlade

## Run Unit Tests

To check if installation was successful, run the unit tests

    $ python test/test_ccblade.py
    $ python test/test_gradients.py

For software issues please use <https://github.com/WISDEM/CCBlade/issues>.  For functionality and theory related questions and comments please use the NWTC forum for [Systems Engineering Software Questions](https://wind.nrel.gov/forum/wind/viewtopic.php?f=34&t=1002).

