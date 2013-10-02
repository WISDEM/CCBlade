# CCBlade

A blade element momentum method for analyzing wind turbine aerodynamic performance that is robust (guaranteed convergence), fast (superlinear convergence rate), and smooth (continuously differentiable).  Analytic gradients are also (optionally) provided for the distributed loads, thrust, torque, and power with respect to design variables of interest.

Author: [S. Andrew Ning](mailto:andrew.ning@nrel.gov)

## Prerequisites

C compiler, Fortran compiler, NumPy, SciPy

## Installation

Install CCBlade with the following command.

    $ python setup.py install

Note that the installation also includes AirfoilPrep.py.  Though not strictly necessary to use with CCBlade, it is convenient when working with AeroDyn input files or doing any aerodynamic preprocessing of airfoil data.

## Run Unit Tests

To check if installation was successful, run the unit tests

    $ python test/test_ccblade.py
    $ python test/test_gradients.py

## Detailed Documentation

Open `docs/_build/html/index.html` in your browser.


## User Information

If you download this software directly from GitHub without going through our website (not yet setup at the moment), we would appreciate it if you could report your user information.  This helps us understand how our software is being used, determine where to direct resources for further development, and notify users of critical updates.  For each user please send an email to [Andrew Ning](mailto:andrew.ning@nrel.gov?subject=CCBlade%20User) with the following information:

- First and last name
- Email address
- Organization name
- Location (Country, if U.S. then State)
- Organization Type (Certifier, Consultant, Individual, Manufacturer, Laboratory, School)
- Intended use (Wind Power, Marine and Hydrokinetic Energy, Other)