A blade element momentum method for analyzing wind turbine aerodynamic performance that is robust (guaranteed convergence), fast (superlinear convergence rate), and smooth (continuously differentiable).  Analytic gradients are also (optionally) provided for the distributed loads, thrust, torque, and power with respect to design variables of interest.

Author: [S. Andrew Ning](mailto:andrew.ning@nrel.gov)

## Detailed Documentation

Open a local copy of the documentation at `docs/_build/html/index.html`.  Or access the online version at <http://wisdem.github.io/CCBlade/>

## Prerequisites

Fortran compiler, NumPy, SciPy, zope.interface

## Installation

### GNU/Linux (Ubuntu 16.04 LTS)
Install CCBlade with the following command. It is recommended to use a clean virtual environment. However, we need to install Numpy manually as CCBlade uses its enhanced distutils in order to compile the Fortran routines. Pip will take care of other dependencies.
```
    $ virtualenv devenv
    $ source devenv/bin/activate
    (devenv)$ pip install numpy==1.12.1
    (devenv)$ pip install -r requirements.txt
    (devenv)$ ccblade
```
Note that the installation also includes AirfoilPrep.py.  Though not strictly necessary to use with CCBlade, it is convenient when working with AeroDyn input files or doing any aerodynamic preprocessing of airfoil data.

#### Windows 10
##### Python 2.7
- [Anaconda2 4.4.0](https://repo.continuum.io/archive/Anaconda2-4.4.0-Windows-x86_64.exe)
- [Microsoft Visual C++ Compiler for Python 2.7](http://aka.ms/vcpython27)

Once both are installed, proceed with an Anaconda Prompt
```prompt
(Anaconda2)> conda create --name devenv python=2.7 scipy=0.18.1 zope.interface=4.4.2 m2w64-gcc-fortran=5.3.0
(Anaconda2)> activate devenv
(devenv)> pip install -r requirements.txt
(devenv)> ccblade
```

##### Python 3.6
- [Anaconda3 4.4.0](https://repo.continuum.io/archive/Anaconda2-4.4.0-Windows-x86_64.exe)
- [Microsoft Visual C++ 2015 Build Tools](http://landinghub.visualstudio.com/visual-cpp-build-tools)

Once both are installed, proceed with an Anaconda Prompt
```prompt
(Anaconda3)> conda create --name ccblade-devenv python=3.6 scipy=0.18.1 zope.interface=4.4.2 m2w64-gcc-fortran=5.3.0
(Anaconda3)> activate ccblade-devenv
(devenv)> pip install -r requirements.txt
(devenv)> ccblade
```

## Run Unit Tests

To check if installation was successful, run the unit tests

    $ python test/test_ccblade.py
    $ python test/test_gradients.py

For software issues please use <https://github.com/WISDEM/CCBlade/issues>.  For functionality and theory related questions and comments please use the NWTC forum for [Systems Engineering Software Questions](https://wind.nrel.gov/forum/wind/viewtopic.php?f=34&t=1002).

