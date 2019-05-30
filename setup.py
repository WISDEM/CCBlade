#!/usr/bin/env python
# encoding: utf-8

from setuptools import setup
from numpy.distutils.core import setup, Extension



setup(
    name='CCBlade',
    version='1.1.1',
    description='Blade element momentum aerodynamics for wind turbines',
    author='NREL WISDEM Team',
    author_email='systems.engineering@nrel.gov',
    package_dir={'': 'src'},
    py_modules=['ccblade'],
    package_data={'ccblade': []},
    packages=['ccblade'],
    install_requires=['airfoilprep>=0.1'],
    # test_suite='test.test_ccblade.py',
    license='Apache License, Version 2.0',
    ext_modules=[Extension('_bem', ['src/ccblade/bem.f90'], extra_compile_args=['-O2','-fPIC','-shared'], extra_link_args=['-shared'])],
    zip_safe=False
)
