#!/usr/bin/env python
# encoding: utf-8

from setuptools import setup, find_packages
from numpy.distutils.core import setup, Extension



setup(
    name='CCBlade',
    version='1.1.1',
    description='Blade element momentum aerodynamics for wind turbines',
    author='S. Andrew Ning',
    author_email='andrew.ning@nrel.gov',
    packages=['ccblade'],
    install_requires=[
        'airfoilprep.py>=0.1',
        'numpy==1.12.1',
        'scipy==0.18.1',
        'zope.interface==4.4.2',
        'matplotlib==2.0.0'
    ],
    package_data={
        'ccblade': ['data/5MW_AFFiles/*']
    },
    # test_suite='test.test_ccblade.py',
    license='Apache License, Version 2.0',
    ext_modules=[Extension('ccblade._bem', ['ccblade/bem.f90'], extra_compile_args=['-O2'])],
    entry_points={
        'console_scripts': ['ccblade=ccblade.ccblade:run']
    },
    zip_safe=False
)
