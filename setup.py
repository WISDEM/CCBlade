#!/usr/bin/env python
# encoding: utf-8

from numpy.distutils.core import setup, Extension
import os

os.environ['NPY_DISTUTILS_APPEND_FLAGS'] = '1'

setup(
    name='CCBlade',
    version='1.3.0',
    description='Blade element momentum aerodynamics for wind turbines',
    author='NREL WISDEM Team',
    author_email='systems.engineering@nrel.gov',
    #package_dir={'': 'src'},
    #py_modules=['ccblade'],
    package_data={'ccblade': []},
    packages=['ccblade'],
    install_requires=[
        "numpy",
        "openmdao>=3.4",
        "scipy",
    ],
    python_requires=">=3.7",
    # test_suite='test.test_ccblade.py',
    license='Apache License, Version 2.0',
    ext_modules=[Extension('ccblade._bem',
                           sources=[os.path.join('ccblade','src','bem.f90')],
                           extra_compile_args=['-O2','-fPIC','-std=c11'])],
    zip_safe=False
)
