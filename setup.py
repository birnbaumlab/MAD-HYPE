#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

## TODO: add README file
#with open('README.md') as f:
#    readme = f.read()
readme = ''

with open('LICENSE') as f:
    license = f.read()

setup(
    name='madhype',
    version='0.1.0',
    author='Patrick Holec, Joseph Berleant',
    author_email='pholec@mit.edu',
    url='https://github.com/OhEvolve/MAD-HYPE',
    description='Multicell Analytical Deconvolution for High Yield Paired-chain Evaluation',
    long_description=readme,
    license=license,
    classifiers=['Programming Language :: Python :: 2'],
    install_requires=[ ## TODO add version numbers for dependencies
        'numpy>=1.13.3',
        'scipy>=0.19.0,<1',
        'matplotlib>=2.0.0,<3',
        'openpyxl>=2.5.4',
        'tqdm>=4.23.4'
    ],
    packages=[ ## TODO replace with find_packages() call when madhype/ is cleaned up
        'madhype',
        'madhype.simulation',
        'madhype.analysis',
        'madhype.postprocessing'
    ]
)

