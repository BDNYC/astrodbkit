#! /usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from setuptools import setup, find_packages
    setup
except ImportError:
    from distutils.core import setup
    setup

from codecs import open
from os import path

from astrodbkit import __version__

setup(
    name='astrodbkit',
    version=__version__,
    description='Astronomy database management using SQL and Python',
    url='https://github.com/BDNYC/astrodbkit.git',
    author='Joe Filippazzo',
    author_email='bdnyc.group@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
    ],
    keywords='astrophysics',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=['numpy','astropy','matplotlib'],

)