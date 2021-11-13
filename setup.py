#!/usr/bin/env python

# setup script for estp
#
# Copyright (C) 2015 Donghan Lee
# http://www.bionmr.org
# License: BSD (See LICENSE for full license)
#
# $Date$

from distutils.core import setup

setup(
    name='estp',
    version='1.0',
    author = 'Donghan Lee and Marta G. Carneiro',
    author_email = 'dlee@bionmr.org',
    packages=['estp'],
    license = 'New BSD License',
    url = 'http://www.bionmr.org', 
    description = 'A program for CEST and DEST experiment',
    long_description = open('README').read(),
    requires = ['numpy','scipy'],
    classifiers = [
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux'
    ]
)
