#!/usr/bin/env python

"""
This is a python DistUtils setup file for the
ISR Blackbox module.

$Id: setup.py 120 2021-07-06 21:09:40Z brideout $
"""

import os, sys

from distutils.core import setup, Extension

setup(name="ISR_Blackbox",
        version="1.0",
        description="Classes to emulate an ISR blackbox",
        author="Bill Rideout",
        author_email="brideout@mit.edu",
        py_modules=['ISR_blackbox'])
