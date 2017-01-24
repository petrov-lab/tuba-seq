#!/usr/bin/env python3
from distutils.core import setup
from Cython.Build import cythonize
#from sys import argv

scripts = ['fastq.pyx', 'nwalign.pyx']

for script in scripts:
    short = script.split('.')[0]
    setup(name=short, ext_modules=cythonize(script))
