#!/usr/bin/env python3
from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(  name='tuba_seq',
        description='Tools for the tuba-seq analysis pipeline',
        author='Christopher McFarland',
        author_email='christopherdmcfarland@gmail.com',
        packages=['tuba_seq'], 
        ext_modules=cythonize(Extension("*", ['tuba_seq/*.pyx'], 
            include_dirs=['seq-align/src', numpy.get_include()])),
        install_requires=[  
            'numpy',
            'scipy',
            'pandas',
            'matplotlib',
            'seaborn',
            'powerlaw',
            'biopython',
            'progressbar2'],
        version=1.0) 

