#!/usr/bin/env python3
from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(  name='tuba_seq',
        description='Tools for the tuba-seq analysis pipeline',
        url='https://github.com/petrov-lab/tuba-seq',
        author='Christopher McFarland',
        author_email='christopherdmcfarland@gmail.com',
        license='MIT',
        classifiers=[
            'Development Status :: 2 - Pre-Alpha',
            'Intended Audience :: Science/Research',
            'Natural Language :: English',
            'Operating System :: POSIX',
            'Topic :: Scientific/Engineering :: Bio-Informatics'],
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
            'progressbar2',
	    'rpy2',
	    'jinja2',
            'regex'],
        scripts=[
            'bin/preprocess.py',
            'bin/postprocess.py',
            'bin/PEAR.py',
            'bin/final_processing.py'],
        python_requires='~=3.2',
        version='2.1',
        ) 

