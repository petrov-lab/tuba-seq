#!/usr/bin/env python3
from distutils.core import setup
from Cython.Build import cythonize

setup_kwargs = dict(name='tuba_seq',
                    description='Tools for tuba_seq analysis pipeline.',
                    author='Christopher McFarland',
                    author_email='christopherdmcfarland@gmail.com',
#                    url=#'https://www.python.org/sigs/distutils-sig/',
                    version='1.0')

setup(packages=['tuba_seq'], **setup_kwargs)
setup(ext_modules=cythonize('tuba_seq/*.pyx'), **setup_kwargs)
