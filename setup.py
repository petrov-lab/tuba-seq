#!/usr/bin/env python3
from distutils.core import setup
from Cython.Build import cythonize

seq_align_dir = '/home/chris/seq-align/src'

setup_kwargs = dict(name='tuba_seq',
                    description='Tools for tuba_seq analysis pipeline.',
                    author='Christopher McFarland',
                    author_email='christopherdmcfarland@gmail.com',
#                    url=#'https://www.python.org/sigs/distutils-sig/',
                    version='1.0')

setup(packages=['tuba_seq'], **setup_kwargs)
setup(ext_modules=cythonize('tuba_seq/*.pyx', include_path=[seq_align_dir]), **setup_kwargs)

#ext = Extension('nwalign_wrapper',
#                sources=['crap/nwalign_wrapper.pyx'], #, "/home/chris/seq-align/src/needleman_wunsch.h", "/home/chris/seq-align/src/alignment_scoring.h"],
#                include_dirs=[src_dir],
#                libraries=[src_dir],
#                library_dirs=[src_dir],
#                runtime_library_dirs=[src_dir])
#'/home/chris/seq-align/src'])

                # libraries=['seq-align'],
               # library_dirs=['/home/chris/seq-align/src'],
#setup(name='nwalign_wrapper', ext_modules=cythonize(ext))


#setup_kwargs = dict(name='crap',
#                    author='Christopher McFarland',
#                    author_email='christopherdmcfarland@gmail.com',
#                    version='1.0')
#


