"""Parallel (multi-threaded) map function for python. 

Uses multiprocessing.Pool with error-resistant importing. There are two map
functions:

1) pmap(function, iterable) -> rapid fork-based multi-threaded map function.

2) low_memory_pmap(function, iterable) -> a more memory-efficient version
    intended for function calls that are individually long & memory-intensive.

"""
import os
from warnings import warn
from pickle import PicklingError
from progressbar import ProgressBar, Bar, Percentage
from time import sleep 

try:
    import multiprocessing
except ImportError:
    print("Cannot import 'multiprocessing' module. Parallelization not possible.")
    pmap = map
    low_memory_pmap = map
    larger_iter_pmap = map
    CPUs = 1
finally:
    CPUs = multiprocessing.cpu_count()
    CHUNKS = 100*CPUs
    def pmap(func, Iter, processes=CPUs):
        with multiprocessing.Pool(processes=processes) as P:
            return P.map(func, Iter)

    def low_memory_pmap(func, Iter, processes=int(round(CPUs/2)), chunksize=1):
        with multiprocessing.Pool(processes=processes) as P:
            return [result for result in P.imap(func, Iter)]
        
    def large_iter_pmap(func, Iter, processes=CPUs, status_bar=True, nice=True, wait_interval=1):
        if nice:
            os.nice(10)
        try:
            with multiprocessing.Pool(processes=processes) as P:
                size = max(1, int(round(len(Iter)/CHUNKS)))
                rs = P.map_async(func, Iter, chunksize=size)
                maxval = rs._number_left
                bar = ProgressBar(maxval=maxval, widgets=[Bar('=', '[', ']'), ' ', Percentage()])
                while not rs.ready():
                    sleep(wait_interval)
                    bar.update(maxval - rs._number_left)
                bar.finish()
                return rs.get()

        except PicklingError:
            warn("Lambda functions cannot be Pickled for Parallelization. Using single Process.", RuntimeWarning)
            return list(map(func, Iter))

