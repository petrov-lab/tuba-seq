"""Parallel (multi-threaded) map function for python. 

Uses multiprocessing.Pool with error-resistant importing. There are two map
functions:

1) pmap(function, iterable) -> rapid fork-based multi-threaded map function.

2) low_memory_pmap(function, iterable) -> a more memory-efficient version
    intended for function calls that are individually long & memory-intensive.

"""

try:
    import multiprocessing
except ImportError:
    print("Cannot import 'multiprocessing' module. Parallelization not possible.")
    pmap = map
    low_memory_pmap = map
    CPUs = 1
finally:
    CPUs = multiprocessing.cpu_count()
    def pmap(func, Iter, processes=CPUs):
        with multiprocessing.Pool(processes=processes) as P:
            return P.map(func, Iter)

    def low_memory_pmap(func, Iter, processes=int(round(CPUs/2)), chunksize=1):
        with multiprocessing.Pool(processes=processes) as P:
            return [result for result in P.imap(func, Iter)]
        
