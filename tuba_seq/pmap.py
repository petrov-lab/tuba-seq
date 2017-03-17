from multiprocessing import Pool, cpu_count

def pmap(func, Iter, processes=cpu_count() - 1):
    with Pool(processes=processes) as P:
        return P.map(func, Iter)

