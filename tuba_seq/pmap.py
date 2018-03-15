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
import numpy as np
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
    CHUNKS = 50*CPUs
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

class _mid_fastq_iter(object):
    def __init__(self, filename, seq_id, start, stop):
        self.filename = filename
        self.seq_id = seq_id
        self.start = start
        self.stop = stop

    def __iter__(self): 
        self.f = open(self.filename, 'rb')
        self.f.seek(self.start)
        # Find the beginning of the next FASTQ header
        for i, line in enumerate(self.f):
            if line.startswith(self.seq_id):
                break
        assert i <= 4, "Took more than 4 lines to find header"
        self.f.seek(self.f.tell() - len(line))
        return self
        
    def __next__(self):
        if self.f.tell() <= self.stop:
            header = self.f.readline()
            dna = self.f.readline()
            self.f.readline()
            qc = self.f.readline()
            return header, dna, qc
        else:
            self.f.close()
            raise StopIteration
        
    def __exit__(self):
        self.f.close()

def fastq_map_sum(in_fastq, out_filenames, func, CPUs=CPUs-1, temp_dir='tmp'):
    """Asynchronously processes an input fastq file.

Processes reads from a single FASTQ file by distributing the analysis work *and*
I/O. A thread is devoted to reading the input FASTQ, threads are devoted to 
writing to every output file, and a number of worker threads communicate with 
these I/O threads to process the reads. Job queues mediate thread interactions.

Inputs:
-------
in_fastq : Input FASTQ filename. Must be uncompressed. 

out_files : List of output filenames. Compression must be smart_open compliant.

func : func(fastq_read_iter, out_filenames) -> tuple of sum-able objects. 

"""
    with open(in_fastq, 'rb') as f:
        seq_id = f.readline().partition(':'.encode('ascii'))[0]
        file_length = f.seek(0, 2)
    
    start_positions = np.linspace(0, file_length, CHUNKS, False).round().astype(int)
    stop_positions = np.r_[start_positions[1:], file_length - 1]
    
    sample = os.path.basename(in_fastq).partition('.fastq')[0]
    # Make temp dir 
    Dir = temp_dir+sample
    Iter = [(_mid_fastq_iter(in_fastq, seq_id, start, stop), [os.path.join(Dir, head, str(start)+tail) for head, tail in map(os.path.split, out_filenames)]) for start, stop in zip(start_positions, stop_positions)]
    
    for head, tail in map(os.path.split, out_filenames):
        os.makedirs(os.path.join(Dir, head))
        os.makedirs(head, exist_ok=True)

    with multiprocessing.Pool(processes=CPUs) as P:
        rs = P.starmap_async(func, Iter, chunksize=1)
        outputs = rs.get()
    intermediate_files = list(zip(*[outfiles for it, outfiles in Iter]))
    for f, i_files in zip(out_filenames, intermediate_files):
        os.system("cat "+' '.join(i_files)+" >> "+f)
        list(map(os.remove, i_files))
        os.rmdir(os.path.join(Dir, os.path.split(f)[0]))
    os.rmdir(Dir)
    return list(map(sum, zip(*outputs))) if type(outputs[0]) == tuple else sum(outputs)

