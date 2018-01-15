"""Shared utilities for this package.

Class logPrint both 'logs' and 'prints' (when verbosity insists) any output
string. Additional info of every command-line script (arguments, runtime) are
also logged. 

"""
from datetime import datetime
import atexit, warnings, functools, os

def smart_open(filename, mode='rb', makedirs=False):
    if makedirs:
        directory = os.path.dirname(filename)
        if directory:
            os.makedirs(directory, exist_ok=True)
    if filename[-5:] == '.gzip' or filename[-3:] == '.gz':
        from gzip import open
    elif filename[-4:] == '.bz2':
        from bz2 import open
    elif filename[-5:] == '.lzma':
        from lzma import open
    else:
        from builtins import open
    return open(filename, mode)

class logPrint(object):
    def line_break(self):
        self.f.write(80*'-'+'\n')

    def __call__(self, line, print_line=False, header=False):
        if self.verbose or print_line:
            print(line)
        if header:
            self.f.write((len(line)+4)*"#"+'\n')
            self.f.write('# '+str(line)+' #\n')        
            self.f.write((len(line)+4)*"#"+'\n')
        else:
            self.f.write(str(line)+'\n')        

    def close_logPrint(self):
        runtime = datetime.now() - self.start_time
        self('Runtime: {:}'.format(str(runtime).split('.')[0]))
        self.line_break()
        self.f.close()
    
    def __init__(self, input_args, filename=None):
        import __main__ as main 
        import os
        self.start_time = datetime.now()
        self.program = os.path.basename(main.__file__).partition('.py')[0]
        self.filename = self.program+'.LOG' if filename is None else filename
        args_dict = input_args.__dict__.copy()
        self.verbose = args_dict.pop('verbose', False) 
        print("Logging output to", self.filename) 
        self.f = open(self.filename, 'a')
        self.f.write('\n')
        self.line_break()
        self.f.write("Output Summary of {0.program}, executed at {0.start_time:%c} with the following input arguments:\n".format(self))
        self.line_break()
        for arg, val in args_dict.items():
            self.f.write("{:}: {:}\n".format(arg, val))
        self.line_break()
        atexit.register(self.close_logPrint)

class sampleMetaData(object):
    def __init__(self, metadata_df, verbose=False):
        self.df = metadata_df

def ignore_warning(warning, count=None):
    """Courtesy of https://gist.github.com/WoLpH/ebebe5d693fe6d0ad1c8"""
    def _ignore_warning(function):
        @functools.wraps(function)
        def __ignore_warning(*args, **kwargs):
            with warnings.catch_warnings(record=True) as ws:
                # Catch all warnings of this type
                warnings.simplefilter('always', warning)
                # Execute the function
                result = function(*args, **kwargs)

            # If we are looking for a specific amount of
            # warnings, re-send all extra warnings
            if count is not None:
                for w in ws[count:]:
                    warnings.showwarning(
                        message=w.message,
                        category=w.category,
                        filename=w.filename,
                        lineno=w.lineno,
                        file=w.file,
                        line=w.line,
                    )

            return result
        return __ignore_warning
    return _ignore_warning

