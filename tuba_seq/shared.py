from datetime import datetime

class logPrint(object):
    def __init__(self, filename=None, verbose=False, hbar=80*'-'+'\n'):
        import __main__ as main 
        import os
        self.start_time = datetime.now()
        self.program = os.path.basename(main.__file__).partition('.py')[0]
        self.filename = self.program+'.LOG' if filename is None else filename
        self.verbose = verbose
        print("Logging output to", self.filename) 
        self.f = open(self.filename, 'a')
        self.f.write(hbar)
        self.f.write("Output Summary of {self.program:}, executed on {self.start_time:%c}:\n".format(self=self))
        self.f.write(hbar)

    def __call__(self, line, print_line=False):
        if self.verbose or print_line:
            print(line)
        self.f.write(line+'\n')        

    def close(self):
        runtime = datetime.now() - self.start_time
        self('Runtime: {:}h:{:}m:{:}s'.format(*str(runtime).split('.')[0].split(':')), print_line=True)
        self.f.close()

class sampleMetaData(object):
    def __init__(self, metadata_df, verbose=False):
        self.df = metadata_df
