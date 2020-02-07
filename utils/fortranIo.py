import subprocess
import os
import numpy as np

# lists of stuff
filetypes = ['shell', 'central']
centralvars = ['Xup', 'Xdown', 'gamma', 'dm']
shellvars = ['v', 'rho', 'r']
filevars = ['dt', 'nmax', 'steps', 'saveevery', 'testcase']
fileline = ['double precision, parameter :: dt =', 'integer, parameter :: nmax =',
            'integer, parameter :: steps =', 'integer, parameter :: saveevery =',
            'character*8, parameter :: testcase =']
lineform = ['{:.10f}d0', '{:d}', '{:d}', '{:d}', '\'{:s}\'']


class FFiles(object):
    """
    Write input files for the fortran simulation of a testcase. Read output files of the fortran
    simulation
    """

    strFFormat = '{:0.12E}'

    def __init__(self, root):
        """
        Constructor

        Parameters
        ----------

        root : str
            The absolute path of where two-stream-tests.py script was ran from
        """

        self.root = os.path.abspath(root)

    def readFFile(self, filename):
        """
        Read a particular fortran output file

        Parameters
        ----------

        filename : str
            The absolute path to the file to be read

        Returns
        -------

        Xup : np.ndarray
            The mass fraction of species X from the simulation in the upstream

        notXup : np.ndarray
            The mass fraction of species (1-X) from the simulation in the upstream

        Xdown : np.ndarray
            The mass fraction of species X from the simulation in the downstream

        notXdown : np.ndarray
            The mass fraction of species (1-X) from the simulation in the downstream
        """

        # since it is formatted simply, it is actually easiest to just use np
        data = np.genfromtxt(filename)

        Xup = np.copy(data[:,0])
        notXup = np.copy(data[:,1])
        Xdown = np.copy(data[:,2])
        notXdown = np.copy(data[:,3])

        return Xup, notXup, Xdown, notXdown

    def relativeRoot(self, path):
        """
        Make the path a relative path with respect to self.root

        Parameters
        ----------

        path : str
            The absolute path

        Returns
        -------

        relpath : str
            The input path that is transformed to a relative path with respect to self.root
        """

        # Confirm that the path can be relative to root, then return relative path
        realpath = os.path.abspath(path)
        if os.path.commonpath([self.root, realpath]) is not '/':
            return os.path.relpath(realpath, self.root)

    def writeFFiles(self, inputFFile, runVars, runConsts):
        """
        Write the shell and central input files. Note that the order of variables in columns are
        based on the `shellvars` and `centralvars` lists

        Parameters
        ----------

        filename : str
            The name of the file, does not include extension

        runVars : dict
            A dictionary that contains all of the relevant run-time data, both shell and central
            quantities

        runConsts : dict
           A dictionary that contains of all of the relevant run-time constants that are hard-coded
           into the fortran code
        """

        # create a list of the list of variables for EACH filetype
        filetype_vars = [shellvars, centralvars]

        for i, filetype in enumerate(filetypes):

            # what is the file being written?
            basename = os.path.basename(inputFFile)
            directory = os.path.dirname(inputFFile)
            basename = basename.split('.')[0]
            
            inputFFile = os.path.join(directory, basename + '.{:s}'.format(filetype))

            # open the file
            ffile = open(inputFFile, 'w')
            
            # This is a list of arrays of all the variable types
            writevars = list(map(runVars.get, filetype_vars[i]))

            # in each array, I want to grab one element and write them all to a line
            for j in range(len(writevars[0])):

                # write a line
                if j != len(runVars.get(filetype_vars[i][0]))-1:
                    strings = list(map(self.strFFormat.format, map(lambda x: x[j], writevars)))
                    ffile.write(' '.join(strings) + '\n')
                else:
                    strings = list(map(self.strFFormat.format, map(lambda x: x[j], writevars)))
                    ffile.write(' '.join(strings))
            ffile.close()
            
            # Now we modify the compiled file
            self.modifyCFile(inputFFile.replace('.{:s}'.format(filetype), '.f'), runConsts)

    def modifyCFile(self, compileFFile, runConsts):
        """
        Modify the appropriate parameters within the fortran file being compiled

        Parameters
        ----------

        compileFFile : str
            The name of the fortran file that we are going to compile

        runConsts : dict
           A dictionary that contains of all of the relevant run-time constants that are hard-coded
           into the fortran code
        """

        cfile = open(compileFFile, 'r')
        lines = cfile.readlines()

        for linenum, line in enumerate(lines):
            for i, val in enumerate(map(line.find, fileline)):
                if val != -1:

                    # split line on an = sign. we need to change fortran float or int portion
                    fnum = lineform[i].format(runConsts[filevars[i]])

                    split = line.split('=')
                    lines[linenum] = split[0] + '= ' + fnum + '\n'

        cfile.close()

        # we have the new lines, we will rewrite the file
        modify = open(compileFFile, 'w')
        modify.writelines(lines)
        modify.close()

        return True
    
    def runFFile(self, compileFFile):
        """
        Use subprocess to run the gfortran code

        Parameters
        ----------

        compileFFile : str
            The name of the fortran file that we are going to compile
        """

        # first compile the file
        binaryFFile = self._compile(compileFFile)

        try:
            subprocess.run(['./{:s}'.format(binaryFFile)],
                           stdout=subprocess.PIPE, check=True)

        except:
            print('The running of the compiled file, {:s}, failed'.format(compileFFile))
            raise 

        # Success! cd back to root
        print('The compiled fortran code succeeded in its run!')
        os.chdir(self.root)

        return True

    def _compile(self, compileFFile):
        """
        Use gfortran to compile the file to be run

        Parameters
        ----------

        compileFFile : str
            The name of the fortran file that we are going to compile
        """

        # this can fail so we try
        try:

            compileFFile = self.relativeRoot(compileFFile)
            
            # for fortran to read the files consistently, we MUST cd into the run directory
            os.chdir(os.path.dirname(compileFFile))

            # Now we use relative to the run directory
            compileFFile = os.path.basename(compileFFile)
            binaryFFile = compileFFile.split('.')[0]+'.a'
            
            subprocess.run(['gfortran', compileFFile, "-o", binaryFFile], stdout=subprocess.PIPE,
                           check=True)

        except:
            print('The compilation of the file, {:s}, failed'.format(compileFFile))
            raise

        return binaryFFile
