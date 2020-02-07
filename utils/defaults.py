import configparser
import os
import numpy as np

class Params(object):
    """
    Create the default parameters file
    """

    paramname = 'two-stream-tests.config'
    # VERY IMPORTANT: each testcase string MUST be 8 characters long!
    # ... fortran
    testcases = ['gaussian', 'betaTest', 'gammaTst']

    # Plotting specifics
    stdRatio = 16./9.
    stdSize = 5.
    dpi = 144
    framerate = 60

    def __init__(self, root):
        """
        Constructor

        Parameters
        ----------

        root : str
            The absolute path of where two-stream-tests.py script was ran from
        """

        # where are the parameter files going?
        self.root = root

        # The param file name
        self.paramfile = os.path.join(self.root, self.paramname)

        # for every run, there are constant parameters
        self.allrunparam = {'C':0.5, 'stdRatio':self.stdRatio, 'stdSize':self.stdSize,
                            'dpi':self.dpi, 'framerate':self.framerate}

        # all testcases in a single dictionary
        self.allTests = {}

        # create each parameter dictionary
        testFunctions = dict(zip(self.testcases, [self.gaussParams, self.betaTestParams,
                                                  self.gammaTstParams]))

        for testcase in self.testcases:
            testFunctions[testcase](testcase)

    def gaussParams(self, testcase):
        """
        Create the gaussian testcase parameters dictionary

        Parameters
        ----------

        testcase : str
            The testcase key
        """

        # the list of keys for the dictionary
        gausskeys = ['dm', 'nconv', 'sigma_m', 'mmax', 'v']

        # create the items
        nconv = 2

        # I will resolve gaussian to 5sigma in all cases. sigma_m gives the size of gaussian
        sigma_m = 2.5
        mmax = 20.

        # what is my velocity
        v = 1.

        # configparser converts everything to string. make csv for easy list making later
        dms = list(map(lambda x: 0.078125*np.power(2.,x), -1*np.arange(0, 5)))
        dms = ','.join(map('{:.14E}'.format, dms))

        self.allTests[testcase] = dict(zip(gausskeys, [dms, nconv, sigma_m, mmax, v]))

    def betaTestParams(self, testcase):
        """
        Create the betaTest testcase parameters dictionary

        Parameters
        ----------

        testcase : str
            The testcase key
        """

        # the list of keys for the dictionary
        betaTestkeys = ['dm', 'nconv', 'sigma_m', 'mmax', 'vmax', 'mstart', 'mstop']

        # convective turnover times assuming v=1
        nconv = 6

        # I will resolve gaussian to 5sigma in all cases. sigma_m gives the size of gaussian
        sigma_m = 1
        mmax = 40.

        # dm will be a sensible resolution so it doesn't take forever
        dm = '{:.10E}'.format(0.078125*np.power(2.,-3.))

        # the maximum velocity and where the quadratic fit will start and end (mass coord)
        vmax = 4.
        mstart = 10.
        mstop = 30.

        self.allTests[testcase] = dict(zip(betaTestkeys, [dm, nconv, sigma_m, mmax, vmax, mstart,
                                                          mstop]))

    def gammaTstParams(self, testcase):
        """
        Create the gammaTst testcase parameters dictionary

        Parameters
        ----------

        testcase : str
            The testcase key
        """

        # the list of keys for the dictionary
        gammaTstkeys = ['dm', 'nconv', 'sigma_m', 'mmax', 'vmax', 'mstart', 'mstop', 'gamma']

        # convective turnover times assuming v=1
        nconv = 6

        # I will resolve gaussian to 5sigma in all cases. sigma_m gives the size of gaussian
        sigma_m = 1
        mmax = 40.

        # dm will be a sensible resolution so it doesn't take forever
        dm = '{:.10E}'.format(0.078125*np.power(2.,-3.))

        # the maximum velocity and where the quadratic fit will start and end (mass coord)
        vmax = 4.
        mstart = 10.
        mstop = 30.

        # Gamma is a constant value of 0.5 the dm
        gamma = '{:.5E}'.format(0.05 * float(dm))

        self.allTests[testcase] = dict(zip(gammaTstkeys, [dm, nconv, sigma_m, mmax, vmax, mstart,
                                                          mstop, gamma]))


    def writeTestcases(self):
        """
        Create the parameters file
        """
        
        # I will use the configparser for writing these parameter files
        config = configparser.ConfigParser()
        config['DEFAULT'] = self.allrunparam

        # Each section is a testcase
        for testcase in self.testcases:
            config[testcase] = self.allTests[testcase]

        # write the param file
        with open(os.path.join(self.paramfile), 'w') as myfile:
            config.write(myfile)
