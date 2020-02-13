import configparser
import os
import shutil
import time
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from utils import fortranIo, gaussian, v_quad

class Testcase(object):
    """
    A base class for all testcases
    """

    advFile = 'two-stream-mixing.f'
    programFile = 'program.f'

    def __init__(self, testcase, paramfile, root, plots, movies):
        """
        Construct the testcase base parameters

        Parameters
        ----------

        testcase : str
            The testcase we are going to compute

        paramfile : str
            The absolute path to the parameters file

        root : str
            The absolute path to where the two-stream-tests.py was called from

        plots : bool
            A boolean as to whether or not we make plots by default

        movies : bool
            A boolean as to whether or not we will create movies with ffmpeg
        """

        # A bunch of strings/bools to hold
        self.mytestcase = testcase
        self.root = root
        self.plots = plots
        self.movies = movies
        self.rundir = os.path.join(self.root, self.mytestcase)
        self.runFile = os.path.join(self.rundir, self.mytestcase + '.f')

        # read in the parameters file and grab appropriate parameters
        configs = configparser.ConfigParser()
        configs.read(paramfile)
        self.testParams = configs[self.mytestcase]

        # grab the defaults
        self.courant = configs['DEFAULT'].getfloat('c')
        self.stdRatio = configs['DEFAULT'].getfloat('stdRatio')
        self.stdSize = configs['DEFAULT'].getfloat('stdSize')
        self.dpi = configs['DEFAULT'].getfloat('dpi')
        self.framerate = configs['DEFAULT'].getint('framerate')

        # Now we can create the run directory and the default fortran file
        if os.path.isdir(self.rundir):
            self.clean()

        # make rundir and composite file
        os.mkdir(self.rundir)

        with open(os.path.join(self.root, self.advFile), 'r') as afile:
            advlines = afile.readlines()

        with open(os.path.join(self.root, self.programFile), 'r') as afile:
            prlines = afile.readlines()

        with open(self.runFile, 'w') as afile:
            self.defaultLines = []
            self.defaultLines.extend(prlines)
            self.defaultLines.extend(advlines)
            afile.writelines(self.defaultLines)

    def clean(self):
        """
        Clean the run directory of the testcase if it already exists
        """

        shutil.rmtree(self.rundir)

    def createVars(self):
        """
        Create all of the shell and central quantities that will be a part of the self.runVars
        dictionary. These variables get written out to files for the fortran code to read-in
        """

        # get the dm for the testcase
        dm = self.sim_dm()

        # other params we need and are common to everything
        mmax = self.testParams.getfloat('mmax')

        self.mass_shell = np.arange(0,int(np.round(mmax / dm))+1) * dm
        self.mass_central = (self.mass_shell + 0.5*dm)[0:-1]

        # Create the shell quantities, radius, rho and v
        radius = self.sim_radius()
        rho = self.sim_rho()
        v = self.sim_v()

        # create the gaussian mass fraction arrays, gamma, and dm_arrays
        Xup, Xdown = self.sim_X()
        gamma = self.sim_gamma()

        # This is always the same for any testcase
        dm_array = np.ones(self.mass_central.shape) * dm

        # Return as tuple for instantiating a dict
        return zip(fortranIo.shellvars + fortranIo.centralvars,
                   [v, rho, radius, Xup, Xdown, gamma, dm_array])

    def createRunConsts(self):
        """
        Create the run constants that are hard coded into the compiled fortran code
        """

        # We need to know how many convective turnovers we are doing for the testcase; always in
        # params file
        nconv = self.testParams.getfloat('nconv')

        # the maximum time step is based on C, dm and the flux (specific to a testcase)
        dt = self.dt_max()

        # we need the number of central mass zones
        nmax = len(self.mass_central)

        # we need the number of steps we are going to take. This is based on nconv, nmax and C
        steps = int(np.round(nconv * nmax / self.courant))

        # we need to know the frequency to which we are outputting.
        saveevery = self.save_every()

        # create the dictionary for later use
        return zip(fortranIo.filevars, [dt, nmax, steps, saveevery, self.mytestcase])

    def run_simulation(self):
        """
        Run the testcase simulation
        """

        # We create the run variables based on what the testcase is.
        self.runVars = dict(self.createVars())

        # Now we get the run variables based on what the testcase is.
        self.runConsts = dict(self.createRunConsts())

        print('Working on {:s} testcase. We are going to run a simulation with {:d} mass zones'
                .format(self.mytestcase, self.runConsts['nmax']))

        # Everything is set, now we create the shell and central files
        self.fortran = fortranIo.FFiles(self.root)
        self.fortran.writeFFiles(self.runFile, self.runVars, self.runConsts)

        # make a txt directory for the output
        os.mkdir(os.path.join(self.rundir, 'output'))

        # make a png directory
        os.mkdir(os.path.join(self.rundir, 'png'))

        # now we compile and run the fortran code
        timeit = time.time()
        self.fortran.runFFile(self.runFile)

        print('That took {:.2f} minutes'.format((time.time() - timeit)/60.))
        print('')

        # find all of the '*.txt' files
        txtfiles = []
        for afile in os.listdir(self.rundir):
            if afile.endswith('.txt'):
                txtfiles.append(os.path.join(self.rundir, afile))

        if self.plots:
            print('Believe it or not, plotting takes time. We need to plot {:d} figures'.format(len(txtfiles)))
            timeit = time.time()
            self.create_plots(txtfiles)

            print('That took {:.2f} minutes'.format((time.time() - timeit)/60.))
            print('')

        # find all of the '*.png' files
        pngfiles = []
        for afile in txtfiles:
            pngfiles.append(afile.replace('.txt', '.png'))

        # move all txt and png files to the folders
        for atxt, apng in zip(txtfiles, pngfiles):
            shutil.move(os.path.join(self.rundir, atxt), os.path.join(self.rundir, 'output'))

            if self.plots:
                shutil.move(os.path.join(self.rundir, apng), os.path.join(self.rundir, 'png'))

        if self.movies:
            print('Creating the movies out of the png files')
            self.create_movies()
            print('')

    def create_movies(self):
        """
        Create movies from the png figures made in self.create_plots
        """

        # it is best to just cd into the directory and run it from there
        os.chdir(os.path.join(self.rundir, 'png'))

        # So, within myfolder there are a bunch of *.png files, make these into a movie
        try:
            subprocess.run(['ffmpeg', '-framerate', '{:d}'.format(self.framerate), '-pattern_type',
                            'glob', '-i', '*.png', '-c:v', 'libx264', '-pix_fmt', 'yuv420p',
                            '{:s}.mp4'.format(self.mytestcase)])
        except subprocess.CalledProcessError as error:
            print('FFmpeg was unable to create a movie')
            raise

        # cd back to root
        os.chdir(self.root)

        # Move the movie to rundir
        shutil.move(os.path.join(self.rundir, 'png', '{:s}.mp4'.format(self.mytestcase)),
                    os.path.join(self.rundir, '{:s}.mp4'.format(self.mytestcase)))

    def create_plots(self, txtfiles):
        """
        Create plots for all of the fortran output data

        Parameters
        ----------

        txtfiles : list
            Contains the absolute paths of all of the output txt files
        """

        # read in every file, create arrays, send information needed to self._plot
        ifig = 0

        for output in txtfiles:

            # get quantities and the output number
            Xup, notXup, Xdown, notXdown = self.fortran.readFFile(output)
            basename = os.path.basename(output)
            num = int((basename.split('.')[0]).split('-')[-1])

            # so I have all of the X's from output file.
            self._plot(Xup, Xdown, ifig, num, len(txtfiles))
            ifig += 1

    ##--------------------+
    ## Overridden methods
    ##--------------------+

    def sim_dm(self):
        """
        Determine the constant dm spacing for the particular simulation

        Returns
        -------

        dm : float
            The spacing of mass cells to be used
        """

        # specified by a child class
        return None

    def sim_v(self):
        """
        The velocity field of the simulation

        Returns
        -------

        v : np.ndarray
            The velocity field defined on the shell
        """

        # specified by a child class
        return None

    def sim_rho(self):
        """
        The density field of the simulation

        Returns
        -------

        rho : np.ndarray
            The density field defined on the shell
        """

        # specified by a child class
        return None

    def sim_radius(self):
        """
        The radius of the shells of the cells in the simulation

        Returns
        -------

        r : np.ndarray
            The radius of the shell
        """

        # specified by a child class
        return None

    def sim_X(self):
        """
        The initial mass fraction distribution of the X species as a function of the central
        mass coordinates within the simulation

        Returns
        -------

        Xup : np.ndarray
            The species X in the upstream

        Xdown : np.ndarray
            The species X in the downstream
        """

        # specified by a child class
        return None

    def sim_gamma(self):
        """
        The gamma coefficient defined at the central mass coordinates within the simulation

        Returns
        -------

        gamma : np.ndarray
            The gamma coefficient
        """

        # specified by a child class
        return None

    def dt_max(self):
        """
        Determine the constant, maximum timestep that the simulation can take with the corresponding
        Courant number

        Returns
        -------

        dt : float
            The maximum timestep
        """

        # specified by a child class
        return None

    def save_every(self):
        """
        Determine how many timesteps are undertaken before a single dump of data is outputted

        Returns
        -------

        saveevery : int
            save every x timesteps
        """

        # specified by a child class
        return None

    def _plot(self, Xup, Xdown, ifig, num, numfiles):
        """
        Make and save a single plot for the testcase

        Parameters
        ----------

        Xup : np.ndarray
            The mass fraction of species X from the simulation in the upstream

        Xdown : np.ndarray
            The mass fraction of species X from the simulation in the downstream

        ifig : int
            The figure number to make

        num : int
            The output txt file number that is to be plotted

        numfiles : list
            The number of files that are being plotted
        """

class Gaussian(Testcase):
    """
    The gaussian testcase. This will test the importance of the resolution of the grid to the
    accuracy of the numerical solution of the advection equation. The initial state is a gaussian
    formed at the bottom of the upstream.
    """

    def __init__(self, testcase, paramfile, root, plots, movies):

        # Run parent initialization
        super().__init__(testcase, paramfile, root, plots, movies)

    # run_simulation
    def run_simulation(self):

        # For gaussian case, we have multiple dm's to run. 
        # we have our gaussian folder but we do multiple compile and compute runs. Make directories
        # to store those files
        folder_str = os.path.join(self.rundir, 'dm_{:.3f}')
        dm_floats = list(map(float, self.testParams.get('dm').split(',')))
        self.folders = list(map(folder_str.format, dm_floats))

        # store the inf norm (max abs error)
        errors_inf = []
        errors_l2 = []
        errors_l2_normalized = []

        # using multiple dm's, we just do our run_simulation multiple times by only updating what
        # self.rundir, self.runFile are each time
        for self.index, folder in enumerate(self.folders):

            # make the directory and update self.rundir and self.runFile
            oldfile = self.runFile
            os.mkdir(folder)
            self.rundir = folder
            self.runFile = os.path.join(self.rundir, os.path.basename(oldfile))

            # copy the fortran file to the appropriate directory
            shutil.copyfile(oldfile, self.runFile)

            # if we are on our first index (0), we delete the oldfile
            if self.index == 0:
                os.remove(oldfile)

            # now we can call super run_simulation
            super().run_simulation()

            # find all of the '*.txt' files
            txtfiles = []
            for afile in os.listdir(os.path.join(self.rundir, 'output')):
                if afile.endswith('.txt'):
                    txtfiles.append(os.path.join(self.rundir, afile))
            numfile = len(txtfiles) - 1

            # only take the Xup mass fraction array, we didn't advect far
            Xup = self.fortran.readFFile(
                    os.path.join(self.rundir, 'output', self.mytestcase+'-{:06d}.txt'.format(numfile)))[0]

            # now we know what the gaussian SHOULD be, it is just shifted by 1 sigma in mass coord
            shifted_gauss = gaussian(self.mass_central, self.mass_shell, self.testParams.getfloat('sigma_m'), 10., limit=False)

            # genfromtxt gives nans for real small numbers. Set these to zero
            Xup = np.nan_to_num(Xup)

            # different norms for the error
            diff = shifted_gauss - Xup
            errors_inf.append(np.max(np.abs(diff)))
            errors_l2.append(np.sqrt(np.sum(diff**2)))
            errors_l2_normalized.append(np.sqrt((1. / float(len(self.mass_central))) * np.sum(diff**2)))

        # using our l2_norms, plot_the simError
        self.plot_simError(errors_l2_normalized, dm_floats)
        

    def plot_simError(self, error, dm_floats):
        """
        Plot the error convergence for the simulations

        Parameters
        ----------
        error : list
            The errors for the particular dm_floats

        dm_floats : list
            List of the simulation dms used
        """

        # create figure
        ifig = 100000
        fig = plt.figure(ifig, figsize=(self.stdRatio*self.stdSize, self.stdSize), dpi=self.dpi)
        ax = fig.add_subplot(111)

        # plot the l2_norm vs dm_floats in log-log space. Add a ~delta_dm^2
        ax.scatter(dm_floats, error, color='k', marker='x', label='2nd-Order Solver With Minmod Limiter')

        dm_line = np.linspace(min(dm_floats), max(dm_floats), 1000.)
        second_order = 0.9*dm_line**2
        first_order = (0.9*dm_line[-1])*dm_line
        ax.plot(dm_line, second_order, color='k', ls='--', label=r'$\mathcal{0}(\delta m^{2})$')
        ax.plot(dm_line, first_order, color='r', ls='--', label=r'$\mathcal{0}(\delta m)$')

        ax.set_yscale('log')
        ax.set_xscale('log')

        # Plot details
        ax.set_xlabel('dm')
        ax.set_ylabel('L2 norm of errors')
        ax.set_ylim([1e-2,1e-5])
        ax.set_xlim([1e-3, 1e-1])
        ax.set_title('Convergence of Gaussian Advection')
        ax.legend()

        # invert both axes
        ax.invert_yaxis()
        ax.invert_xaxis()

        fig.savefig(os.path.join(self.root, self.mytestcase, 'L2-norm.png'), dpi=self.dpi)
        plt.close()
        
    def _plot(self, Xup, Xdown, ifig, num, numfiles):

        # create figure
        fig = plt.figure(ifig, figsize=(self.stdRatio*self.stdSize, self.stdSize), dpi=self.dpi)
        ax = fig.add_subplot(111)

        # plot Xup and Xdown as a function of self.mass_central
        ax.plot(self.mass_central, Xup, label='Upstream')
        ax.plot(self.mass_central, Xdown, label='Downstream')

        # if on last number
        if num == numfiles - 1:
            og_gauss = gaussian(self.mass_central, self.mass_shell, self.testParams.getfloat('sigma_m'), 10., limit=False)
            ax.plot(self.mass_central, og_gauss, color='k', ls='--', label='Gaussian at t=0')

        # Plot details
        ax.set_xlim([self.mass_central.min(), self.mass_central.max()])
        ax.set_xlabel('M')
        ax.set_ylim([0, 1])
        ax.set_ylabel('X')

        # add some quantities of interest in title
        dm = self.sim_dm()
        step = int(num * (0.1/dm))

        ax.set_title('dm={:.3f}, C={:.1f}, step={:d}'.format(dm, self.courant, step))

        ax.legend(frameon=False)

        fig.savefig(os.path.join(self.rundir, self.mytestcase + '-{:06d}'.format(num)),
                    dpi=self.dpi)
        plt.close()

    ##--------------------+
    ## Overridden methods
    ##--------------------+

    def sim_dm(self):

        # based on our index, return the appropriate float
        dm_floats = list(map(float, self.testParams.get('dm').split(',')))
        return dm_floats[self.index]

    def sim_v(self):
        
        return np.ones(self.mass_shell.shape) * self.testParams.getfloat('v')

    def sim_rho(self):

        return np.power(2.*np.pi,-1) * np.ones(self.mass_shell.shape)

    def sim_radius(self):

        return np.ones(self.mass_shell.shape)

    def sim_X(self):

        Xup = gaussian(self.mass_central, self.mass_shell, self.testParams.getfloat('sigma_m'),
                       10., limit=False)
        Xdown = gaussian(self.mass_central, self.mass_shell, self.testParams.getfloat('sigma_m'),
                         10., limit=False)

        return Xup, Xdown

    def sim_gamma(self):

        return np.zeros(self.mass_central.shape)

    def dt_max(self):

        # 2 * pi * r^2 * rho is a flux of 1 everywhere in this testcase
        return self.sim_dm() * self.courant / (1. * np.mean(self.sim_v()))

    def save_every(self):
        
        # We want to have approximately one output for every 0.1 of a mass coordinate output
        return int(np.round(0.1 / (self.courant * self.sim_dm())))

class BetaTest(Testcase):
    """
    The beta testcase. This simulation will use a non-constant velocity field to enforce horizontal
    mixing between the two streams. The initial state is a gaussian formed at the bottom of the
    upstream.
    """

    def __init__(self, testcase, paramfile, root, plots, movies):

        # call super init
        super().__init__(testcase, paramfile, root, plots, movies)

    ##--------------------+
    ## Overridden methods
    ##--------------------+

    def sim_dm(self):

        # simply
        return self.testParams.getfloat('dm')

    def sim_v(self):

        # Here we will construct a quadratic function defined with 3 points
        vmax = self.testParams.getfloat('vmax')
        mstart = self.testParams.getfloat('mstart')
        mstop = self.testParams.getfloat('mstop')

        # indices where I apply quadratic formula
        starti = np.argmin(abs(self.mass_shell - mstart))
        stopi = np.argmin(abs(self.mass_shell - mstop))

        v = np.ones(self.mass_shell.shape)
        vquad = v_quad(self.mass_shell, vmax, starti, stopi)
        v[starti:stopi+1] = vquad

        return v

    def sim_rho(self):

        return np.power(2.*np.pi,-1) * np.ones(self.mass_shell.shape)

    def sim_radius(self):

        return np.ones(self.mass_shell.shape)

    def sim_X(self):

        # 5 sigma is the cutoff, this is where the mean is located
        # note, m is the central mass coordinate, where X is defined
        Xup = gaussian(self.mass_central, self.mass_shell, self.testParams.getfloat('sigma_m'), 5.,
                       limit=True, extent=5.)
        Xdown = np.zeros(self.mass_central.shape)

        return Xup, Xdown

    def sim_gamma(self):

        return np.zeros(self.mass_central.shape)

    def dt_max(self):

        # 2 * pi * r^2 * rho * v plus beta is the flux in this case
        # note 2 * pi * r^2 * rho is 1 always!
        return self.sim_dm() * self.courant / np.max((self.runVars['v'][0:-1] +
                                                      np.abs(np.diff(self.runVars['v']))))
                                
    def save_every(self):
        
        # We want to have approximately one output for every 0.4 of a mass coordinate output.
        # This assumes v=1 everywhere which is not the case here but is good enough
        return int(np.round(0.4 / (self.courant * self.sim_dm())))
    
    def _plot(self, Xup, Xdown, ifig, num, numfiles):
        """
        Make a single plot
        """

        # create figure
        fig = plt.figure(ifig, figsize=(self.stdRatio*self.stdSize, self.stdSize), dpi=self.dpi)
        ax = fig.add_subplot(111)

        # plot Xup and Xdown as a function of self.mass_central
        ln1 = ax.plot(self.mass_central, Xup, label='Upstream')
        ln2 = ax.plot(self.mass_central, Xdown, label='Downstream')

        # We also want to plot the velocity field
        ax2 = ax.twinx()
        ln3 = ax2.plot(self.mass_shell, self.runVars['v'], label='v', color='k', ls=':')

        lns = ln1+ln2+ln3
        labs = [l.get_label() for l in lns]

        # Plot details
        ax.set_xlim([self.mass_central.min(), self.mass_central.max()])
        ax.set_xlabel('M')
        ax.set_ylim([0, 1])
        ax.set_ylabel('X')

        ax2.set_ylim([self.runVars['v'].min(), self.runVars['v'].max()])
        ax2.set_ylabel('v')

        # add some quantities of interest in title
        dm = self.sim_dm()
        step = int(num * self.runConsts['saveevery'])

        ax.set_title('dm={:.3f}, C={:.1f}, step={:d}'.format(dm, self.courant, step))

        ax.legend(lns, labs, frameon=False)

        fig.savefig(os.path.join(self.rundir, self.mytestcase + '-{:06d}'.format(num)),
                    dpi=self.dpi)
        plt.close()


class GammaTest(Testcase):
    """
    The gamma testcase. This will test the enforced horizontal mixing, beta, due to the non-constant
    velocity field as well as the additional horizontal mixing due to gamma. Gamma is constant
    throughout the entire domain. The initial state is a gaussian formed at the bottom of the
    upstream. 
    """

    def __init__(self, testcase, paramfile, root, plots, movies):

        # call super init!
        super().__init__(testcase, paramfile, root, plots, movies)

    ##--------------------+
    ## Overridden methods
    ##--------------------+

    def sim_dm(self):

        # simply
        return self.testParams.getfloat('dm')

    def sim_v(self):

        # Here we will construct a quadratic function defined with 3 points
        vmax = self.testParams.getfloat('vmax')
        mstart = self.testParams.getfloat('mstart')
        mstop = self.testParams.getfloat('mstop')

        # indices where I apply quadratic formula
        starti = np.argmin(abs(self.mass_shell - mstart))
        stopi = np.argmin(abs(self.mass_shell - mstop))

        v = np.ones(self.mass_shell.shape)
        vquad = v_quad(self.mass_shell, vmax, starti, stopi)
        v[starti:stopi+1] = vquad

        return v

    def sim_rho(self):

        return np.power(2.*np.pi,-1) * np.ones(self.mass_shell.shape)

    def sim_radius(self):

        return np.ones(self.mass_shell.shape)

    def sim_X(self):

        # 5 sigma is the cutoff, this is where the mean is located
        # note, m is the central mass coordinate, where X is defined
        Xup = gaussian(self.mass_central, self.mass_shell, self.testParams.getfloat('sigma_m'), 5.,
                       limit=True, extent=5.)
        Xdown = np.zeros(self.mass_central.shape)

        return Xup, Xdown

    def sim_gamma(self):

        # we have a factor for gamma
        gamma_factor = self.testParams.getfloat('gamma')
        return np.ones(self.mass_central.shape) * gamma_factor

    def dt_max(self):

        # 2 * pi * r^2 * rho * v plus beta AND gamma is the flux in this case
        # note 2 * pi * r^2 * rho is 1 always!
        return self.sim_dm() * self.courant / np.max((self.runVars['v'][0:-1] +
                                                      np.abs(np.diff(self.runVars['v'])) +
                                                      self.runVars['gamma']))
                                
    def save_every(self):
        
        # We want to have approximately one output for every 0.4 of a mass coordinate output.
        # This assumes v=1 everywhere which is not the case here but is good enough
        return int(np.round(0.4 / (self.courant * self.sim_dm())))

    def _plot(self, Xup, Xdown, ifig, num, numfiles):

        # create figure
        fig = plt.figure(ifig, figsize=(self.stdRatio*self.stdSize, self.stdSize), dpi=self.dpi)
        ax = fig.add_subplot(111)

        # plot Xup and Xdown as a function of self.mass_central
        ln1 = ax.plot(self.mass_central, Xup, label='Upstream')
        ln2 = ax.plot(self.mass_central, Xdown, label='Downstream')

        # We also want to plot the velocity field
        ax2 = ax.twinx()
        ln3 = ax2.plot(self.mass_shell, self.runVars['v'], label='v', color='k', ls=':')

        lns = ln1+ln2+ln3
        labs = [l.get_label() for l in lns]

        # Plot details
        ax.set_xlim([self.mass_central.min(), self.mass_central.max()])
        ax.set_xlabel('M')
        ax.set_ylim([0,1])
        ax.set_ylabel('X')

        ax2.set_ylim([self.runVars['v'].min(),self.runVars['v'].max()])
        ax2.set_ylabel('v')

        # add some quantities of interest in title
        dm = self.sim_dm()
        step = int(num * self.runConsts['saveevery'])

        ax.set_title(r'dm={:.3f}, C={:.1f}, $\gamma$={:.1e}, step={:d}'
                     .format(dm, self.courant, self.runVars['gamma'].mean(), step))

        ax.legend(lns, labs, frameon=False)

        fig.savefig(os.path.join(self.rundir, self.mytestcase + '-{:06d}'.format(num)),
                    dpi=self.dpi)
        plt.close()

