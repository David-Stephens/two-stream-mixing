from utils import defaults, testcases
import shutil
import argparse
import os

class TwoStreams(object):
    """
    Create and choose the corresponding testcase. Run the simulation
    """

    # simulation objects that can be made with their testcase key
    all_sims = dict(zip(defaults.Params.testcases, [testcases.Gaussian, testcases.BetaTest,
                                                    testcases.GammaTest]))

    # the actual testcases have the limit of 8 characters (due to fortran....) and so for a cleaner
    # presentation we use the testcase_search strings as a front to those actual testcases
    testcase_search = dict(zip(['GaussianTest', 'BetaTest', 'GammaTest'], defaults.Params.testcases))

    def __init__(self, root, testcase, plots, movies):
        """
        Constructor

        Parameters
        ----------

        root : str
            The absolute path of where this python script was ran from

        testcase : str
            The testcase chosen by the user

        plots : bool
            A boolean as to whether or not we make plots by default

        movies : bool
            A boolean as to whether or not we make movies by default
        """

        # create the parameters object to have access to the parameters file
        self.root = root
        self.params = defaults.Params(self.root)

        # are we doing movies and/or plots?
        self.plots = plots
        self.movies = movies

        # check if the parameters file exists in "home"
        if not os.path.isfile(self.params.paramfile):
            self.params.writeTestcases()

        # ok, our fortran string testcase is instead
        self.testcase = self.testcase_search[testcase]

    def finishedTestcase(self):
        """
        Print text to notify user that the testcase was finished correctly
        """

        print('We have now finished the testcase, {:s}'.format(self.testcase))

        if self.movies:
            print('You can find the png and movie output of the simulation within ./{0:s}/png and ./{0:s}'
                  .format(self.testcase))
        else:
            print('You can find the png output of the simulation within ./{:s}/png'.format(self.testcase))

    def run_simulation(self):
        """
        Based on the testcase, create the appropriate object and run the simulation
        """

        # create appropriate object
        simulation = self.all_sims[self.testcase](self.testcase, self.params.paramfile, self.root,
                                                  self.plots, self.movies)

        simulation.run_simulation()
        self.finishedTestcase()


def cleanRuns(root):
    """
    Clean the repository directory

    Parameters
    ----------

    root : str
        The absolute path of where this python script was ran from
    """

    for testcase in defaults.Params.testcases:

        # if the directory exists, we will delete it
        if os.path.isdir(os.path.join(root, testcase)):
            shutil.rmtree(os.path.join(root, testcase))

# if run..
if __name__ == "__main__":

    # The arguments to the script
    shparams = argparse.ArgumentParser()
    shparams.add_argument("-a", '--allTests', help="Run all testcases",
                          action="store_true")
    shparams.add_argument("-t", "--testcase", type=str, help="Testcase(s) to run",
                          choices=TwoStreams.testcase_search.keys(), nargs='+')
    shparams.add_argument("-np", '--noPlots', help="Do not create any plots",
                          action="store_false")
    shparams.add_argument('-m', '--movies', help="Do we create movies from the plotted outputs?",
                          action="store_true")
    shparams.add_argument('-c', '--clean', help="Clean the directory of all runfiles",
                          action="store_true")

    # parse the args
    args = shparams.parse_args()

    # where did we run this?
    root = os.path.dirname(os.path.realpath(__file__))

    allTests = args.allTests
    someTests = args.testcase
    plots = args.noPlots
    movies = args.movies
    clean = args.clean

    # we can't make movies if we aren't making plots
    if not plots:
        movies = False

    # are we cleaning the directory?
    if clean:
        cleanRuns(root)
        print('All of the run directories in this repository have been removed')

    else:

        # If we run all...
        if allTests:
            testsToRun = list(TwoStreams.testcase_search.keys())
        else:
            testsToRun = someTests

        # loop through the testcases we need to run
        for testcase in testsToRun:

            # create the twostream object
            testing = TwoStreams(root, testcase, plots, movies)
            testing.run_simulation()

