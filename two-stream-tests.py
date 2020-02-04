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

    def __init__(self, root, testcase, movies):
        """
        Constructor

        Parameters
        ----------

        root : str
            The absolute path of where this python script was ran from

        testcase : str
            The testcase chosen by the user

        movies : bool
            A boolean as to whether or not we make movies by default
        """

        # create the parameters object to have access to the parameters file
        self.root = root
        self.params = defaults.Params(self.root)

        # are we doing movies?
        self.movies = movies

        # check if the parameters file exists in "home"
        if not os.path.isfile(self.params.paramfile):
            self.params.writeTestcases()

        if testcase not in defaults.Params.testcases:
            print('Error: The {:s} testcase is not in the parameters file {:s}'
                .format(self.params.paramfile))

        # ok, our testcase is
        self.testcase = testcase

    def run_simulation(self):
        """
        Based on the testcase, create the appropriate object and run the simulation
        """

        # create appropriate object
        simulation = self.all_sims[self.testcase](self.testcase, self.params.paramfile, self.root,
                                            self.movies)

        simulation.run_simulation()


def clean(root):
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
    shparams.add_argument("--testcase", '-t', type=str, help="The testcase to run",
                          choices=defaults.Params.testcases)
    shparams.add_argument('--movie', '-m', help="Do we create movies from the plotted outputs?",
                          action="store_true")
    shparams.add_argument('--clean', '-c', help="Clean the directory of all runfiles",
                          action="store_true")

    # parse the args
    args = shparams.parse_args()

    # where did we run this?
    root = os.path.dirname(os.path.realpath(__file__))

    if args.clean:
        clean(root)
        print('All of the run directories in this repository have been removed')

    else:
        # create the twostream object
        testing = TwoStreams(root, args.testcase, args.movie)
        testing.run_simulation()

        print('We have now finished the testcase, {:s}'.format(args.testcase))

        if args.movie:
            print('You can find the png and movie output of the simulation within ./{:s}'
                  .format(args.testcase))
        else:
            print('You can find the png output of the simulation within ./{:s}'.format(args.testcase))
