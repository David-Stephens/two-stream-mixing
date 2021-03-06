# two-stream-mixing

This repository contains fortran code within the `two-stream-mixing.f` which has
a subroutine, `advect`, that will *advect* nuclear species over a single time
step. The main application of this subroutine is for the mixing of species
within a convection zone when modeling the nucleosynthesis that occurs during
neutron capture processes such as the *i-process*. The original application of
this advective mixing model was on a Rapidly Accreting White Dwarf (RAWD) in this
[paper](https://ui.adsabs.harvard.edu/abs/2020arXiv200110969S/abstract) where
the algorithm is mentioned in detail.

## Getting Started ##

  * In order to compile the fortran code, it is assumed that you have installed
    *gfortran*. If you are using a Mac you can easily install this with `brew
    install gcc` which will include *gcc* and *g++* as well. To install
    *Homebrew* take a look at their [webpage](https://brew.sh/).
  
  * To run the test cases that are included in this repository, through
    `two-stream-tests.py`, you will need to have python 3.5 or later installed.
    External modules that need to be installed include `numpy` and `matplotlib`
  
  * If you want to be able to make movies of the output python plots
    automatically ensure that you have `ffmpeg` installed. On Mac it can be
    easily installed with `brew install ffmpeg`
  
  * To get help on how to run `two-stream-tests.py` just type in the terminal
    `python two-stream-tests.py -h`
  
## Test Cases ##

The entire suite of test cases can be run by simply running the command:

`python two-stream-tests.py -a`

Png's of the output text files of all tests are made by default. This option can be turned off with the `-np` flag

Movies of the output matplotlib plots, if there is no `-np` flag, are not made by default. You need to
include the `-m` or `--movie` flag like so

`python two-stream-tests.py -a -m`

All of the standard *Test Cases* have their movies, the code that was compiled
and their configuration files under the `./reference-output` folder. The fortran
code found there can be compiled and run directly rather than using
`two-stream-tests.py` however the fortran code will only compute the advection
and output text files, all other operations are done within
`two-stream-tests.py`

### 1. GaussianTest ###

The purpose of these tests is to show the advection algorithms convergence as a
function of the grid size. This test case has a constant flux,
![flux](https://latex.codecogs.com/svg.latex?2%20%5Cpi%20r%5E2%20%5Crho%20v),
throughout the domain such that
![beta](https://latex.codecogs.com/svg.latex?%5Cbeta) is always 0. There is also
no additional horizontal mixing,
![gamma](https://latex.codecogs.com/svg.latex?%5Cgamma), in this test case. This
simulation can be run simply with the command:

`python two-stream-tests.py -t GaussianTest`

The convergence results of the `GaussianTest` are shown in the figure below

![](https://github.com/David-Stephens/two-stream-mixing/blob/master/reference-output/gaussian/L2-norm.png)

### 2. BetaTest, ![beta](https://latex.codecogs.com/svg.latex?%5Cbeta)  ###

The purpose of this test is to showcase the influence of the
![beta](https://latex.codecogs.com/svg.latex?%5Cbeta) coefficient, the enforced
horizontal mixing, on the advection of a gaussian. Here there is a quadratic
velocity field that is shown in the output figures. This simulation can be run
simply with the command:

`python two-stream-tests.py -t BetaTest`

### 3. GammaTest, ![gamma](https://latex.codecogs.com/svg.latex?%5Cgamma) ###

The purpose of this test is to show the influence of the
![gamma](https://latex.codecogs.com/svg.latex?%5Cgamma) coefficient on the
advection of a gaussian as well as the inclusion of
![beta](https://latex.codecogs.com/svg.latex?%5Cbeta). The
![gamma](https://latex.codecogs.com/svg.latex?%5Cgamma) coefficient is constant
throughout the entire simulation and its value can be seen in the
`two-stream-tests.config` parameters file. That file contains all of the
specifications required for the test cases and they can be changed to create new
test cases. Checkout the python scripts under `./utils` for more information. This
simulation can be run simply with the command:

`python two-stream-tests.py -t gammaTest`

## Applying to Post-Processing Code ##

In principle, this subroutine can be used as a replacement to any diffusive
mixing routines assuming that the burn and mix solvers are computed separately.
The main limiting feature of this code is the fact that the number of grid cells
per stream must remain constant throughout the entire simulation. This was
enforced due to the significant increase in code complexity when allowing for
the streams changing size on a time step basis. An overview of how the
subroutine can be functionally used, as well as extensive doc strings for each
subroutine, can be found in any of the `./reference-output` fortran files.

Another issue which has not been addressed by this code release is that the code
is not as efficient as it could be and it wasn't written in a 1-to-1
correspondence with the algorithm described in the
[paper](https://ui.adsabs.harvard.edu/abs/2020arXiv200110969S/abstract). This
may be changed in the near future as a new version however it will require
considerable re-writing of the code. Stay tuned!
