import numpy as np

# Generic functions used throughout
def gaussian(m, sigma, offset=5., extent=5.):
    """
    Create a gaussian distribution of the mass fraction of species X in the mass coordinates given

    Parameters
    ----------

    m : np.ndarray
        The mass coordinates to determine

    sigma : float
        The sigma of the gaussian distribution.

    offset : float
        The offset, in a multiple of sigma,  where the mass coordinate of the mean is located

    extent : float
        The extent in a multiple of sigma to which the gaussian is resolved to

    Returns
    -------
    gauss : np.ndarray
        The mass fraction of the initialized X species with a gaussian model
    """

    gauss = np.exp(-0.5*np.power(((m - offset*sigma)/sigma),2.0))

    # note, argmin turns the first occurence, i.e one below where I want
    sigma5_index = np.argmin(abs(offset*sigma + extent*sigma - m)) + 1
    gauss[sigma5_index:] = 0.

    return gauss

def v_quad(mass, vmax, starti, stopi):
    """
    Create a quadratic velocity field between the mass indices starti and stopi with a maximum
    velocity of vmax.

    Parameters
    ----------

    mass : np.ndarray
        The mass coordinate of the velocity field (shell)

    vmax : float
        The maximum velocity

    starti : int
        The starting index to where the parabola will be fit in mass coordinates

    stopi : int
        The ending index  to where the parabola will be fit in mass coordinates

    Returns
    -------

    v : np.ndarray
        The velocity field that is quadratically fit between mass[starti:stopi+1]
    """

    # let m = 0 be at mass[starti] for simpler a,b,c coefficient determinations
    shift_axis = mass[starti:stopi+1] - mass[starti]
    midpointi = int(len(shift_axis)/2) - 1

    # trivially...
    c = 1.
    a = (vmax - c) / (shift_axis[midpointi]**2 -
                                    2*shift_axis[midpointi]*shift_axis[midpointi])
    b = -2*a*shift_axis[midpointi]

    return a*shift_axis**2 + b*shift_axis + c

