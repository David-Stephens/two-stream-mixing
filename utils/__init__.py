import numpy as np

# Generic functions used throughout
def gaussian(mc, m_shell, sigma, offset=5., limit=True, extent=5.):
    """
    Create a cell averaged gaussian distribution of the mass fraction of species X in the mass coordinates given

    Parameters
    ----------

    mc : np.ndarray
        The mass coordinates of the center of the cell

    m_shell : np.ndarray
        The mass coordinates at the interfaces of the cell (shell)

    sigma : float
        The sigma of the gaussian distribution.

    offset : float
        The offset where the mass coordinate of the mean is located

    limit : bool
        Are we going to artificially limit the gaussian for clarity in images?

    extent : float
        The extent in a multiple of sigma to which the gaussian is resolved to when limited

    Returns
    -------
    gauss : np.ndarray
        The mass fraction of the initialized X species with a gaussian model
    """

    # Take an average based on the borders of the gaussian
    gauss = np.zeros(len(mc))

    ml = m_shell[0:-1]
    mr = m_shell[1:]
    
    gaussl = np.exp(-0.5*np.power(((ml - offset)/sigma),2.0))
    gaussr = np.exp(-0.5*np.power(((mr - offset)/sigma),2.0))
    gaussc = np.exp(-0.5*np.power(((mc - offset)/sigma),2.0))

    gauss = 1/6.0 * (gaussl + 4*gaussc + gaussr)

    # Make it zero above and below the 5 sigma limits
    if limit:
        sigma5_uindex = np.argmin(abs(offset + extent*sigma - mc)) + 1
        gauss[sigma5_uindex:] = 0.

        sigma5_dindex = np.argmin(abs(offset - extent*sigma - mc))
        if sigma5_dindex != 0:
            gauss[0:sigma5_dindex] = 0.

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

