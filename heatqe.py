import numpy as np
#  import matplotlib.pyplot as plt
import scipy.fftpack as scifft


def dutdt(u, t, D, k):
    d2u_hat = -np.power(k, 2)*scifft.fft(u)
    d2u = scifft.fft(d2u_hat).real
    d2u[0] = d2u[-1] = 0
    return D*d2u


def initiate(D=23, time=0.3,
             xsize=10, dx=0.1, dt=1.0e-2,
             temp_src=373., temp_amb=296.,
             silent=False):

    """Initialize the parameters

    The default

    Args:
        u0: 1D numpy array containing the initial condition.
        D: the heat coefficient. D=23 mm²/2 for iron.
        time: real simulation time. time=0.3 s.
        xsize: the length of the simulated system. size=10 mm.
        dx: the length spacing. dx=0.1 mm.
        dt: time step.
        temp_src: hottest temperature for the initial condition.
        temp_amp: coldest temperature for the initial condition.
        silent: whether to print messages.
    returns:
        u0: numpy array containing solution. Shape (time/dt, xscale/dx).
        t: numpy array containing all instants. Shape(time/dt,)
        xscale: numpy array containg space grid. Shape (xscale/dx, )"""

    xscale = np.arange(-0.5*xsize, 0.5*xsize, dx)
    ntime = int(np.round(time/dt))
    t = np.linspace(0, time, int(ntime))
    Nx = len(xscale)

    u0 = (temp_src - temp_amb)*np.exp(-2*xscale**2/xsize) + temp_amb

    if not silent:
        print(f"Number of time steps: {ntime}")
        print(f"Number of points: {Nx}")
        print(f"Time step: {dt:.3e}")

    return u0, t, xscale
