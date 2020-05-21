import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as scifft
import scipy.integrate as sciint


def initiate(time=0.3, xsize=10, dx=0.1, dt=1.0e-2,
             temp_src=373., temp_amb=296.,
             silent=False):
    """Initialize the parameters

    The default

    Args:
        u0 {numpy.ndarray} -- 1D numpy array containing the initial condition.
        D {float} -- the heat coefficient. D=23 mm²/2 for iron.
        time {numpy.ndarray} -- real simulation time. time=0.3 s.
        xsize {float} -- the length of the simulated system. size=10 mm.
        dx {float} -- the length spacing. dx=0.1 mm.
        dt {float} -- time step.
        temp_src {float} -- hottest temperature for the initial condition.
        temp_amp {float} -- coldest temperature for the initial condition.
        silent {float} -- whether to print messages.
    returns:
        u0 {numpy.ndarray} -- the initial condition
            as a gaussian distribution od temperature.
            Shape (time/dt, xscale/dx).
        t {numpy.ndarray} -- instants of time. Shape(time/dt,)
        xscale {numpy.ndarray} -- space grid. Shape (xscale/dx, )"""

    xscale = np.arange(-0.5*xsize, 0.5*xsize, dx)
    ntime = int(np.round(time/dt))
    t = np.linspace(0, time, int(ntime))
    Nx = len(xscale)

    u0 = (temp_src - temp_amb)*np.exp(-xscale**2/5) + temp_amb

    if not silent:
        print(f"Number of time steps: {ntime}")
        print(f"Number of points: {Nx}")
        print(f"Time step: {dt:.3e}")

    return u0, t, xscale


def dudt(u, D, k):
    """The righthand side of the heat equation

    Arguments:
        u {numpy.ndarray} -- The temperature distribution
        D {float} -- Diffusive coeffitient
        k {numpy.ndarray} -- fft frequencies

    Returns:
        numpy.ndarray -- rhs of the heat equation derived
    """
    d2u_hat = -np.power(k, 2)*scifft.fft(u)
    d2u = scifft.ifft(d2u_hat).real
    d2u[0] = d2u[-1] = 0
    return D*d2u


def kappa(u0, dx):
    return 2*np.pi*scifft.fftfreq(len(u0), d=dx)


def solve(rhs, u0, t, D, kappa):
    """Integrate the heat equation

    Arguments:
        u0 {numpy.ndarray} -- initial condition
        rhs {function} -- righthand side of the heat equation
        t {numpy.ndarray} -- instants of time

    Returns:
        numpy.ndarray -- set of solutions for every instant of time
    """
    return sciint.odeint(lambda u, t: rhs(u, D, kappa), u0, t)


if __name__ == '__main__':

    D = 23  # for iron
    dx = 0.1
    u0, t, xscale = initiate(dx=dx)
    k = kappa(u0, dx)

    u = solve(dudt, u0, t, D, k)

    plt.figure(figsize=(8, 8))
    plt.contourf(xscale, t, u, 100, cmap="magma")
    plt.xlabel(r"$x(mm)$", fontsize=16)
    plt.ylabel(r"time (s)", fontsize=16)
    cbar = plt.colorbar()
    cbar.set_label("Temperature", fontsize=16, rotation=270,
                    labelpad=25, y=0.45)
    plt.show()
