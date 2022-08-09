import numpy as np
import time
from numba import njit
from pulse import generate_pulse

"""
IF U1 IS PROVIDED: For solving a wave given two input waves with a time step between them
IF U1 IS NOT PROVIDED: For solving a static wave given an input wave (initial wave velocity is zero)
Endpoints are held constant at zero displacement
Momentum and energy are conserved
"""


def solve(t_bound, dx, dt, c):
    """
    Creates an array of the numerical solution to the acoustic wave equation in one dimension

    :param t_bound: the end of the time domain (the beginning is assumed to be t = 0)
    :param dx: the density of the x mesh
    :param dt: timestep
    :param c: an array representing the speed of sound of the media
    :return: a, an array that holds the behavior of the wave
    """

    # Initialize array and special numbers for calculations
    length = c.size
    desired_frames = 500
    a = []
    nx = length - 1
    nt = int(t_bound // dt)
    C2 = (dt * dt) / (dx * dx)

    # Calculate u0
    u0 = generate_pulse(length)

    # Calculate u1
    u1 = initialize_u1(u0, nx, C2, c)

    # Calculate u2 and generate a
    t0 = time.time()
    print("Beginning Calculation")
    for n in range(1, nt + 1):
        # Update mesh points for all points in u2
        # Assumes ideally reflecting bounds
        u2 = update_inner_mesh_points(u0, u1, nx, C2, c)

        # Update a
        if n % round((nt + 1) / desired_frames) == 0:
            a.append(np.copy(u2))

        # Change Variables
        u0[:], u1[:] = u1, u2

    print("Done Calculating")
    print("Real time taken:", round(time.time() - t0, 2), "seconds")

    return a


@njit
def initialize_u1(u0, nx, C2, c):
    u1 = np.zeros(u0.size, np.float32)
    for i in range(1, nx - 1):
        u1[i] = u0[i] + 0.5 * C2 * (u0[i + 1] * (c[i + 1] ** 2)
                                    - 2 * u0[i] * (c[i] ** 2)
                                    + u0[i - 1] * (c[i - 1] ** 2))
    return u1


@njit
def update_inner_mesh_points(u0, u1, nx, C2, c):
    u2 = np.zeros(u0.size, np.float32)
    for i in range(0, nx + 1):
        ip1 = i + 1 if i < nx else i - 1
        im1 = i - 1 if i > 0 else i + 1
        u2[i] = 2 * u1[i] - u0[i] + C2 * (u1[ip1] * (c[ip1] ** 2)
                                          - 2 * u1[i] * (c[i] ** 2)
                                          + u1[im1] * (c[im1] ** 2))
    return u2
