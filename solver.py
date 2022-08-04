import numpy as np
import time
import math

"""
IF U1 IS PROVIDED: For solving a wave given two input waves with a time step between them
IF U1 IS NOT PROVIDED: For solving a static wave given an input wave (initial wave velocity is zero)
Endpoints are held constant at zero displacement
Momentum and energy are conserved
"""


def solve(t_bound, dx, dt, c):
    # Initialize array and special numbers for calculations
    length = c.size
    u0 = np.zeros(length, np.double)
    u1 = np.zeros(length, np.double)
    u2 = np.zeros(length, np.double)
    desired_frames = 500
    a = []
    nx = length - 1
    nt = int(t_bound // dt)
    C2 = (dt * dt) / (dx * dx)

    # Calculate u1
    for i in range(1, nx - 1):
        u1[i] = u0[i] + 0.5 * C2 * (u0[i + 1] * (c[i + 1] ** 2)
                                    - 2 * u0[i] * (c[i] ** 2)
                                    + u0[i - 1] * (c[i - 1] ** 2))

    # Calculate u2 and generate a
    t0 = time.time()
    print("Beginning Calculation")
    for n in range(1, nt + 1):
        # Update mesh points for all points in u2
        # Assumes ideally reflecting bounds
        for i in range(0, nx + 1):
            ip1 = i + 1 if i < nx else i - 1
            im1 = i - 1 if i > 0 else i + 1
            u2[i] = 2 * u1[i] - u0[i] + C2 * (u1[ip1] * (c[ip1] ** 2)
                                              - 2 * u1[i] * (c[i] ** 2)
                                              + u1[im1] * (c[im1] ** 2))

        # Send out wave
        if n < 1000:
            u2[500] = math.exp(-(5 * (n / 1000 - 0.5)) ** 2)

        # Calculate time to completion
        if n == round(nt / 100):
            print("Estimated time from completion:",
                  round((time.time() - t0) * 100, 2),
                  "seconds")

        # Update a
        if n % round((nt + 1) / desired_frames) == 0:
            a.append(np.copy(u2))

        # Change Variables
        u0[:], u1[:] = u1, u2

        # Display progress
        if n % (round(nt / 10)) == 0:
            print(round(n / (nt + 1) * 100), "% Done Calculating")

    print("Done Calculating")
    print("Real time taken:", round(time.time() - t0, 2), "seconds")

    return a
