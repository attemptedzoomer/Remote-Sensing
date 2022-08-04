import numpy as np
from solver import solve


def generate_scene(lxb, rxb, n):
    # N-point uniform mesh
    dx = (rxb - lxb) / (n - 1)
    x0 = np.linspace(lxb, rxb, n)  # returns an array of n values evenly spaced between 0 and 1

    # Sound speed setup
    c = np.ones(n, np.double)
    damping_coefficient = .99
    for i in range(0, c.size):
        if i < int(round(n * 0.41)):
            c[i] *= 3
        if i < int(round(n * 0.59)):
            c[i] *= 3
    # c_temp = np.copy(c)
    # for i in range(20, n - 20):
    #     c_temp[i] = np.average(c[i - 20:i + 20])
    # c = c_temp
    for i in range(n-1, 0, -1):
        if x0[i] < 0:
            c[i] = c[i+1] * damping_coefficient
    for i in range(0, n-1, 1):
        if x0[i] > 1:
            c[i] = c[i-1] * damping_coefficient
        round(c[i], 4)

    # Calculate based on C
    C = 0.5
    dt = C / np.max(c) * dx

    # Quality check
    print("CFL constant is {0} (should be < 1 for stability)".format(C))

    # Space for time steps
    t_bound = 3
    a = solve(t_bound, dx, dt, c)

    return a
