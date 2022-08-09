import numpy as np


def generate_pulse(array_length=1000,
                   starting_location='left',
                   pulse_type='gauss',
                   sigma=0.05,
                   left_domain_bound=0,
                   right_domain_bound=1,
                   ):
    """
    Generates a pulse of a specific type and returns the pulse for further propagation

    :param array_length: length of the target array within the VISIBLE DOMAIN
    :param starting_location: string value of where the center of the wave should be
    :param pulse_type: shape of the wave
    :param sigma: width of the generated wave
    :param left_domain_bound: the x value of the left bound of the domain
    :param right_domain_bound: the x value of the right bound of the domain
    :return: returns an u0 with a pre-initialized pulse based on args
    """

    L = right_domain_bound - left_domain_bound

    if starting_location == 'center':
        xc = L / 2
    elif starting_location == 'right':
        xc = L
    else:
        xc = 0

    if pulse_type == 'gauss':
        def I(x):
            return np.exp(-0.5 * ((x - xc) / sigma) ** 2)
    elif pulse_type == 'plug':
        def I(x):
            return 0 if abs(x - xc) > sigma else 1
    elif pulse_type == 'cosinehat':
        def I(x):
            a = 2 * sigma
            return 0.5 * (1 + np.cos(np.pi * (x - xc))) if xc - a <= x <= xc + a else 0
    elif pulse_type == 'half-cosinehat':
        def I(x):
            a = 4 * sigma
            return np.cos(np.pi * (x - xc) / a) if xc - 0.5 * a <= x <= xc + 0.5 * a else 0
    else:
        raise ValueError(f'Invalid pulse type, entered {pulse_type}')

    u0 = np.zeros(array_length, np.float32)
    index_set = np.linspace(left_domain_bound, right_domain_bound, array_length)

    for i in range(array_length - 1):
        u0[i] = I(index_set[i])

    return u0
