"""
This module implements flux space sampling for kinetic models.

Author: Daniel Machado
"""
from ..kinetic.simulation import find_steady_state
import numpy as np

import warnings


def sample_kinetic_model(model, size, parameters=None, distribution='normal', dist_args=(0, 1), log_scale=True):
    """ Random flux sampling using a kinetic model.

    Args:
        model (ODEModel): kinetic model
        size (int): sample size
        parameters (list): list of parameters that can vary (default: all)
        distribution (str): parameter sampling distribution ('normal' (default) or 'uniform')
        dist_args (tuple): distribution specific arguments:
                            (mean, std) for normal distribution (default: (0,1))
                            (min, max) for uniform distribution
        log_scale (bool): perform variation in log scale (default: True)

    Returns:
        tuple: parameter samples, flux samples

    """

    model_params = model.get_parameters(exclude_compartments=True)
    if parameters:
        p0 = np.array([model_params[key] for key in parameters])
    else:
        parameters = model_params.keys()
        p0 = np.array(model_params.values())

    p_sample = []
    v_sample = []
    for i in range(int(size)):
        p = parameter_perturbation(p0, distribution, dist_args, log_scale)
        new_params = dict(zip(parameters, p))
        _, v = find_steady_state(model, parameters=new_params)
        if v:
            p_sample.append(p)
            v_sample.append(v.values())

    fail_rate = 100 * (size - len(v_sample))/float(size)

    if fail_rate > 10:
        warnings.warn('Warning: {}% of simulations failed.'.format(int(fail_rate)), RuntimeWarning)

    return p_sample, v_sample


def parameter_perturbation(p0, distribution, dist_args, log_scale=True):

    if distribution == 'normal':
        d = dist_args[0] + dist_args[1] * np.random.randn(len(p0))
    elif distribution == 'uniform':
        d = dist_args[0] + (dist_args[1] - dist_args[0]) * np.random.rand(len(p0))
    else:
        warnings.warn("Invalid sampling distribution '{}'".format(distribution))
        return

    if log_scale:
        d = 10 ** d

    return p0*d
