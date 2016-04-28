__author__ = 'daniel'

from ..kinetic.simulation import find_steady_state

from numpy import array
from numpy.random import rand, randn

def sample(model, size, parameters=None, scale='log', dist='normal', dist_args=(0, 1)):

    model_params = model.get_parameters(exclude_compartments=True)
    if parameters:
        p0 = array([model_params[key] for key in parameters])
    else:
        parameters = model_params.keys()
        p0 = array(model_params.values())

    p_sample = []
    v_sample = []
    for i in range(size):
        p = parameter_perturbation(p0, scale, dist, dist_args)
        new_params = dict(zip(parameters, p))
        v = find_steady_state(model, parameters=new_params)
        if v:
            p_sample.append(p)
            v_sample.append(v.values())

    fail_rate = 100 * (size - len(v_sample))/float(size)

    if fail_rate > 10:
        print 'Warning: {}% of simulations failed.'.format(int(fail_rate))

    return p_sample, v_sample


def parameter_perturbation(p0, scale, dist, dist_args):

    if dist == 'normal':
        d = dist_args[0] + dist_args[1]*randn(len(p0))
    else:
        d = dist_args[0] + (dist_args[1] - dist_args[0]) * rand(len(p0))

    if scale == 'log10':
        d = 10 ** d
    elif scale == 'log2':
        d = 2 ** d

    return p0*d
