"""
This module implements time-course and steady-state simulation for kinetic models.

Author: Daniel Machado
"""


from numpy import linspace, array, dot, isnan
from numpy.linalg import norm
from scipy.integrate import odeint
from collections import OrderedDict
from warnings import warn


def time_course(model, time=0, steps=100, t_steps=None, parameters=None, compute_rates=False, integrator_args=None):
    """ Perform a time-course simulation using a kinetic model.

    Args:
        model (ODEModel): kinetic model
        time (float): final simulation time (optional if t_steps is used instead)
        steps (int): number of simulations steps (default: 100)
        t_steps (list): list of exact time steps to evaluate (optional)
        parameters (dict): override model parameters (optional)
        integrator_args (dict): additional parameters to pass along to scipy odeint integrator

    Returns:
        tuple: time steps, metabolite concentrations[, reaction rates (optional)]

    """
    r = OrderedDict()
    f = model.get_ode(r, parameters)
    f2 = lambda x, t: f(t, x)
    x0 = model.concentrations.values()

    if t_steps is None:
        t_steps = linspace(0, time, steps)

    if not integrator_args:
        integrator_args = {'mxstep': 10000}

    X = odeint(f2, x0, t_steps, **integrator_args)

    if compute_rates:
        return t_steps, X, r
    else:
        return t_steps, X


# def simulate2(model, time, n_steps=100, method='vode', integrator_args=None):
#     #simulation based on scipy.integrate.ode (for some reason much slower than odeint)
#
#     f = model.get_ode()
#     x0 = model.concentrations.values()
#
#     if not integrator_args:
#         if method == 'vode':
#             integrator_args = {'method': 'bdf', 'order': 15}
#         else:
#             integrator_args = {}
#
#     integrator = ode(f).set_integrator(method, **integrator_args)
#     integrator.set_initial_value(x0)
#
#     dt = float(time)/n_steps
#     t = [0]
#     X = [x0]
#     while integrator.successful() and integrator.t < time:
#         integrator.integrate(integrator.t + dt)
#         t.append(integrator.t)
#         X.append(integrator.y)
#
#     return t, array(X)


def find_steady_state(model, parameters=None, endtime=1e9, abstol=1e-6):
    """ Determine steady-state

    Args:
        model (ODEModel): kinetic model
        parameters (dict): override model parameters
        endtime (float): final integration time (default: 1e9)
        abstol (float): maximum tolerance for norm(S*v)

    Returns:
        tuple: steady-state concentrations, steady-state fluxes

    """

    _, X, v_ss = time_course(model, t_steps=[0, endtime], parameters=parameters, compute_rates=True)

    x_ss = OrderedDict(zip(model.metabolites.keys(), X[-1, :]))

    S = array(model.stoichiometric_matrix())
    error = norm(dot(S, v_ss.values()))

    if error > abstol or isnan(error):
        warn('Simulation did not reach a steady state.')
        v_ss = None

    return x_ss, v_ss

