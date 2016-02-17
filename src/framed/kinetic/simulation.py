__author__ = 'daniel'

from numpy import linspace, array
from numpy.linalg import norm
from scipy.integrate import odeint, ode
from collections import OrderedDict
from warnings import warn


def simulate(model, time=0, steps=100, t_steps=None, parameters=None, integrator_args=None):
    f = model.get_ODEs(parameters)
    f2 = lambda x, t: f(t, x)
    x0 = model.concentrations.values()
    if t_steps is None:
        t_steps = linspace(0, time, steps)
    if not integrator_args:
        integrator_args = {'mxstep': 10000}
    X = odeint(f2, x0, t_steps, **integrator_args)
    return t_steps, X


def simulate2(model, time, n_steps=100, method='vode', integrator_args=None):
    #simulation based on scipy.integrate.ode (for some reason much slower than odeint)

    f = model.get_ODEs()
    x0 = model.concentrations.values()

    if not integrator_args:
        if method == 'vode':
            integrator_args = {'method': 'bdf', 'order': 15}
        else:
            integrator_args = {}

    integrator = ode(f).set_integrator(method, **integrator_args)
    integrator.set_initial_value(x0)

    dt = float(time)/n_steps
    t = [0]
    X = [x0]
    while integrator.successful() and integrator.t < time:
        integrator.integrate(integrator.t + dt)
        t.append(integrator.t)
        X.append(integrator.y)

    return t, array(X)


def find_steady_state(model, parameters=None, endtime=1e9, abstol=1e-6):

    _, X = simulate(model, t_steps=[0, endtime], parameters=parameters)
    x_ss = X[-1,:]
    v_ss = None
    f = model.get_ODEs(parameters)
    error = norm(f(0, x_ss))

    if error < abstol:
        p = model.get_parameters().values()
        v_ss = OrderedDict([(r_id, rate(x_ss, p)) for r_id, rate in model.rates.items()])
    else:
        warn('Simulation did not reach a steady state.')

    return x_ss, v_ss

