"""
This module implements kinetic model calibration methods.

Author: Daniel Machado
"""


from framed.kinetic.simulation import time_course
from numpy import array, sum, isfinite
from scipy.optimize import minimize


def fit_from_metabolomics(model, t_steps, data, parameters=None, bounds=None, method=None, update_model=False):
    """ Fit model parameters using time-course metabolomics data

    Args:
        model (ODEModel): kinetic model
        t_steps (list): measured time steps
        data (dict): metabolomics data in dict format (metabolite id to meausured values)
        parameters (list): specify list of parameters to be calibrated (optional, default: all)
        bounds (list): list of bounds for each parameter (optional, default: (0, None))
        method (str): optimization method (optional, see `scipy.optimize.minimize` for details)
        update_model (bool): automatically update model with new parameters (default: False)

    Returns:
        dict: fitted parameters

    """

    model_params = model.get_parameters(exclude_compartments=True)
    if parameters:
        p0 = [model_params[key] for key in parameters]
    else:
        parameters = model_params.keys()
        p0 = model_params.values()

    if not bounds:
        bounds = [(0, None)]*len(p0)

    mets = [model.metabolites.keys().index(m_id) for m_id in data.keys()]
    X_exp = array(data.values()).T

    def fit_distance(p):
        new_params = dict(zip(parameters, p))
        _, X = time_course(model, t_steps=t_steps, parameters=new_params)
        error = sum((X[:,mets] - X_exp)**2)
        if not isfinite(error):
            print 'warning: error = ', error
        return error

    res = minimize(fit_distance, p0, method=method, bounds=bounds)
    fitted_params = dict(zip(parameters, res.x))

    if update_model:
        model.set_parameters(fitted_params)

    return fitted_params
