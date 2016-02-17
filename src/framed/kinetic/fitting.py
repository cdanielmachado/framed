from framed.kinetic.simulation import simulate
from numpy import array, sum, isfinite
from scipy.optimize import minimize

__author__ = 'daniel'


def fit_from_metabolomics(model, t_steps, data, parameters=None, bounds=None, method=None, update_model=False):

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
        _, X = simulate(model, t_steps=t_steps, parameters=new_params)
        error = sum((X[:,mets] - X_exp)**2)
        if not isfinite(error):
            print 'warning: error = ', error
        return error

    res = minimize(fit_distance, p0, method=method, bounds=bounds)
    fitted_params = dict(zip(parameters, res.x))

    if update_model:
        model.set_parameters(fitted_params)

    return fitted_params
