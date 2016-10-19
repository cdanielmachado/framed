""" This module implements flux variability analysis methods.

@author: Daniel Machado, Kai Zhuang

"""

from collections import OrderedDict
from ..solvers import solver_instance
from ..solvers.solver import Status
from simulation import FBA, looplessFBA
from numpy import linspace


def FVA(model, obj_percentage=0, reactions=None, constraints=None, loopless=False, internal=None, solver=None):
    """ Run Flux Variability Analysis (FVA).
    
    Arguments:
        model : CBModel -- a constraint-based model
        obj_percentage : float -- minimum percentage of growth rate (default 0.0, max: 1.0)
        reactions : list (of str) -- list of reactions to analyze (default: all)
        
    Returns:
        dict (of str to (float, float)) -- flux variation ranges
    """

    _constraints = {}
    if constraints:
        _constraints.update(constraints)

    if not solver:
        solver = solver_instance(model)

    if obj_percentage > 0:
        target = model.detect_biomass_reaction()
        solution = FBA(model, objective={target: 1}, constraints=constraints, solver=solver)
        _constraints[target] = (obj_percentage * solution.fobj, None)

    if not reactions:
        reactions = model.reactions.keys()

    variability = OrderedDict([(r_id, [None, None]) for r_id in reactions])

    for r_id in reactions:
        if loopless:
            solution = looplessFBA(model, {r_id: 1}, True, constraints=_constraints, internal=internal,
                                   solver=solver, get_values=False)
        else:
            solution = FBA(model, {r_id: 1}, True, constraints=_constraints, solver=solver, get_values=False)

        if solution.status == Status.OPTIMAL:
            variability[r_id][0] = solution.fobj
        elif solution.status == Status.UNBOUNDED:
            variability[r_id][0] = None
        elif solution.status == Status.INF_OR_UNB: #taking a wild guess here!
            variability[r_id][0] = None
        else:
            variability[r_id][0] = 0

    for r_id in reactions:
        if loopless:
            solution = looplessFBA(model, {r_id: 1}, False, constraints=_constraints, internal=internal,
                                   solver=solver, get_values=False)
        else:
            solution = FBA(model, {r_id: 1}, False, constraints=_constraints, solver=solver, get_values=False)
             
        if solution.status == Status.OPTIMAL:
            variability[r_id][1] = solution.fobj
        elif solution.status == Status.UNBOUNDED:
            variability[r_id][1] = None
        elif solution.status == Status.INF_OR_UNB: #taking a wild guess here!
            variability[r_id][1] = None
        else:
            variability[r_id][1] = 0

    return variability


def blocked_reactions(model, abstol=1e-9):
    """ Find all blocked reactions in a model
    
    Arguments:
        model : CBModel -- a constraint-based model
        
    Returns:
        list (of str) -- blocked reactions
    """

    variability = FVA(model)

    return [r_id for r_id, (lb, ub) in variability.items()
            if lb is not None and ub is not None and abs(lb) < abstol and abs(ub) < abstol]


def flux_envelope(model, r_x, r_y, steps=10, constraints=None):
    """ Calculate the flux envelope for a pair of reactions.

    Arguments:
        model : CBModel -- the model
        r_x : str -- reaction on x-axis
        r_y : str -- reaction on y-axis
        steps : int -- number of steps to compute (default: 10)
        constraints : dict -- custom constraints to the FBA problem

    Returns:
        list (of float), list (of float), list (of float) -- x values, y min values, y max values
    """

    x_range = FVA(model, reactions=[r_x], constraints=constraints)
    xmin, xmax = x_range[r_x]
    xvals = linspace(xmin, xmax, steps).tolist()
    ymins, ymaxs = [None] * steps, [None] * steps

    if constraints is None:
        _constraints = {}
    else:
        _constraints = {}
        _constraints.update(constraints)

    for i, xval in enumerate(xvals):
        _constraints[r_x] = xval
        y_range = FVA(model, reactions=[r_y], constraints=_constraints)
        ymins[i], ymaxs[i] = y_range[r_y]

    return xvals, ymins, ymaxs


def production_envelope(model, r_target, r_biomass=None, steps=10, constraints=None):
    """ Calculate the production envelope of the target reaction

    Arguments:
        model : CBModel -- the model
        r_target: str -- the target reaction id
        steps: int -- number of steps along the envelope to be calculated (default: 10)
        r_biomass: str -- the biomass reaction id (default: None)
        constraints : dict -- custom constraints to the FBA problem

    Returns:
        list (of float), list (of float), list (of float) -- biomass values, target minimum values, target maximum values
    """
    if not r_biomass:
        r_biomass = model.detect_biomass_reaction()

    return flux_envelope(model, r_x=r_biomass, r_y=r_target, steps=steps, constraints=constraints)


def flux_envelope_3d(model, r_x, r_y, r_z, steps=10, constraints=None):
    """ Calculate the flux envelope for a triplet of reactions.

    Arguments:
        model : CBModel -- the model
        r_x : str -- reaction on x-axis
        r_y : str -- reaction on y-axis
        r_z : str -- reaction on z-axis
        steps : int -- number of steps to compute along both x and y axis (default: 10)
        constraints : dict -- custom constraints to the FBA problem

    Returns:
        list (of float), list (of float), list (of float), list (of float) -- z min values, z max values, x coordinates, y coordinates
    """

    xvals, ymins, ymaxs = flux_envelope(model, r_x, r_y, steps, constraints)

    yvals=[None]*steps
    zmins, zmaxs = [None] * steps, [None] * steps
    x_coors, y_coors = [None] * steps, [None] * steps

    for i, xval in enumerate(xvals):

        zmins[i], zmaxs[i] = [None] * steps, [None] * steps
        x_coors[i], y_coors[i] = [None] * steps, [None] * steps

        yvals[i] = linspace(ymins[i], ymaxs[i], steps).tolist()

        for j, yval in enumerate(yvals[i]):
            x_coors[i][j] = xval
            y_coors[i][j] = yval
            x_constraint = {r_x: xval}
            y_constraint = {r_y: yval}
            constraints.update(x_constraint)
            constraints.update(y_constraint)
            z_range = FVA(model, reactions=[r_z], constraints=constraints)

            zmins[i][j], zmaxs[i][j] = z_range[r_z]
    return zmins, zmaxs, x_coors, y_coors
