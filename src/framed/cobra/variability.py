""" This module implements a few flux variability analysis methods.

Authors: Daniel Machado, Kai Zhuang

"""

from collections import OrderedDict
from ..solvers import solver_instance
from ..solvers.solver import Status
from simulation import FBA
from thermodynamics import looplessFBA
from numpy import linspace
from warnings import warn


def FVA(model, obj_percentage=0, reactions=None, constraints=None, loopless=False, internal=None, solver=None):
    """ Run Flux Variability Analysis (FVA).
    
    Arguments:
        model (CBModel): a constraint-based model
        obj_percentage (float): minimum percentage of growth rate (default 0.0, max: 1.0)
        reactions (list): list of reactions to analyze (default: all)
        constraints (dict): additional constraints (optional)
        loopless (bool): run looplessFBA internally (very slow) (default: false)
        internal (list): list of internal reactions for looplessFBA (optional)
        solver (Solver): pre-instantiated solver instance (optional)
        
    Returns:
        dict: flux variation ranges
    """

    _constraints = {}
    if constraints:
        _constraints.update(constraints)

    if not solver:
        solver = solver_instance(model)

    if obj_percentage > 0:
        target = model.biomass_reaction
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
            pass
        elif solution.status == Status.INF_OR_UNB:
            pass
        elif solution.status == Status.INFEASIBLE:
            warn('Infeasible solution status')
        else:
            warn('Unknown solution status')

    for r_id in reactions:
        if loopless:
            solution = looplessFBA(model, {r_id: 1}, False, constraints=_constraints, internal=internal,
                                   solver=solver, get_values=False)
        else:
            solution = FBA(model, {r_id: 1}, False, constraints=_constraints, solver=solver, get_values=False)
             
        if solution.status == Status.OPTIMAL:
            variability[r_id][1] = solution.fobj
        elif solution.status == Status.UNBOUNDED:
            pass
        elif solution.status == Status.INF_OR_UNB:
            pass
        elif solution.status == Status.INFEASIBLE:
            warn('Infeasible solution status')
        else:
            warn('Unknown solution status')

    return variability


def blocked_reactions(model, constraints=None, reactions=None, abstol=1e-9):
    """ Find all blocked reactions in a model
    
    Arguments:
        model (CBModel): a constraint-based model
        constraints (dict): additional constraints (optional)
        reactions (list): List of reactions which will be tested (default: None, test all reactions)
        abstol (float): absolute tolerance (default: 1e-9)
        
    Returns:
        list: blocked reactions
    """

    variability = FVA(model, reactions=reactions, constraints=constraints)

    return [r_id for r_id, (lb, ub) in variability.items()
            if lb is not None and ub is not None and abs(lb) < abstol and abs(ub) < abstol]


def flux_envelope(model, r_x, r_y, steps=10, constraints=None):
    """ Calculate the flux envelope for a pair of reactions.

    Arguments:
        model (CBModel): the model
        r_x (str): reaction on x-axis
        r_y (str): reaction on y-axis
        steps (int): number of steps to compute (default: 10)
        constraints (dict): custom constraints to the FBA problem

    Returns:
        tuple: x values, y min values, y max values
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
        model (CBModel): the model
        r_target (str): the target reaction id
        steps (int): number of steps along the envelope to be calculated (default: 10)
        r_biomass (str): the biomass reaction id (optional)
        constraints (dict): custom constraints to the FBA problem

    Returns:
        tuple: biomass values, target minimum values, target maximum values
    """
    if not r_biomass:
        r_biomass = model.biomass_reaction

    return flux_envelope(model, r_x=r_biomass, r_y=r_target, steps=steps, constraints=constraints)


def flux_envelope_3d(model, r_x, r_y, r_z, steps=10, constraints=None):
    """ Calculate the flux envelope for a triplet of reactions.

    Arguments:
        model (CBModel): the model
        r_x (str): reaction on x-axis
        r_y (str): reaction on y-axis
        r_z (str): reaction on z-axis
        steps (int): number of steps to compute along both x and y axis (default: 10)
        constraints (dict): custom constraints to the FBA problem

    Returns:
        tuple: z min values, z max values, x coordinates, y coordinates
    """

    xvals, ymins, ymaxs = flux_envelope(model, r_x, r_y, steps, constraints)

    yvals = [None]*steps
    zmins, zmaxs = [None] * steps, [None] * steps
    x_coors, y_coors = [None] * steps, [None] * steps

    if constraints is None:
        _constraints = {}
    else:
        _constraints = {}
        _constraints.update(constraints)

    for i, xval in enumerate(xvals):

        zmins[i], zmaxs[i] = [None] * steps, [None] * steps
        x_coors[i], y_coors[i] = [None] * steps, [None] * steps

        yvals[i] = linspace(ymins[i], ymaxs[i], steps).tolist()

        for j, yval in enumerate(yvals[i]):
            x_coors[i][j] = xval
            y_coors[i][j] = yval
            x_constraint = {r_x: xval}
            y_constraint = {r_y: yval}
            _constraints.update(x_constraint)
            _constraints.update(y_constraint)
            z_range = FVA(model, reactions=[r_z], constraints=_constraints)

            zmins[i][j], zmaxs[i][j] = z_range[r_z]
    return zmins, zmaxs, x_coors, y_coors
