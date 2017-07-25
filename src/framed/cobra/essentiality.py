"""
This module implements methods to compute gene and reaction essentiality.

@author: Daniel Machado
   
"""

from ..solvers import solver_instance
from ..solvers.solver import Status
from simulation import FBA
from deletion import gene_deletion, reaction_deletion


def essential_genes(model, min_growth=0.01, constraints=None, solver=None):
    """ Compute the set of essential genes.
    
    Arguments:
        model (CBModel): model
        min_growth (float): minimum fraction of growth rate to consider a deletion non-letal (default: 0.01)
        constraints (dict): environmental or additional constraints (optional)

    Returns:
        list: essential genes
    """
    return essentiality(model, 'genes', min_growth, constraints, solver)


def essential_reactions(model, min_growth=0.01, constraints=None, solver=None):
    """ Compute the set of essential reactions.
    
    Arguments:
        model (CBModel): model
        min_growth (float): minimum fraction of growth rate to consider a deletion non-letal (default: 0.01)
        constraints (dict): environmental or additional constraints (optional)

    Returns:
        list: essential reactions
    """

    return essentiality(model, 'reactions', min_growth, constraints, solver)


def essentiality(model, kind='reactions', min_growth=0.01, constraints=None, solver=None):
    """ Generic interface for computing gene or reaction essentiality.
    
    Arguments:
        model (CBModel): model
        kind (str): genes or reactions (default: reactions)
        min_growth (float): minimum fraction of growth rate to consider a deletion non-letal (default: 0.01)
        constraints (dict): environmental or additional constraints (optional)

    Returns:
        list: essential elements
    """

    if solver is None:
        solver = solver_instance(model)

    wt_solution = FBA(model, constraints=constraints, solver=solver)
    wt_growth = wt_solution.fobj

    if kind == 'genes':
        elements = model.genes
    else:
        kind = 'reactions'
        elements = model.reactions

    essential = []

    for elem in elements:
        if kind == 'genes':
            solution = gene_deletion(model, [elem], constraints=constraints, solver=solver)
        else:
            solution = reaction_deletion(model, [elem], constraints=constraints, solver=solver)

        if solution and (solution.status == Status.OPTIMAL
                         and solution.fobj < min_growth * wt_growth
                         or solution.status == Status.INFEASIBLE):
            essential.append(elem)

    return essential