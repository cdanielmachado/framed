""" This module implements gene and reaction deletion methods.

Author: Daniel Machado
   
"""

from simulation import FBA, pFBA, MOMA, lMOMA, ROOM


def gene_deletion(model, genes, method='FBA', reference=None, constraints=None, solver=None, compute_silent_deletions=False):
    """ Simulate the deletion of a set of genes.
    
    Arguments:
        model (CBModel): model
        genes (list): genes to delete
        method (str): simulation method: FBA (default), pFMA, MOMA, lMOMA, ROOM
        reference (dict): reference flux distribution for MOMA, lMOMA or ROOM (optional)
        constraints (dict): additional constraints
        solver (Solver): solver instance instantiated with the model, for speed (optional)
        compute_silent_deletions (bool): don't compute gene deletion if no reactions are affected (optional, default: True)

    Returns:
        Solution: solution
    """

    inactive_reactions = deleted_genes_to_reactions(model, genes)

    if inactive_reactions or compute_silent_deletions:
        solution = reaction_deletion(model, inactive_reactions, method, reference, constraints, solver)
    else:
        solution = None

    return solution


def deleted_genes_to_reactions(model, genes):
    """ Convert a set of deleted genes to the respective deleted reactions.
    
    Arguments:
        model (CBModel): model
        genes (list): genes to delete

    Returns:
        list: list of deleted reactions
    """
    active_genes = set(model.genes) - set(genes)
    active_reactions = model.evaluate_gprs(active_genes)
    inactive_reactions = set(model.reactions) - set(active_reactions)

    return inactive_reactions


def reaction_deletion(model, reactions, method='FBA', reference=None, constraints=None, solver=None):
    """ Simulate the deletion of a set of reactions.
    
    Arguments:
        model (CBModel): model
        reactions (list): reactions to delete
        method (str): simulation method: FBA (default), pFMA, MOMA, lMOMA, ROOM
        reference (dict): reference flux distribution for MOMA, lMOMA or ROOM (optional)
        constraints (dict): additional constraints
        solver (Solver): solver instance instantiated with the model, for speed (optional)

    Returns:
        Solution: solution
    """

    _constraints = {}

    if constraints:
        _constraints.update(constraints)

    for r_id in reactions:
        _constraints[r_id] = 0

    if method == 'FBA':
        solution = FBA(model, constraints=_constraints, solver=solver)
    elif method == 'pFBA':
        solution = pFBA(model, constraints=_constraints, solver=solver)
    elif method == 'MOMA':
        solution = MOMA(model, reference, constraints=_constraints, solver=solver)
    elif method == 'lMOMA':
        solution = lMOMA(model, reference, constraints=_constraints, solver=solver)
    elif method == 'ROOM':
        solution = ROOM(model, reference, constraints=_constraints, solver=solver)

    return solution


def deletion(model, elements, kind='reactions', method='FBA', reference=None, constraints=None, solver=None):
    """ Generic interface for gene or reaction deletion.
    
    Arguments:
        model (CBModel): model
        elements (list): elements to delete
        kind (str): genes or reactions (default)
        method (str): simulation method: FBA (default), pFMA, MOMA, lMOMA, ROOM
        reference (dict): reference flux distribution for MOMA, lMOMA or ROOM (optional)
        constraints (dict): additional constraints
        solver (Solver): solver instance instantiated with the model, for speed (optional)

    Returns:
        Solution: solution
    """

    if kind == 'genes':
        solution = gene_deletion(model, elements, method, reference, constraints, solver)
    else:
        solution = reaction_deletion(model, elements, method, reference, constraints, solver)

    return solution