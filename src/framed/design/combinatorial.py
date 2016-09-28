""" This module implements methods for combinatorial deletions.

@author: Daniel Machado

"""

from itertools import combinations
from collections import OrderedDict
from ..analysis.deletion import deletion
from ..analysis.simulation import pFBA
from ..solvers.solver import Status
from ..solvers import solver_instance
from ..analysis.essentiality import essentiality


def combinatorial_gene_deletion(model, objective, max_dels, targets=None, method='FBA', reference=None, min_growth=0.01,
                                abstol=1e-3):
    """ Compute solutions for a set of combinatorial gene deletions.
    
    Arguments:
        model : CBModel -- model
        fobj : function -- optimization objective function
        max_dels : maximum number of deletions
        method : str -- simulation method: FBA (default) or MOMA
        targets : list (of str) -- deletion targets (default: all)
        min_growth : float -- minimum fraction of growth rate to consider a deletion non-letal (default: 0.01)
        abstol : float -- minimum objective function value (default: 1e-4)

    Returns:
        list (of (list of str, float)) -- valid solutions
    """

    return combinatorial_deletion(model, objective, max_dels, 'genes', targets, method, reference, min_growth, abstol)


def combinatorial_reaction_deletion(model, objective, max_dels, targets=None, method='FBA', reference=None,
                                    min_growth=0.01, abstol=1e-3):
    """ Compute solutions for a set of combinatorial reaction deletions.
    
    Arguments:
        model : CBModel -- model
        fobj : function -- optimization objective function
        max_dels : maximum number of deletions
        method : str -- simulation method: FBA (default) or MOMA
        targets : list (of str) -- deletion targets (default: all)
        min_growth : float -- minimum fraction of growth rate to consider a deletion non-letal (default: 0.01)
        abstol : float -- minimum objective function value (default: 1e-4)

    Returns:
        list (of (list of str, float)) -- valid solutions
    """

    return combinatorial_deletion(model, objective, max_dels, 'reactions', targets, method, reference, min_growth,
                                  abstol)


def combinatorial_deletion(model, fobj, max_dels, kind='reactions', targets=None, method='FBA', reference=None,
                           min_growth=0.01, abstol=1e-3):
    """ Generic interface for computing for a set of combinatorial gene or reaction deletions.
    
    Arguments:
        model : CBModel -- model
        fobj : function -- optimization objective function
        max_dels : maximum number of deletions
        kind : str -- genes or reactions (default)
        method : str -- simulation method: FBA (default) or MOMA
        targets : list (of str) -- deletion targets (default: all)
        min_growth : float -- minimum fraction of growth rate to consider a deletion non-letal (default: 0.01)
        abstol : float -- minimum objective function value (default: 1e-4)

    Returns:
        list (of (list of str, float)) -- valid solutions
    """

    if kind == 'genes':
        targets = model.genes if not targets else targets
    else:
        kind = 'reactions'
        targets = model.reactions if not targets else targets

    solver = solver_instance(model)

    if not reference:
        #don't reuse solver here, we don't want the temp variables from pFBA to be persistent
        wt_solution = pFBA(model)
        reference = wt_solution.values

    biomass = model.detect_biomass_reaction()
    wt_growth = reference[biomass]
    wt_fval = fobj(reference)

    # for single deletions there is no gain in computing all essential genes/reactions
    if max_dels > 1:
        essential = essentiality(model, kind, min_growth)
        targets = [target for target in targets if target not in essential]

    #    del_sets = [x for x in combinations(targets, max_dels)]
    del_sets = [del_set for i in range(max_dels) for del_set in combinations(targets, i + 1)]

    solutions = dict()

    for del_set in del_sets:
        solution = deletion(model, del_set, kind, method, reference, solver)

        if solution and solution.status == Status.OPTIMAL:
            fval = fobj(solution.values)
            if fval > wt_fval + abstol and solution.values[biomass] >= min_growth * wt_growth and not _redundant(
                    del_set, fval, solutions, abstol):
                solutions[del_set] = fval

    solutions = OrderedDict(sorted(solutions.items(), key=lambda (_, fval): fval, reverse=True))
    return solutions


def _redundant(del_set, fval, solutions, abstol):
    redundant = False
    del_set = set(del_set)
    for previous, fval0 in solutions.items():
        if del_set.issuperset(set(previous)) and fval < fval0 + abstol:
            redundant = True
            break
    return redundant
