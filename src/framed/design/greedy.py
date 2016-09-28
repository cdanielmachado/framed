""" This module implements a greedy approach for combinatorial deletions.

@author: Daniel Machado
   
"""

from ..analysis.deletion import deletion
from ..analysis.simulation import FBA
from ..solvers.solver import Status
from ..solvers import solver_instance
from ..analysis.essentiality import essentiality


def greedy_deletion(model, fobj, max_dels, kind='reactions', targets=None, method='FBA', reference=None, min_growth=0.1,
                    pop_size=10, abstol=1e-3):
    """ Generic interface for finding an optimal set of gene or reaction deletions using a greedy approach.
    
    Arguments:
        model : CBModel -- model
        objective : dict (of str to float) -- optimization objective (reaction ids and coefficients)
        max_dels : maximum number of deletions
        kind : str -- genes or reactions (default)
        method : str -- simulation method: FBA (default) or MOMA
        targets : list (of str) -- deletion targets (default: all)
        min_growth : float -- minimum percentage of growth rate to consider a deletion non-letal (default: 0.01)
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
        wt_solution = FBA(model, solver=solver)
        reference = wt_solution.values

    biomass = model.detect_biomass_reaction()
    wt_growth = reference[biomass]
    wt_fval = fobj(reference)

    essential = essentiality(model, kind, min_growth)
    targets = [target for target in targets if target not in essential]

    candidates = [{target} for target in targets]
    solutions = []
    visited = []

    while candidates:
        for del_set in candidates:
            visited.append(del_set)
            solution = deletion(model, del_set, kind, method, reference, solver)
            if solution and solution.status == Status.OPTIMAL:
                fval = fobj(solution.values)
                if fval > wt_fval + abstol and solution.values[biomass] >= min_growth * wt_growth and not _redundant(
                        del_set, fval, solutions, abstol):
                    solutions.append((del_set, fval))
        solutions.sort(key=lambda (_, fval): fval, reverse=True)
        candidates = [del_set | {target} for del_set, _ in solutions[:pop_size] for target in targets
                      if len(del_set | {target}) <= max_dels and (del_set | {target}) not in visited]
    #        print 'solutions: ', len(solutions)
    #        print 'new candidates: ', len(candidates)
    #        print 'best so far: ', solutions[0]

    return solutions


def _redundant(del_set, fval, solutions, abstol):
    redundant = False
    for previous, fval0 in solutions:
        if del_set.issuperset(previous) and fval < fval0 + abstol:
            redundant = True
            break
    return redundant
