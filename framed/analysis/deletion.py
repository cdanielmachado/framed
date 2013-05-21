'''
Created on May 15, 2013

@author: daniel
'''

from simulation import FBA, MOMA


def gene_deletion(model, genes, method='FBA', reference=None, solver=None):
    """ Simulate deletion of a set of genes. """
    
    active_genes = set(model.genes) - set(genes)
    active_reactions = model.eval_GPR(active_genes)
    inactive_reactions = set(model.reactions) - set(active_reactions)
    
    return reaction_deletion(model, inactive_reactions, method, reference, solver)


def reaction_deletion(model, reactions, method='FBA', reference=None, solver=None):
    """ Simulate deletion of a set of reactions. """
    
    if method == 'MOMA' and not reference:
        wt_solution = FBA(model, solver=solver)
        ref_fluxes = wt_solution.values
        
    constraints = {r_id: (0, 0) for r_id in reactions}
    
    if method == 'FBA':
        solution = FBA(model, constraints=constraints, solver=solver)        
    elif method == 'MOMA':
        solution = MOMA(model, ref_fluxes, constraints=constraints, solver=solver)

    return solution


def deletion(model, elements, kind='reactions', method='FBA', reference=None, solver=None):
    
    if kind == 'genes':    
        solution = gene_deletion(model, elements, method, reference, solver)
    else:
        solution = reaction_deletion(model, elements, method, reference, solver)
        
    return solution