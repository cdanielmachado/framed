'''
Created on May 15, 2013

@author: Daniel Machado
'''

from simulation import FBA, MOMA


def gene_deletion(model, genes, method='FBA', reference=None, solver=None, compute_silent_deletions=False):
    """ Simulate deletion of a set of genes. """
    
    inactive_reactions = deleted_genes_to_reactions(model, genes)
    
    if inactive_reactions or compute_silent_deletions:
        solution = reaction_deletion(model, inactive_reactions, method, reference, solver)
    else:
        solution = None
    
    return solution

def deleted_genes_to_reactions(model, genes):
    active_genes = set(model.genes) - set(genes)
    active_reactions = model.eval_GPR(active_genes)
    inactive_reactions = set(model.reactions) - set(active_reactions)

    return inactive_reactions
    

def reaction_deletion(model, reactions, method='FBA', reference=None, solver=None):
    """ Simulate deletion of a set of reactions. """
    
    if method == 'MOMA' and not reference:
        wt_solution = FBA(model, solver=solver)
        reference = wt_solution.values
        
    constraints = {r_id: (0, 0) for r_id in reactions}
    
    if method == 'FBA':
        solution = FBA(model, constraints=constraints, solver=solver)        
    elif method == 'MOMA':
        solution = MOMA(model, reference, constraints=constraints, solver=solver)

    return solution


def deletion(model, elements, kind='reactions', method='FBA', reference=None, solver=None):
    
    if kind == 'genes':    
        solution = gene_deletion(model, elements, method, reference, solver)
    else:
        solution = reaction_deletion(model, elements, method, reference, solver)
        
    return solution