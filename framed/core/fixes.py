'''
Fixes to clean up common models from different sources/groups.

@author: Daniel Machado
'''

def fix_bigg_model(model, boundary_metabolites=True, reversibility=True, bounds=True):
    """ Fix models from BiGG.
    
    Arguments:
        boundary_metabolites : bool -- remove boundary metabolites (ending with '_b') (default: True) 
        reversibility : bool -- make reaction reversibility consistent with the bounds (default: True) 
        clean_bounds : bool -- remove artificially large bounds (unbounded = no bounds) (default: True) 
    """
    if boundary_metabolites:
        remove_boundary_metabolites(model)
    if reversibility:
        fix_reversibility(model)
    if bounds:
        clean_bounds(model)


def remove_boundary_metabolites(model):
    """ Remove remove boundary metabolites (ending with '_b'). """
    
    drains = filter(lambda m_id: m_id.endswith('_b'), model.metabolites)
    model.remove_metabolites(drains)
    
def fix_reversibility(model):
    """ Make reaction reversibility consistent with the bounds. """
    
    for r_id, (lb, _) in model.bounds.items():
        model.reactions[r_id].reversibility = True if lb is None else lb < 0

def clean_bounds(model, threshold=1000):
    """ Remove artificially large bounds (unbounded = no bounds). """
    
    for r_id, (lb, ub) in model.bounds.items():
        model.bounds[r_id] = (lb if lb > -threshold else None, ub if ub < threshold else None)
