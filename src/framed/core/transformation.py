"""
Module for model transformation operations.

@author: Daniel Machado
   
"""

from model import Metabolite, Reaction
from framed.core.cbmodel import CBModel
from ..analysis.variability import blocked_reactions


def simplify(model):
    """ Removes all blocked reactions in a constraint based model
    
    Arguments:
        model : CBModel
        
    Returns:
        (list (of str), list (of str), list (of str)) : lists of removed reactions, metabolites, and genes
    """

    del_reactions = blocked_reactions(model)
    model.remove_reactions(del_reactions)
    del_metabolites = disconnected_metabolites(model)
    model.remove_metabolites(del_metabolites)
    del_genes = disconnected_genes(model)
    model.remove_genes(del_genes)
    return del_reactions, del_metabolites, del_genes


def make_irreversible(model):
    """ Splits all reversible reactions into forward and backward directions.
    For efficiency the given model is converted. To keep a copy use deepcopy first.
    
    Arguments:
        model : Model (or CBmodel)
        
    Returns:
        dictionary (str to (str, str)): mapping of old reaction ids to splitted reaction ids
    """

    mapping = dict()

    for r_id, reaction in model.reactions.items():
        if reaction.reversible:
            fwd_id = reaction.id + '_f'
            bwd_id = reaction.id + '_b'
            mapping[r_id] = (fwd_id, bwd_id)
            bwd_stoichiometry = [(m_id, -coeff) for m_id, coeff in reaction.stoichiometry.items()]
            model.add_reaction(Reaction(fwd_id, reaction.name, False, reaction.stoichiometry, reaction.regulators))
            model.add_reaction(Reaction(bwd_id, reaction.name, False, bwd_stoichiometry, reaction.regulators))

            if isinstance(model, CBModel):
                lb, ub = model.bounds[r_id]
                lb_fwd = max(0, lb) if lb is not None else 0
                ub_fwd = max(0, ub) if ub is not None else None
                lb_bwd = max(-ub, 0) if ub is not None else 0
                ub_bwd = max(-lb, 0) if lb is not None else None
                model.set_flux_bounds(fwd_id, lb_fwd, ub_fwd)
                model.set_flux_bounds(bwd_id, lb_bwd, ub_bwd)
                obj = model.objective[r_id]
                model.set_reaction_objective(fwd_id, obj if obj >= 0 else 0)
                model.set_reaction_objective(bwd_id, -obj if obj < 0 else 0)
                model.set_gpr_association(fwd_id, model.gpr_associations[r_id])
                model.set_gpr_association(bwd_id, model.gpr_associations[r_id])

            model.remove_reaction(r_id)

    return mapping


def disconnected_metabolites(model):
    m_r_table = model.metabolite_reaction_lookup_table()
    return [m_id for m_id, edges in m_r_table.items() if not edges]


def disconnected_genes(model):
    disconnected = set(model.genes)
    for gpr in model.gpr_associations.values():
        if gpr:
            disconnected -= set(gpr.get_genes())
    return disconnected
