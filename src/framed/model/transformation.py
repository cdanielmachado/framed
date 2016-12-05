"""
Module for model transformation operations.

Author: Daniel Machado
   
"""

from .model import Reaction
from .cbmodel import CBModel, CBReaction
from ..cobra.variability import blocked_reactions


def simplify(model, clean_compartments=True, inplace=True):
    """ Removes all blocked reactions in a constraint based model
    
    Arguments:
        model (CBModel): model
        clean_compartments (bool): remove empty compartments (default: True)
        inplace (bool): change model in place (default), otherwise create a copy first
        
    Returns:
        CBModel: simplified model (if not in place)
    """

    if not inplace:
        model = model.copy()

    del_reactions = blocked_reactions(model)
    model.remove_reactions(del_reactions)

    del_metabolites = disconnected_metabolites(model)
    model.remove_metabolites(del_metabolites, safe_delete=False)

    del_genes = disconnected_genes(model)
    model.remove_genes(del_genes)

    if clean_compartments:
        empty = empty_compartments(model)
        model.remove_compartments(empty, delete_metabolites=False)

    if not inplace:
        return model


def make_irreversible(model, inplace=True, reactions=None):
    """ Splits all reversible reactions into forward and backward directions.
    
    Arguments:
        model : Model (or CBmodel)
        inplace (bool): change model inplace (default), otherwise create a copy first
        reactions (list) : split only reactions in this list (optional)

    Returns:
        dict: mapping of old reaction ids to splitted reaction ids
    """

    if not inplace:
        model = model.copy()

    if reactions is None:
        reactions = model.reactions.keys()

    mapping = dict()

    for r_id, reaction in model.reactions.items():
        if reaction.reversible and r_id in reactions:
            fwd_id = reaction.id + '_f'
            bwd_id = reaction.id + '_b'
            mapping[r_id] = (fwd_id, bwd_id)
            bwd_stoichiometry = [(m_id, -coeff) for m_id, coeff in reaction.stoichiometry.items()]

            if isinstance(model, CBModel):
                lb, ub = reaction.lb, reaction.ub
                lb_fwd = max(0, lb) if lb is not None else 0
                ub_fwd = max(0, ub) if ub is not None else None
                lb_bwd = max(-ub, 0) if ub is not None else 0
                ub_bwd = max(-lb, 0) if lb is not None else None
                obj = reaction.objective
                obj_fwd = obj if obj >= 0 else 0
                obj_bwd = -obj if obj < 0 else 0
                r_fwd = CBReaction(fwd_id, reaction.name, False, reaction.stoichiometry, reaction.regulators,
                                   lb_fwd, ub_fwd, obj_fwd, reaction.gpr)
                r_bwd = CBReaction(bwd_id, reaction.name, False, bwd_stoichiometry, reaction.regulators,
                                   lb_bwd, ub_bwd, obj_bwd, reaction.gpr)
            else:
                r_fwd = Reaction(fwd_id, reaction.name, False, reaction.stoichiometry, reaction.regulators)
                r_bwd = Reaction(bwd_id, reaction.name, False, bwd_stoichiometry, reaction.regulators)

            model.add_reaction(r_fwd)
            model.add_reaction(r_bwd)
            model.remove_reaction(r_id)

    return mapping


def disconnected_metabolites(model):
    m_r_table = model.metabolite_reaction_lookup()
    return [m_id for m_id, edges in m_r_table.items() if not edges]


def disconnected_genes(model):
    disconnected = set(model.genes)
    for reaction in model.reactions.values():
        disconnected -= set(reaction.get_associated_genes())
    return disconnected


def empty_compartments(model):
    used = [met.compartment for met in model.metabolites.values()]
    empty = set(model.compartments) - set(used)
    return empty