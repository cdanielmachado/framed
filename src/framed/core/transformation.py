'''
Module for model transformation operations.

@author: Daniel Machado

   Copyright 2013 Novo Nordisk Foundation Center for Biosustainability,
   Technical University of Denmark.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
   
'''

from models import Metabolite, Reaction, CBModel
from ..analysis.variability import blocked_reactions


def simplify(model):
    """ Removes all blocked reactions in a constraint based model
    
    Arguments:
        model : ConstraintBasedModel
        
    Returns:
        (list (of str), list (of str)) : lists of removed reactions and metabolites
    """

    del_reactions = blocked_reactions(model)
    model.remove_reactions(del_reactions)
    del_metabolites = _disconnected_metabolites(model)
    model.remove_metabolites(del_metabolites)

    if isinstance(model, CBModel):
        del_genes = _disconnected_genes(model)
        model.remove_genes(del_genes)
        return del_reactions, del_metabolites, del_genes
    else:
        return del_reactions, del_metabolites


def make_irreversible(model):
    """ Splits all reversible reactions into forward and backward directions.
    For efficiency the given model is converted. To keep a copy use deepcopy first.
    
    Arguments:
        model : StoichiometricModel
        
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
            model.add_reaction(Reaction(fwd_id, reaction.name, False, reaction.stoichiometry, reaction.modifiers))
            model.add_reaction(Reaction(bwd_id, reaction.name, False, bwd_stoichiometry, reaction.modifiers))

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
                model.set_rule(fwd_id, model.rules[r_id])
                model.set_rule(bwd_id, model.rules[r_id])

            model.remove_reaction(r_id)

    return mapping


def add_ratio_constraint(model, r_id_num, r_id_den, ratio):
    """ Add a flux ratio constraint to a model.

    Arguments:
        model : StoichiometricModel
        r_id_num : str -- id of the numerator
        r_id_num : str -- id of the denominator
        ratio : float -- ratio value

    Returns:
        str : identifier of the pseudo-metabolite
    """

    if r_id_num in model.reactions and r_id_den in model.reactions:
        m_id = 'ratio_{}_{}'.format(r_id_num, r_id_den)
        model.add_metabolite(Metabolite(m_id))
        model.set_stoichiometry(m_id, r_id_num, 1)
        model.set_stoichiometry(m_id, r_id_den, -ratio)
        return m_id


def _disconnected_metabolites(model):
    m_r_table = model.metabolite_reaction_lookup_table()
    return [m_id for m_id, edges in m_r_table.items() if not edges]

def _disconnected_genes(model):
    connected_genes = reduce(set.__or__, model.reaction_genes.values())
    return set(model.genes) - connected_genes