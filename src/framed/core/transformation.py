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

    blocked = blocked_reactions(model)
    model.remove_reactions(blocked)
    disconnected = _disconnected_metabolites(model)
    model.remove_metabolites(disconnected)

    return blocked, disconnected


def make_irreversible(model):
    """ Splits all reversible reactions into forward and backward directions.
    For efficiency the given model is converted. To keep a copy use deepcopy first.
    
    Arguments:
        model : StoichiometricModel
        
    Returns:
        dictionary (str to (str, str)): mapping of old reaction ids to splitted reaction ids
    """

    mapping = dict()
    table = model.reaction_metabolite_lookup_table()

    for r_id, reaction in model.reactions.items():
        if reaction.reversible:
            fwd_id = reaction.id + '_f'
            bwd_id = reaction.id + '_b'
            mapping[r_id] = (fwd_id, bwd_id)

            model.add_reaction(Reaction(fwd_id, reaction.name, False))
            model.add_reaction(Reaction(bwd_id, reaction.name, False))

            for m_id, coeff in table[r_id].items():
                model.stoichiometry[(m_id, fwd_id)] = coeff
                model.stoichiometry[(m_id, bwd_id)] = -coeff

            if isinstance(model, CBModel):
                lb, ub = model.bounds[r_id]
                model.set_flux_bounds(fwd_id, 0, ub)
                model.set_flux_bounds(bwd_id, 0, -lb if lb != None else None)
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
        model.add_stoichiometry([(m_id, r_id_num, 1), (m_id, r_id_den, -ratio)])
        return m_id


def _disconnected_metabolites(model):
    m_r_table = model.metabolite_reaction_lookup_table()
    return [m_id for m_id, edges in m_r_table.items() if not edges]