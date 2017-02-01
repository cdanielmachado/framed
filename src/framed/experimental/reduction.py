""" Implementation of the conjunctive reduction method proposed in:

Machado, D, et al. (2010) "Model transformation of metabolic networks using a Petri net based framework."
International Workshop on Biological Processes & Petri Nets (BioPPN).
"""

from framed.model.model import Reaction
from framed.model.cbmodel import CBModel
from uuid import uuid4
import warnings


def balanced_model_reduction(model, metabolites, fluxes, must_keep=None, max_degree=None, clean_null_fluxes=True,
                             clean_disconnected=True, abstol=1e-9):
    if clean_null_fluxes:
        model.remove_reactions([r_id for r_id, val in fluxes.items() if abs(val) < abstol])

    if max_degree:
        m_r_table = model.metabolite_reaction_lookup()
        metabolites = [m_id for m_id in metabolites if len(m_r_table[m_id]) <= max_degree]

    for m_id in metabolites:
        remove_balanced_metabolite(model, m_id, fluxes, must_keep, abstol)

    if clean_disconnected:
        model.remove_metabolites(_disconnected_metabolites(model), safe_delete=False)


def remove_balanced_metabolite(model, m_id, fluxes, must_keep=None, abstol=1e-9):
    neighbours = _metabolite_neighbours(model, [m_id])

    balance = sum([model.stoichiometry[(m_id, r_id)] * fluxes[r_id] for r_id in neighbours])
    turnover = sum([abs(model.stoichiometry[(m_id, r_id)] * fluxes[r_id]) for r_id in neighbours]) / 2.0

    #   print 'removing {}\t balance {}\t turnover {}'.format(m_id, balance, turnover)

    assert abs(balance) < abstol

    if abs(turnover) > abstol:

        new_neighbours = _reaction_neighbours(model, neighbours)
        new_coeffs = dict()

        for m_id2 in new_neighbours:
            coeff = sum([model.stoichiometry[(m_id2, r_id)] * fluxes[r_id] for r_id in neighbours if
                         (m_id2, r_id) in model.stoichiometry]) / turnover
            flow = sum([abs(model.stoichiometry[(m_id2, r_id)]) * fluxes[r_id] for r_id in neighbours if
                        (m_id2, r_id) in model.stoichiometry]) / 2
            if abs(coeff) > abstol:
                new_coeffs[m_id2] = coeff
            else:
                if must_keep and m_id2 in must_keep and flow > abstol:
                #                    print 'removing {} violated {} turnover {} coeff {} flow {}'.format(m_id, m_id2, turnover, coeff, flow)
                    return

        if new_coeffs:
            new_id = 'R_' + str(uuid4())[:8]
            reversible = all([model.reactions[r_id].reversible for r_id in neighbours])
            model.add_reaction(Reaction(new_id, new_id, reversible))

            if not reversible and isinstance(model, CBModel):
                model.set_lower_bound(new_id, 0)

            for m_id2, coeff in new_coeffs.items():
                model.stoichiometry[(m_id2, new_id)] = coeff

                fluxes[new_id] = turnover

        model.remove_reactions(neighbours)
    else:
        model.remove_reactions([r_id for r_id in neighbours if abs(fluxes[r_id]) < abstol])

    model.remove_metabolite(m_id)


def _metabolite_neighbours(model, metabolites):
    return _get_neighbours(model, metabolites, 'metabolites')


def _reaction_neighbours(model, reactions):
    return _get_neighbours(model, reactions, 'reactions')


def _get_neighbours(model, elements, kind):
    if kind == 'metabolites':
        table = model.metabolite_reaction_lookup()
    elif kind == 'reactions':
        table = model.reaction_metabolite_lookup_table()
    neighbours = []
    for elem in elements:
        for neighbour in table[elem].keys():
            if neighbour not in neighbours:
                neighbours.append(neighbour)
    return neighbours


def _verify_balance(model, metabolites, fluxes, abstol=1e-9):
    m_r_table = model.metabolite_reaction_lookup()

    success = True

    for m_id in metabolites:
        neighbours = m_r_table[m_id]
        balance = sum([coeff * fluxes[r_id] for r_id, coeff in neighbours.items()])
        if abs(balance) > abstol:
            success = False
            warnings.warn('{} balance {}'.format(m_id, balance), UserWarning)
    return success


def _disconnected_metabolites(model):
    m_r_table = model.metabolite_reaction_lookup()
    return [m_id for m_id, edges in m_r_table.items() if not edges]


def _disconnected_reactions(model):
    return [r_id for r_id, rxn in model.reactions.items() if len(rxn.stoichiometry) == 0]
