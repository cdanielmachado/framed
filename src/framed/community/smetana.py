from framed.community.model import Community
from framed.experimental.medium import minimal_medium
from framed.model.cbmodel import Environment, CBReaction

from itertools import combinations
from random import sample
from warnings import warn
import re


def mip_score(community, exchange_pattern="^R_EX_", direction=-1, extracellular_id="C_e", min_mass_weight=False,
              min_growth=1, max_uptake=10):
    """
    Implements the metabolic interaction potential (MIP) score as defined in (Zelezniak et al, 2015).

    Args:
        community (Community): microbial community model
        exchange_pattern (str): regex patter for guessing exchange reactions
        direction (int): direction of uptake reactions (negative or positive, default: -1)
        extracellular_id (str): extracellular compartment id
        min_mass_weight (bool): minimize by molecular weight of nutrients (default: False)
        min_growth (float): minimum growth rate (default: 1)
        max_uptake (float): maximum uptake rate (default: 100)

    Returns:
        float: MIP score
    """

    for model in community.organisms.values():
        Environment.complete(model, exchange_pattern=exchange_pattern, inplace=True)

    merged_model = community.merge_models(extracellular_id, merge_extracellular=False, common_biomass=True)
    complete = Environment.complete(merged_model, exchange_pattern='^EX')
    exchange_rxns = complete.keys()

    community_medium, sol1 = minimal_medium(merged_model, exchange_rxns,
                                             direction=direction,
                                             min_mass_weight=min_mass_weight,
                                             min_growth=min_growth,
                                             max_uptake=max_uptake, validate=True)

    block_interactions(merged_model, exchange_pattern, direction, inplace=True)

    noninteracting_medium, sol2 = minimal_medium(merged_model, exchange_rxns,
                                                 direction=direction,
                                                 min_mass_weight=min_mass_weight,
                                                 min_growth=min_growth,
                                                 max_uptake=max_uptake, validate=True)

    score = len(noninteracting_medium) - len(community_medium)
    extras = (noninteracting_medium, community_medium, sol1, sol2)

    return score, extras


def mro_score(community, exchange_pattern="^R_EX_", direction=-1, extracellular_id="C_e", min_mass_weight=False,
              min_growth=1, max_uptake=10):
    """
    Implements the metabolic resource overlap (MRO) score as defined in (Zelezniak et al, 2015).

    Args:
        community (Community): microbial community model
        exchange_pattern (str): regex patter for guessing exchange reactions
        direction (int): direction of uptake reactions (negative or positive, default: -1)
        extracellular_id (str): extracellular compartment id
        min_mass_weight (bool): minimize by molecular weight of nutrients (default: False)
        min_growth (float): minimum growth rate (default: 1)
        max_uptake (float): maximum uptake rate (default: 100)

    Returns:
        float: MRO score
    """

    for model in community.organisms.values():
        Environment.complete(model, exchange_pattern=exchange_pattern, inplace=True)

    non_interacting = community.merge_models(extracellular_id, merge_extracellular=False, common_biomass=True)

    block_interactions(non_interacting, exchange_pattern, direction)
    complete = Environment.complete(non_interacting, exchange_pattern='^EX')
    exchange_rxns = complete.keys()

    noninteracting_medium, sol = minimal_medium(non_interacting, exchange_rxns,
                                                    direction=direction,
                                                    min_mass_weight=min_mass_weight,
                                                    min_growth=min_growth,
                                                    max_uptake=max_uptake, validate=True)

    interacting = community.merge_models(extracellular_id, merge_extracellular=False, common_biomass=False)
    Environment.empty(interacting, exchange_pattern='^EX', inplace=True)

    for r_id in noninteracting_medium:
        interacting.set_lower_bound(r_id, -max_uptake)

    individual_media = {}
    re_pattern = re.compile(exchange_pattern)

    for org_id, organism in community.organisms.items():

        exchange_rxns = [r_id for r_id in interacting.reactions
                         if r_id.rsplit('_', 1)[1] == org_id and re_pattern.search(r_id)]

        biomass = organism.biomass_reaction
        interacting.biomass_reaction = '{}_{}'.format(biomass, org_id)

        medium, _ = minimal_medium(interacting, exchange_rxns,
                                   direction=direction,
                                   min_mass_weight=min_mass_weight,
                                   min_growth=min_growth,
                                   max_uptake=max_uptake, validate=True)

        individual_media[org_id] = {r_id.rsplit('_', 1)[0] for r_id in medium}

    pairwise = {(org1, org2): individual_media[org1] & individual_media[org2]
                for i, org1 in enumerate(community.organisms)
                for j, org2 in enumerate(community.organisms) if i < j}

    numerator = sum(map(len, pairwise.values())) / float(len(pairwise))
    denominator = sum(map(len, individual_media.values())) / float(len(individual_media))

    score = numerator / denominator if denominator != 0 else None
    extras = (noninteracting_medium, individual_media, pairwise)

    return score, extras


def block_interactions(model, exchange_pattern="^R_EX_", direction=-1, inplace=True):

    if not inplace:
        model = model.copy()

    pool_mets = model.get_metabolites_by_compartment('pool')
    re_pattern = re.compile(exchange_pattern)

    for m_id in pool_mets:
        if direction < 0:
            rxns = model.get_metabolite_producers(m_id)
        else:
            rxns = model.get_metabolite_consumers(m_id)

        for r_id in rxns:
            if re_pattern.search(r_id):
                if direction < 0:
                    model.set_upper_bound(r_id, 0)
                    mets = model.reactions[r_id].get_substrates()
                else:
                    model.set_lower_bound(r_id, 0)
                    mets = model.reactions[r_id].get_products()

                for met in mets:
                    sink = CBReaction('Sink_{}'.format(met), reversible=False, stoichiometry={met: -1})
                    model.add_reaction(sink)

    model._clear_temp()

    if not inplace:
        return model


def apply_metric_to_subsamples(models, n, k, metric, **kwargs):

    scores = {}

    metric_map = {
        'MRO': mro_score,
        'MIP': mip_score
    }

    if metric not in metric_map:
        raise RuntimeError('Unsupported metric: {}. Currently supported {}'.format(metric, ','.join(metric_map.keys())))

    subsamples = sample(list(combinations(models, k)), n)

    for subsample in subsamples:
        comm_id = ','.join(model.id for model in subsample)
        comm = Community(comm_id, subsample, copy=False)
        function = metric_map[metric]
        try:
            result = function(comm, **kwargs)
            scores[comm_id] = result[0]
        except:
            warn('{} calculation failed for {}'.format(metric, comm_id))
            continue

    return scores
