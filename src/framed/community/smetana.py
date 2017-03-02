from framed.community.model import Community
from framed.experimental.medium import minimal_medium
from framed.model.cbmodel import Environment

from itertools import combinations
from random import sample


def mip_score(community, exchange_pattern="^R_EX_", direction=-1, extracellular_id="C_e", min_mass_weight=False,
              min_growth=1, max_uptake=100):
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

    individual_media = {}

    for org_id, organism in community.organisms.items():

        complete = Environment.complete(organism, exchange_pattern, max_uptake)
        complete.apply(organism)
        exchange_rxns = complete.keys()

        medium, _ = minimal_medium(organism, exchange_rxns,
                                   direction=direction,
                                   min_mass_weight=min_mass_weight,
                                   min_growth=min_growth,
                                   max_uptake=max_uptake)
        individual_media[org_id] = set(medium)

    media_union = reduce(set.__or__, individual_media.values())

    merged_model = community.merge_models(extracellular_id, common_biomass=True)
    complete = Environment.complete(merged_model, exchange_pattern, max_uptake)
    exchange_rxns = complete.keys()

    community_medium, _ = minimal_medium(merged_model, exchange_rxns,
                                         direction=direction,
                                         min_mass_weight=min_mass_weight,
                                         min_growth=min_growth,
                                         max_uptake=max_uptake)

    score = len(media_union) - len(community_medium)

    return score, individual_media, media_union, community_medium


def mro_score(community, exchange_pattern="^R_EX_", direction=-1, min_mass_weight=False,
              min_growth=1, max_uptake=100):
    """
    Implements the metabolic resource overlap (MRO) score as defined in (Zelezniak et al, 2015).

    Args:
        community (Community): microbial community model
        exchange_pattern (str): regex patter for guessing exchange reactions
        direction (int): direction of uptake reactions (negative or positive, default: -1)
        min_mass_weight (bool): minimize by molecular weight of nutrients (default: False)
        min_growth (float): minimum growth rate (default: 1)
        max_uptake (float): maximum uptake rate (default: 100)

    Returns:
        float: MRO score
    """

    individual_media = {}

    for org_id, organism in community.organisms.items():

        complete = Environment.complete(organism, exchange_pattern, max_uptake)
        complete.apply(organism)
        exchange_rxns = complete.keys()

        medium, _ = minimal_medium(organism, exchange_rxns,
                                   direction=direction,
                                   min_mass_weight=min_mass_weight,
                                   min_growth=min_growth,
                                   max_uptake=max_uptake)
        individual_media[org_id] = set(medium)

    combinations = {(org1, org2): individual_media[org1] & individual_media[org2]
                    for i, org1 in enumerate(community.organisms)
                    for j, org2 in enumerate(community.organisms) if i < j}

    numerator = sum(map(len, combinations.values())) / float(len(combinations))
    denominator = sum(map(len, individual_media.values())) / float(len(individual_media))

    score = numerator / denominator

    return score, individual_media, combinations


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
        comm = Community('', subsample, copy=False)
        function = metric_map[metric]
        result = function(comm, **kwargs)
        scores[comm_id] = result[0]

    return scores
