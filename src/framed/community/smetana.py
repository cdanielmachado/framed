from framed.experimental.medium import minimal_medium
from framed.model.cbmodel import Environment


def mip_score(community, exchange_pattern="^R_EX_", extracellular_id="C_e", direction=-1, min_mass_weight=False,
              min_growth=1, max_uptake=100):
    """
    Implements the metabolic interaction potential (MIP) score as defined in (Zelezniak et al, 2015).

    Args:
        community (Community): microbial community model

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


def mro_score(community, exchange_pattern="^R_EX_", extracellular_id="C_e", direction=-1, min_mass_weight=False,
              min_growth=1, max_uptake=100):
    """
    Implements the metabolic resource overlap (MRO) score as defined in (Zelezniak et al, 2015).

    Args:
        community (Community): microbial community model

    Returns:
        float: MRO score
    """

    pass

