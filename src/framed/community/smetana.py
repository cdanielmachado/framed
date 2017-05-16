from framed.community.model import Community
from framed.experimental.medium import minimal_medium
from framed.solvers import solver_instance
from framed.solvers.solver import VarType, Status
from framed.model.cbmodel import Environment, CBReaction

from collections import Counter
from itertools import combinations, chain
from random import sample
from warnings import warn
from scipy.misc import comb


from framed.solvers.solver import Status



def species_coupling_score(community, environment, min_growth=1, max_uptake=100):
    if community.merge_extracellular_compartments:
        raise KeyError("For MIP calculation <Community.merge_extracellular_compartments> should be disabled")

    if not community.create_biomass_reaction:
        raise KeyError("For MIP calculation <Community.create_biomass_reaction> should be enabled")

    if not community.interacting:
        raise KeyError("For MIP calculation <Community.interacting> should be enabled")

    interacting_community = community.copy(copy_models=True, interacting=True)
    environment.apply(interacting_community.merged)

    biomass2model = {v: k for k, v in interacting_community.organisms_biomass_reactions.iteritems()}
    n_solutions = int(comb(len(biomass2model), len(biomass2model)-1, exact=False))

    scores = {}
    extras = {'dependencies': {}}
    for biomass_id, org_id in biomass2model.iteritems():
        other_biomass_reactions = set(biomass2model) - set([biomass_id])
        interacting_community.merged.biomass_reaction = biomass_id

        medium_list, solutions = minimal_medium(interacting_community.merged, exchange_reactions=other_biomass_reactions,
                                     direction=1, min_growth=min_growth, n_solutions=n_solutions,
                                     max_uptake=max_uptake, validate=True)
        medium_list_n = float(len(medium_list))

        scores[org_id] = {biomass2model[b_id]: count/medium_list_n for b_id, count in Counter(chain(*medium_list)).iteritems()}
        extras['dependencies'][org_id] = medium_list

    return scores, extras

def species_uptake_score(community, environment, min_mass_weight=True, min_growth=1.0, max_uptake=100.0):
    solutions = []
    interacting_community = community.copy(copy_models=True, interacting=True, create_biomass=True)
    environment.apply(interacting_community.merged, inplace=True)
    n_solutions = int(comb(len(environment), len(environment)-1, exact=False))
    rxn2met = {ex.organism_reaction: ex.original_metabolite
               for org_exchanges in interacting_community.organisms_exchange_reactions.itervalues()
               for ex in org_exchanges.itervalues()}

    scores = {}
    extras = {'dependencies': {}}
    for org_id, exchange_rxns in community.organisms_exchange_reactions.iteritems():
        biomass_reaction = interacting_community.organisms_biomass_reactions[org_id]
        interacting_community.merged.biomass_reaction = biomass_reaction

        medium_list, sol = minimal_medium(interacting_community.merged,
                                   exchange_reactions=exchange_rxns,
                                   direction=-1,
                                   min_mass_weight=min_mass_weight,
                                   min_growth=min_growth,
                                   n_solutions=n_solutions,
                                   max_uptake=max_uptake, validate=True)
        solutions.append(sol)

        medium_list_n = float(len(medium_list))
        scores[org_id] = {rxn2met[rxn_id]: count / medium_list_n
                          for rxn_id, count in Counter(chain(*medium_list)).iteritems()}
        extras['dependencies'][org_id] = medium_list


    return scores, extras


def metabolite_production_score(community, environment, max_uptake=100, abstol=1e-6):
    interacting_community = community.copy(copy_models=True, interacting=True, create_biomass=True)
    environment.apply(interacting_community.merged, inplace=True)

    scores = {}
    extras = {'solutions': {}}

    solver = solver_instance(interacting_community.merged)
    for org_id, exchange_rxns in community.organisms_exchange_reactions.iteritems():
        objective = {}
        for r_id in exchange_rxns:
            objective['y_' + r_id] = 1
            solver.add_variable('y_' + r_id, 0, 1, vartype=VarType.BINARY, update_problem=False)
            solver.add_constraint('c_' + r_id, {r_id: 1, 'y_' + r_id: -max_uptake}, '<', 0, update_problem=False)

        solver.update()
        solution = solver.solve(objective, minimize=False)

        if solution.status != Status.OPTIMAL:
            warn('No solution found')

        scores[org_id] = {r_id: float(solution.values[r_id] < -abstol) for r_id in exchange_rxns}
        extras['solutions'][org_id] = solution

    return scores, extras


def mip_score(community, min_mass_weight=False, min_growth=1, direction=-1, max_uptake=100):
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
    if community.merge_extracellular_compartments:
        raise KeyError("For MIP calculation <Community.merge_extracellular_compartments> should be disabled")

    if not community.create_biomass_reaction:
        raise KeyError("For MIP calculation <Community.create_biomass_reaction> should be enabled")

    if not community.interacting:
        raise KeyError("For MIP calculation <Community.interacting> should be enabled")

    community_medium, sol1 = minimal_medium(community.merged, direction=direction,
                                             min_mass_weight=min_mass_weight,
                                             min_growth=min_growth,
                                             max_uptake=max_uptake, validate=True)

    #
    # Calculate minimal media for non-interacting community
    #
    noninteracting_community = community.copy(copy_models=False, interacting=False)
    noninteracting_medium, sol2 = minimal_medium(noninteracting_community.merged,
                                                 direction=direction,
                                                 min_mass_weight=min_mass_weight,
                                                 min_growth=min_growth,
                                                 max_uptake=max_uptake, validate=True)

    score = len(noninteracting_medium) - len(community_medium)
    extras = {'noninteracting_medium': noninteracting_medium, 'interacting_medium': community_medium,
              'noninteracting_solution': sol2, 'interacting_solution': sol1}

    return score, extras


def mro_score(community, direction=-1, min_mass_weight=False, min_growth=1, max_uptake=100):
    """
    Implements the metabolic resource overlap (MRO) score as defined in (Zelezniak et al, 2015).

    Args:
        community (Community): microbial community model
        direction (int): direction of uptake reactions (negative or positive, default: -1)
        extracellular_id (str): extracellular compartment id
        min_mass_weight (bool): minimize by molecular weight of nutrients (default: False)
        min_growth (float): minimum growth rate (default: 1)
        max_uptake (float): maximum uptake rate (default: 100)

    Returns:
        float: MRO score
    """
    if community.merge_extracellular_compartments:
        raise KeyError("For MRO calculation <Community.merge_extracellular_compartments> should be disabled")

    interacting_community = community.copy(copy_models=True, interacting=True, create_biomass=True)
    # TODO: currently this is not needed. Should we change the algorithm not to touch the initial bounds
    # on exchange reactions?
    for model in community.organisms.values():
        Environment.complete(model, inplace=True)

    noninteracting_community = community.copy(copy_models=False, interacting=False, create_biomass=True)
    noninteracting_medium, sol = minimal_medium(noninteracting_community.merged,
                                                    direction=direction,
                                                    min_mass_weight=min_mass_weight,
                                                    min_growth=min_growth,
                                                    max_uptake=max_uptake, validate=True)

    solutions = [sol]

    if sol.status != Status.OPTIMAL:
        raise RuntimeError('Failed to find a valid solution')

    # Allow uptake of only metabolites that both of species need to reduce number of binary variables
    Environment.empty(interacting_community.merged, inplace=True)
    for r_id in noninteracting_medium:
        interacting_community.merged.set_lower_bound(r_id, -max_uptake)

    individual_media = {}
    for org_id, exchange_rxns in interacting_community.organisms_exchange_reactions.iteritems():
        biomass_reaction = interacting_community.organisms_biomass_reactions[org_id]
        interacting_community.merged.biomass_reaction = biomass_reaction

        medium, sol = minimal_medium(interacting_community.merged,
                                   direction=direction,
                                   min_mass_weight=min_mass_weight,
                                   min_growth=min_growth,
                                   max_uptake=max_uptake, validate=True)
        solutions.append(sol)

        if sol.status != Status.OPTIMAL:
            raise RuntimeError('Failed to find a valid solution')

        individual_media[org_id] = {r_id[:-(1+len(org_id))] for r_id in medium}

    pairwise = {(org1, org2): individual_media[org1] & individual_media[org2]
                for i, org1 in enumerate(community.organisms)
                for j, org2 in enumerate(community.organisms) if i < j}

    numerator = sum(map(len, pairwise.values())) / float(len(pairwise))
    denominator = sum(map(len, individual_media.values())) / float(len(individual_media))

    score = numerator / denominator if denominator != 0 else None
    extras = {'noninteracting_medium': noninteracting_medium, 'individual_media': individual_media, 'pairwise': pairwise,
              'solutions': solutions}

    return score, extras

def score_subcommunities(models, metric, n=None, k=2,  **kwargs):
    """ Apply a given score to subcommunities generated from a list of organisms.

    Args:
        models (list): list of models (CBModel)
        metric (str): metric to apply (currently available: 'MRO', 'MIP')
        n (int): number of samples to generate (optional, all combinations by default)
        k (int): subcommunity size (default: 2)

    Returns:
        dict: score results indexed by subcommunity id

    """

    scores = {}

    metric_map = {
        'MRO': mro_score,
        'MIP': mip_score
    }

    if metric not in metric_map:
        raise RuntimeError('Unsupported metric: {}. Currently supported {}'.format(metric, ','.join(metric_map.keys())))

    subsamples = list(combinations(models, k))

    if n is not None:
        subsamples = sample(subsamples, n)

    for subsample in subsamples:
        comm_id = ','.join(model.id for model in subsample)
        comm = Community(comm_id, subsample, copy_models=False)
        function = metric_map[metric]
        try:
            scores[comm_id], _ = function(comm, **kwargs)
        except:
            warn('{} calculation failed for {}'.format(metric, comm_id))
            continue

    return scores
