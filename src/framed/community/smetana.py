from framed.community.model import Community
from framed.experimental.medium import minimal_medium
from framed.solvers import solver_instance
from framed.solvers.solver import VarType, Status
from framed import Environment

from collections import Counter
from itertools import combinations, chain
from warnings import warn


def species_coupling_score(community, environment=None, min_growth=0.1, n_solutions=100):
    """
    Calculate frequency of community species dependency on each other

    Zelezniak A. et al, Metabolic dependencies drive species co-occurrence in diverse microbial communities (PNAS 2015)

    Args:
        community (Community): microbial community
        environment (Environment): metabolic environment (optional)
        min_growth (float): minimum growth rate (default: 0.1)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)
        n_solutions (int): number of alternative solutions to calculate (default: 100)

    Returns:
        dict: Keys are dependent organisms, values are dictionaries with required organism frequencies
        dict: Extra information
    """
    interacting_community = community.copy(copy_models=False, interacting=True, create_biomass=False,
                                           merge_extracellular_compartments=False)

    if environment:
        environment.apply(interacting_community.merged, inplace=True, warning=False)

    for b in interacting_community.organisms_biomass_reactions.values():
        interacting_community.merged.reactions[b].lb = 0

    solver = solver_instance(interacting_community.merged)

    for org_id, rxns in interacting_community.organisms_reactions.items():
        org_var = 'y_{}'.format(org_id)
        solver.add_variable(org_var, 0, 1, vartype=VarType.BINARY, update_problem=False)

    solver.update()

    bigM = 100
    for org_id, rxns in interacting_community.organisms_reactions.items():
        org_var = 'y_{}'.format(org_id)
        for r_id in rxns:
            if r_id == interacting_community.organisms_biomass_reactions[org_id]:
                continue
            solver.add_constraint('c_{}_lb'.format(r_id), {r_id: 1, org_var: bigM}, '>', 0, update_problem=False)
            solver.add_constraint('c_{}_ub'.format(r_id), {r_id: 1, org_var: -bigM}, '<', 0, update_problem=False)

    solver.update()

    scores = {}
    extras = {'dependencies': {}}

    for org_id, biomass_id in interacting_community.organisms_biomass_reactions.items():
        other_biomasses = {o for o in interacting_community.organisms if o != org_id}
        solver.add_constraint('SMETANA_Biomass', {interacting_community.organisms_biomass_reactions[org_id]: 1}, '>', min_growth)
        objective = {"y_{}".format(o): 1.0 for o in other_biomasses}

        previous_constraints = []
        donors_list = []
        failed = False

        for i in range(n_solutions):
            sol = solver.solve(objective, minimize=True, get_values=True)

            if sol.status != Status.OPTIMAL:
                failed = i == 0
                break

            donors = [o for o in other_biomasses if sol.values["y_{}".format(o)]]
            donors_list.append(donors)

            previous_con = 'iteration_{}'.format(i)
            previous_constraints.append(previous_con)
            previous_sol = {"y_{}".format(o): 1 for o in donors}
            solver.add_constraint(previous_con, previous_sol, '<', len(previous_sol) - 1)

        solver.remove_constraint('SMETANA_Biomass')
        for con in previous_constraints:
            solver.remove_constraint(con)

        if not failed:
            donors_list_n = float(len(donors_list))
            donors_counter = Counter(chain(*donors_list))
            scores[org_id] = {o: donors_counter[o]/donors_list_n for o in other_biomasses}
            extras['dependencies'][org_id] = donors_list
        else:
            scores[org_id] = None
            extras['dependencies'][org_id] = None

    return scores, extras


def metabolite_uptake_score(community, environment=None, min_mass_weight=False, min_growth=0.1,
                            max_uptake=10.0, abstol=1e-6, validate=False, n_solutions=100):
    """
    Calculate frequency of metabolite requirement for species growth

    Zelezniak A. et al, Metabolic dependencies drive species co-occurrence in diverse microbial communities (PNAS 2015)

    Args:
        community (Community): microbial community
        environment (Environment): metabolic environment
        min_mass_weight (bool): Prefer smaller compounds (default: False)
        min_growth (float): minimum growth rate (default: 0.1)
        max_uptake (float): maximum uptake rate (default: 10)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)
        validate (bool): validate solution using FBA (for debugging purposes, default: False)
        n_solutions (int): number of alternative solutions to calculate (default: 100)

    Returns:
        dict: Keys are organism names, values are dictionaries with metabolite frequencies 
        dict: Extra information
    """
    interacting_community = community.copy(copy_models=False, interacting=True, create_biomass=False,
                                           merge_extracellular_compartments=False)

    if environment:
        environment.apply(interacting_community.merged, inplace=True, warning=False)

    solutions = []
    scores = {}
    extras = {'dependencies': {}}
    
    for org_id, exchange_rxns in community.organisms_exchange_reactions.items():
        biomass_reaction = interacting_community.organisms_biomass_reactions[org_id]
        interacting_community.merged.biomass_reaction = biomass_reaction

        medium_list, sol = minimal_medium(interacting_community.merged,
                                   exchange_reactions=exchange_rxns.keys(),
                                   min_mass_weight=min_mass_weight,
                                   min_growth=min_growth,
                                   n_solutions=n_solutions,
                                   max_uptake=max_uptake, validate=validate, abstol=abstol)
        solutions.append(sol)

        if medium_list:
            counter = Counter(chain(*medium_list))
            scores[org_id] = {cnm.original_metabolite: counter[ex] / float(len(medium_list))
                              for ex, cnm in exchange_rxns.items()}
            extras['dependencies'][org_id] = medium_list
        else:
            scores[org_id] = None
            extras['dependencies'][org_id] = None

    return scores, extras


def metabolite_production_score(community, environment=None):
    """
    Discover metabolites which species can produce in community

    Zelezniak A. et al, Metabolic dependencies drive species co-occurrence in diverse microbial communities (PNAS 2015)

    Args:
        community (Community): community object
        environment (Environment): Metabolic environment in which the SMETANA score is colulated
        min_growth (float): minimum growth rate (default: 0.1)
        max_uptake (float): maximum uptake rate (default: 10)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)

    Returns:
        dict: Keys are model names, values are list with produced compounds
        dict: Extra information
    """

    interacting_community = community.copy(copy_models=False, interacting=True, create_biomass=False,
                                           merge_extracellular_compartments=False)
    if environment:
        environment.apply(interacting_community.merged, inplace=True, warning=False)

    solver = solver_instance(interacting_community.merged)

    scores = {}
    extras = {}
    for org_id, exchange_rxns in community.organisms_exchange_reactions.items():
        scores[org_id] = {}
        scores[org_id] = {}

        for r_id, cnm in exchange_rxns.items():
            sol = solver.solve(linear={r_id: 1}, minimize=False, get_values=False)
            scores[org_id][cnm.original_metabolite] = 1 if sol.fobj > 0 else 0
            extras[org_id][cnm.original_metabolite] = sol

    return scores, extras


def mip_score(community, environment=None, min_mass_weight=False, min_growth=1, direction=-1, max_uptake=100, validate=False):
    """
    Implements the metabolic interaction potential (MIP) score as defined in (Zelezniak et al, 2015).

    Args:
        community (Community): microbial community model
        environment (Environment): Metabolic environment in which the SMETANA score is calculated
        direction (int): direction of uptake reactions (negative or positive, default: -1)
        extracellular_id (str): extracellular compartment id
        min_mass_weight (bool): minimize by molecular weight of nutrients (default: False)
        min_growth (float): minimum growth rate (default: 1)
        max_uptake (float): maximum uptake rate (default: 100)
        validate (bool): validate solution using FBA (for debugging purposes, default: False)

    Returns:
        float: MIP score
    """

    interacting_community = community.copy(copy_models=False, interacting=True, merge_extracellular_compartments=False, create_biomass=True)
    noninteracting = community.copy(copy_models=False, interacting=False)
    exch_reactions = set(interacting_community.merged.get_exchange_reactions())

    if environment:
        environment.apply(interacting_community.merged, inplace=True, warning=False)
        environment.apply(noninteracting.merged, inplace=True, warning=False)
        exch_reactions &= set(environment)

    interacting_medium, sol1 = minimal_medium(interacting_community.merged, direction=direction,
                                             exchange_reactions=exch_reactions,
                                             min_mass_weight=min_mass_weight,
                                             min_growth=min_growth,
                                             max_uptake=max_uptake, validate=validate)

    noninteracting_medium, sol2 = minimal_medium(noninteracting.merged,
                                                 exchange_reactions=exch_reactions,
                                                 direction=direction,
                                                 min_mass_weight=min_mass_weight,
                                                 min_growth=min_growth,
                                                 max_uptake=max_uptake, validate=validate)

    if noninteracting_medium is None:
        warn('MIP: Failed to find a valid solution for non-interacting community')
        return None, None
    else:
        score = len(noninteracting_medium) - len(interacting_medium)

    extras = {'noninteracting_medium': noninteracting_medium, 'interacting_medium': interacting_medium,
              'noninteracting_solution': sol2, 'interacting_solution': sol1}

    return score, extras


def mro_score(community, environment=None, direction=-1, min_mass_weight=False, min_growth=1, max_uptake=100, validate=False):
    """
    Implements the metabolic resource overlap (MRO) score as defined in (Zelezniak et al, 2015).

    Args:
        community (Community): microbial community model
        environment (Environment): Metabolic environment in which the SMETANA score is colulated
        direction (int): direction of uptake reactions (negative or positive, default: -1)
        extracellular_id (str): extracellular compartment id
        min_mass_weight (bool): minimize by molecular weight of nutrients (default: False)
        min_growth (float): minimum growth rate (default: 1)
        max_uptake (float): maximum uptake rate (default: 100)

    Returns:
        float: MRO score
    """

    inter_community = community.copy(copy_models=False, interacting=True, merge_extracellular_compartments=False, create_biomass=False)
    indep_community = inter_community.copy(copy_models=False, interacting=False, create_biomass=True)
    exch_reactions = set(inter_community.merged.get_exchange_reactions())

    if environment:
        environment.apply(inter_community.merged, inplace=True, warning=False)
        environment.apply(indep_community.merged, inplace=True, warning=False)
        exch_reactions &= set(environment)

    noninteracting_medium, sol = minimal_medium(indep_community.merged,
                                                    exchange_reactions=exch_reactions,
                                                    direction=direction,
                                                    min_mass_weight=min_mass_weight,
                                                    min_growth=min_growth,
                                                    max_uptake=max_uptake, validate=validate)

    solutions = [sol]

    if sol.status != Status.OPTIMAL:
        warn('MRO: Failed to find a valid solution for non-interacting community')
        return None, None

    # anabiotic environment is limited to non-interacting community minimal media
    noninteracting_exch = set(noninteracting_medium)
    indep_environment = Environment.from_reactions(noninteracting_exch, max_uptake=max_uptake)
    indep_environment.apply(inter_community.merged, inplace=True)

    individual_media = {}

    # TODO: why not optimize the individual models instead ?

    for org_id in inter_community.organisms:
        biomass_reaction = inter_community.organisms_biomass_reactions[org_id]
        inter_community.merged.biomass_reaction = biomass_reaction

        org_noninteracting_exch = inter_community.organisms_exchange_reactions[org_id]

        medium, sol = minimal_medium(inter_community.merged,
                                   exchange_reactions=org_noninteracting_exch,
                                   direction=direction,
                                   min_mass_weight=min_mass_weight,
                                   min_growth=min_growth,
                                   max_uptake=max_uptake, validate=validate)
        solutions.append(sol)

        if sol.status != Status.OPTIMAL:
            warn('MRO: Failed to find a valid solution for: ' + org_id)
            return None, None

        individual_media[org_id] = {org_noninteracting_exch[r].original_metabolite for r in medium}

    pairwise = {(o1, o2): individual_media[o1] & individual_media[o2] for o1, o2 in combinations(community.organisms, 2)}

    numerator = len(individual_media) * sum(map(len, pairwise.values()))
    denominator = float(len(pairwise) * sum(map(len, individual_media.values())))

    score = numerator / denominator if denominator != 0 else None
    extras = {'noninteracting_medium': noninteracting_medium, 'individual_media': individual_media, 
              'pairwise': pairwise, 'solutions': solutions}

    return score, extras
