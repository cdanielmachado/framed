from framed.community.model import Community
from framed.experimental.medium import minimal_medium
from framed.solvers import solver_instance
from framed.solvers.solver import VarType, Status
from framed.model.cbmodel import CBReaction
from framed import FBA, write_model_to_file, Environment

from collections import Counter
from itertools import combinations, chain
from random import sample
from warnings import warn
import itertools


from framed.solvers.solver import Status

class GrowthStatus:
    NOGROWTH = "no growth"
    DEPENDENT = "dependant"
    INDEPENDENT = "independant"


class SmetanaScore(object):
    def __init__(self, donor_organism, receiver_organism, metabolite, metabolite_production_score, metabolite_uptake_score, species_coupling_score):
        self.donor_organism = donor_organism
        self.receiver_organism = receiver_organism
        self.metabolite = metabolite
        self.metabolite_production_score = metabolite_production_score
        self.metabolite_uptake_score = metabolite_uptake_score
        self.species_coupling_score = species_coupling_score
        self.score = metabolite_production_score * metabolite_uptake_score * species_coupling_score

    def __repr__(self):
        return "<{}/{}/{}:{}>".format(self.donor_organism, self.metabolite, self.receiver_organism, self.score)


def calculate_smetana_score(community, scscores, mpscores, muscores, report_zero_scores=False):
    metabolites = {met for mets in mpscores.itervalues() if mets for met in mets}

    scores = []
    for org_donor in community.organisms:
        for org_receiver in community.organisms:
            if org_donor == org_receiver:
                continue

            for met in metabolites:
                mpscore = float(mpscores[org_donor] is not None and met in mpscores[org_donor])
                muscore = muscores[org_receiver].get(met, 0.0) if muscores[org_receiver] else 0.0
                scscore = scscores[org_receiver].get(org_donor, 0.0) if scscores[org_receiver] else 0.0

                s = SmetanaScore(donor_organism=org_donor, receiver_organism=org_receiver, metabolite=met,
                                 metabolite_production_score=mpscore, metabolite_uptake_score=muscore,
                                 species_coupling_score=scscore)

                if not report_zero_scores and not s.score:
                    continue

                scores.append(s)

    return scores

def smetana_score(community, environment, report_zero_scores=False, min_mass_weight=False, min_growth=1, max_uptake=100, abstol=1e-6, validate=False, n_solutions=100):
    """
    SMETANA value scores likelyhood of metabolite exchange from species A to species B
    
    Zelezniak A. et al, Metabolic dependencies drive species co-occurrence in diverse microbial communities (PNAS 2015)

    Args:
        community (Community): community object
        environment (Environment): Metabolic environment in which the SMETANA score is colulated
        report_zero_scores (bool): Report zero SMETANA scores, where there is no metabolite exchange predicted (default false)
        min_mass_weight (bool): minimize by molecular weight of nutrients (default: False)
        min_growth (float): minimum growth rate (default: 1)
        max_uptake (float): maximum uptake rate (default: 100)
        validate (bool): validate solution using FBA (for debugging purposes, default: False)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)
        n_solutions (int): How many unique solutions to calculate for Metabolite Uptake and Species Coupling scores  (default: 100)

    Returns:
        list: Species --> Metabolite --> Species SMETANA scores
        dict: Extra information
    """
    scscores, scextras = species_coupling_score(community, environment, min_growth=min_growth, max_uptake=max_uptake, n_solutions=n_solutions)
    mpscores, mpextras = metabolite_production_score(community, environment)
    muscores, muextras = metabolite_uptake_score(community, environment, min_mass_weight=min_mass_weight, min_growth=min_growth, max_uptake=max_uptake, abstol=abstol, validate=validate, n_solutions=n_solutions)

    scores = calculate_smetana_score(community=community, scscores=scscores, mpscores=mpscores, muscores=muscores, report_zero_scores=report_zero_scores)
    extras = {'status': {},
              "metabolite_production": {'scores': mpscores, 'extras': mpextras},
              "metabolite_uptake": {'scores': muscores, 'extras': muextras},
              "species_coupling": {'scores': scscores, 'extras': scextras}
}
    for org_receiver in community.organisms:
        if scscores[org_receiver] is None:
            extras['status'][org_receiver] = GrowthStatus.NOGROWTH
        elif len(scscores[org_receiver]) == 0:
            extras['status'][org_receiver] = GrowthStatus.INDEPENDENT
        else:
            extras['status'][org_receiver] = GrowthStatus.DEPENDENT

    return scores, extras


def species_coupling_score(community, environment, min_growth=1.0, max_uptake=100, n_solutions=100):
    """
    Calculate frequency of community species dependency on each other

    Zelezniak A. et al, Metabolic dependencies drive species co-occurrence in diverse microbial communities (PNAS 2015)

    Args:
        community (Community): community object
        environment (Environment): Metabolic environment in which the SMETANA score is colulated
        min_growth (float): minimum growth rate (default: 1)
        max_uptake (float): maximum uptake rate (default: 100)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)
        n_solutions (int): How many unique solutions to calculate (default: 100)

    Returns:
        dict: Keys are dependant model names, values are dictionaries with required species frequencies 
        dict: Extra information
    """
    interacting_community = community.copy(copy_models=False, interacting=True, create_biomass=False,
                                           merge_extracellular_compartments=False)
    environment.apply(interacting_community.merged, inplace=True) # other values are copied from previous copy

    for b in interacting_community.organisms_biomass_reactions.itervalues():
        interacting_community.merged.reactions[b].lb = 0

    solver = solver_instance(interacting_community.merged)
    for org_id, rxns in interacting_community.organisms_reactions.iteritems():
        org_var = 'y_{}'.format(org_id)
        solver.add_variable(org_var, 0, 1, vartype=VarType.BINARY, update_problem=False)
        for r_id in rxns:
            lb = min_growth if r_id == interacting_community.organisms_biomass_reactions[org_id] else -max_uptake
            solver.add_constraint('c_{}_lb'.format(r_id), {r_id: 1, org_var: -lb}, '>', 0, update_problem=False)
            solver.add_constraint('c_{}_ub'.format(r_id), {r_id: 1, org_var: -max_uptake}, '<', 0, update_problem=False)

    solver.update()
    scores = {}
    extras = {'dependencies': {}}

    for org_id, biomass_id in interacting_community.organisms_biomass_reactions.iteritems():
        other_biomasses = {o for o in interacting_community.organisms if o != org_id}
        solver.add_constraint('SMETANA_Biomass', {interacting_community.organisms_biomass_reactions[org_id]: 1}, '>', min_growth)
        objective = {"y_{}".format(o): 1.0 for o in other_biomasses}

        previous_constraints = []
        donors_list = []
        for i in xrange(n_solutions):
            sol = solver.solve(objective, minimize=True, get_values=True)

            if sol.status != Status.OPTIMAL:
                if i == 0: donors_list = None # species can not grow
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


        if donors_list:
            donors_list_n = float(len(donors_list))
            donors_counter = Counter(chain(*donors_list))
            scores[org_id] = {o: donors_counter[o]/donors_list_n for o in other_biomasses}
            extras['dependencies'][org_id] = donors_list
        else:
            scores[org_id] = None
            extras['dependencies'][org_id] = donors_list

    return scores, extras


def metabolite_uptake_score(community, environment, min_mass_weight=True, min_growth=1.0, max_uptake=100.0, abstol=1e-6, validate=False, n_solutions=100):
    """
    Calculate frequency of metabolite requirement for species growth

    Zelezniak A. et al, Metabolic dependencies drive species co-occurrence in diverse microbial communities (PNAS 2015)

    Args:
        community (Community): community object
        environment (Environment): Metabolic environment in which the SMETANA score is colulated
        min_mass_weight: Prefer simpler compounds 
        min_growth (float): minimum growth rate (default: 1)
        max_uptake (float): maximum uptake rate (default: 100)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)
        validate (bool): validate solution using FBA (for debugging purposes, default: False)
        n_solutions (int): How many unique solutions to calculate (default: 100)

    Returns:
        dict: Keys are dependant model names, values are dictionaries with required compounds frequencies 
        dict: Extra information
    """
    interacting_community = community.copy(copy_models=False, interacting=True, create_biomass=False,
                                           merge_extracellular_compartments=False)

    environment.apply(interacting_community.merged, inplace=True)
    rxn2met = {ex.organism_reaction: ex.original_metabolite
               for org_exchanges in interacting_community.organisms_exchange_reactions.itervalues()
               for ex in org_exchanges.itervalues()}

    m_rxns = interacting_community.merged.reactions
    media_metabolites = {met
           for exch_id in interacting_community.merged.get_exchange_reactions()
           for met in interacting_community.merged.reactions[exch_id].stoichiometry
           if exch_id in m_rxns and m_rxns[exch_id].lb < 0}

    solutions = []
    scores = {}
    extras = {'dependencies': {}}
    for org_id, exchange_rxns in community.organisms_exchange_reactions.iteritems():
        biomass_reaction = interacting_community.organisms_biomass_reactions[org_id]
        interacting_community.merged.biomass_reaction = biomass_reaction

        # Remove metabolites present in the medium from the list of uptake candidates
        exchange_rxns = {rxn_id for rxn_id, cnm in exchange_rxns.iteritems()
                         if cnm.extracellular_metabolite not in media_metabolites}

        medium_list, sol = minimal_medium(interacting_community.merged,
                                   exchange_reactions=exchange_rxns,
                                   direction=-1,
                                   min_mass_weight=min_mass_weight,
                                   min_growth=min_growth,
                                   n_solutions=n_solutions,
                                   max_uptake=max_uptake, validate=validate, abstol=abstol)
        solutions.append(sol)

        if medium_list:
            medium_list = map(lambda medium: [rxn2met[rxn_id] for rxn_id in medium], medium_list)
            medium_list_n = float(len(medium_list))
            scores[org_id] = {m: count / medium_list_n for m, count in Counter(chain(*medium_list)).iteritems()}
        else:
            scores[org_id] = {} if sol.status == Status.OPTIMAL else None
        extras['dependencies'][org_id] = medium_list

    return scores, extras


def metabolite_production_score(community, environment=None, max_uptake=100, min_growth=1.0, abstol=1e-6):
    """
    Discover metabolites which species can produce in community

    Zelezniak A. et al, Metabolic dependencies drive species co-occurrence in diverse microbial communities (PNAS 2015)

    Args:
        community (Community): community object
        environment (Environment): Metabolic environment in which the SMETANA score is colulated
        min_growth (float): minimum growth rate (default: 1)
        max_uptake (float): maximum uptake rate (default: 100)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)

    Returns:
        dict: Keys are model names, values are list with produced compounds
        dict: Extra information
    """
    interacting_community = community.copy(copy_models=False, interacting=True, create_biomass=False,
                                           merge_extracellular_compartments=False)
    if environment:
        environment.apply(interacting_community.merged, inplace=True)

    reactions = interacting_community.merged.reactions
    rxn2met = {ex.organism_reaction: ex.original_metabolite
           for org_exchanges in interacting_community.organisms_exchange_reactions.itervalues()
           for ex in org_exchanges.itervalues()}
    media_metabolites = {met
           for exch_id in interacting_community.merged.get_exchange_reactions()
           for met in interacting_community.merged.reactions[exch_id].stoichiometry
           if exch_id in reactions and reactions[exch_id].lb < 0}

    solver = solver_instance(interacting_community.merged)

    # Binary constraints that forces biomass production of any model that activates exchanges
    for org_id, exchanges in interacting_community.organisms_exchange_reactions.iteritems():
        org_var = 'y_{}'.format(org_id)
        solver.add_variable(org_var, 0, 1, vartype=VarType.BINARY, update_problem=False)
        for r_id in exchanges:
            if r_id == interacting_community.organisms_biomass_reactions[org_id]:
                lb = min_growth
            else:
                lb = -max_uptake if reactions[r_id].lb is None else reactions[r_id].lb

            ub = max_uptake if reactions[r_id].ub is None else reactions[r_id].ub
            solver.add_constraint('c_{}_lb'.format(r_id), {r_id: 1, org_var: -lb}, '>', 0, update_problem=False)
            solver.add_constraint('c_{}_ub'.format(r_id), {r_id: 1, org_var: -ub}, '<', 0, update_problem=False)

    scores = {}
    for org_id, exchange_rxns in community.organisms_exchange_reactions.iteritems():
        # Remove metabolites present in the medium from the list of product candidates
        exchange_rxns = {rxn_id for rxn_id, cnm in exchange_rxns.iteritems()
                         if cnm.extracellular_metabolite not in media_metabolites}

        org_biomass = community.organisms_biomass_reactions[org_id]
        solver.add_constraint('SMETANA_Biomass', {org_biomass: 1}, '>', min_growth, update_problem=False)
        solver.update()

        org_products = set()
        for i in xrange(30000):
            if not exchange_rxns:
                break

            objective = {r_id: 1.0 for r_id in exchange_rxns}
            solution = solver.solve(objective, minimize=False)
            if solution.status != Status.OPTIMAL:
                if i == 0: org_products = None
                break

            i_products = {r_id for r_id in exchange_rxns if solution.values[r_id] > abstol}
            org_products = org_products.union(i_products)
            exchange_rxns = exchange_rxns - i_products

            if not i_products:
                break

        if org_products is not None:
            scores[org_id] = {rxn2met[r_id] for r_id in org_products}
        else:
            scores[org_id] = None
        solver.remove_constraint('SMETANA_Biomass')

    return scores, {}


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
    # TODO: 1_program_cloneModels.prof
    interacting_community = community.copy(copy_models=False, interacting=True, merge_extracellular_compartments=False, create_biomass=True)
    noninteracting = community.copy(copy_models=False, interacting=False)

    exch_reactions = interacting_community.merged.get_exchange_reactions()

    if environment:
        environment.apply(interacting_community.merged, inplace=True)
        environment.apply(noninteracting.merged, inplace=True)
        exch_reactions = set(exch_reactions) - set(environment)
        
    interacting_medium, sol1 = minimal_medium(interacting_community.merged, direction=direction,
                                             exchange_reactions=exch_reactions,
                                             min_mass_weight=min_mass_weight,
                                             min_growth=min_growth,
                                             max_uptake=max_uptake, validate=validate)

    #
    # Calculate minimal media for non-interacting community
    #
    noninteracting_medium, sol2 = minimal_medium(noninteracting.merged,
                                                 exchange_reactions=exch_reactions,
                                                 direction=direction,
                                                 min_mass_weight=min_mass_weight,
                                                 min_growth=min_growth,
                                                 max_uptake=max_uptake, validate=validate)

    if noninteracting_medium is None:
        score = None
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
    # TODO: 1_program_cloneModels.prof
    inter_community = community.copy(copy_models=False, interacting=True, merge_extracellular_compartments=False, create_biomass=False)
    indep_community = inter_community.copy(copy_models=False, interacting=False, create_biomass=True)

    exch_reactions = set(inter_community.merged.get_exchange_reactions()) - set([inter_community.merged.biomass_reaction])
    if environment:
        environment.apply(inter_community.merged, inplace=True)
        environment.apply(indep_community.merged, inplace=True)
        exch_reactions = exch_reactions - set(environment)

    noninteracting_medium, sol = minimal_medium(indep_community.merged,
                                                    exchange_reactions=exch_reactions,
                                                    direction=direction,
                                                    min_mass_weight=min_mass_weight,
                                                    min_growth=min_growth,
                                                    max_uptake=max_uptake, validate=validate)

    solutions = [sol]

    if sol.status != Status.OPTIMAL:
        raise RuntimeError('Failed to find a valid solution')

    # anabiotic environment is limited to non-interacting community minimal media
    noninteracting_exch = set(noninteracting_medium)

    minimal_medium_set = noninteracting_medium | set(environment)
    indep_environment = Environment.from_reactions(minimal_medium_set, max_uptake=max_uptake)
    indep_environment.apply(inter_community.merged, inplace=True)

    individual_media = {}
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
            raise RuntimeError('Failed to find a valid solution')

        individual_media[org_id] = {org_noninteracting_exch[r].original_metabolite for r in medium}


    pairwise = {(o1, o2): individual_media[o1] & individual_media[o2] for o1, o2 in combinations(community.organisms, 2)}

    numerator = len(individual_media) * sum(map(len, pairwise.values()))
    denominator = float(len(pairwise) * sum(map(len, individual_media.values())))

    score = numerator / denominator if denominator != 0 else None
    extras = {'noninteracting_medium': noninteracting_medium, 'individual_media': individual_media, 
              'pairwise': pairwise, 'solutions': solutions}

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
