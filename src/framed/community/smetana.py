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

def smetana_score(community, environment, report_zero_scores=False, min_mass_weight=True, min_growth=1, max_uptake=100, abstol=1e-6):
    mpscores, mpextras = metabolite_production_score(community, environment)
    muscores, muextras = metabolite_uptake_score(community, environment, min_mass_weight=min_mass_weight, min_growth=min_growth, max_uptake=max_uptake, abstol=abstol)
    scscores, scextras = species_coupling_score(community, environment, min_growth=min_growth, max_uptake=max_uptake, abstol=abstol)
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


def species_coupling_score(community, environment, min_growth=1, max_uptake=100, abstol=1e-6):
    if community.merge_extracellular_compartments:
        raise KeyError("For MIP calculation <Community.merge_extracellular_compartments> should be disabled")

    if not community.create_biomass_reaction:
        raise KeyError("For MIP calculation <Community.create_biomass_reaction> should be enabled")

    if not community.interacting:
        raise KeyError("For MIP calculation <Community.interacting> should be enabled")

    interacting_community = community.copy(copy_models=True, interacting=True, create_biomass=False)
    environment.apply(interacting_community.merged, exclusive=True)

    for b in interacting_community.organisms_biomass_reactions.itervalues():
        interacting_community.merged.reactions[b].lb = 0

    solver = solver_instance(interacting_community.merged)
    for org_id, rxns in interacting_community.organisms_reactions.iteritems():
        org_var = 'y_{}'.format(org_id)
        solver.add_variable(org_var, 0, 1, vartype=VarType.BINARY, update_problem=False)
        for r_id in rxns:
            solver.add_constraint('c_{}_lb'.format(r_id), {r_id: 1, org_var: 999}, '>', -1e-6, update_problem=False)
            solver.add_constraint('c_{}_ub'.format(r_id), {r_id: 1, org_var: -999}, '<', 1e-6, update_problem=False)

    solver.update()

    n_solutions = 30000

    scores = {}
    extras = {'dependencies': {}}

    for org_id, biomass_id in interacting_community.organisms_biomass_reactions.iteritems():
        other_biomasses = {o: b for o, v in interacting_community.organisms.iteritems() if o != org_id}
        solver.add_constraint('SMETANA_Biomass', {interacting_community.organisms_biomass_reactions[org_id]: 1}, '>', 1)
        objective = {"y_{}".format(o): 1.0 for o in other_biomasses}

        previous_constraints = []
        donors_list = []
        for i in xrange(n_solutions):
            sol = solver.solve(objective, minimize=True, get_values=True)
            #solver.problem.write("_{}.lp".format(org_id))
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
            extras['dependencies'][org_id] = []

    return scores, extras

def metabolite_uptake_score(community, environment, min_mass_weight=True, min_growth=1.0, max_uptake=100.0, abstol=1e-6):
    solutions = []
    interacting_community = community.copy(copy_models=True, interacting=True, create_biomass=False)
    environment.apply(interacting_community.merged, inplace=True)
    rxn2met = {ex.organism_reaction: ex.original_metabolite
               for org_exchanges in interacting_community.organisms_exchange_reactions.itervalues()
               for ex in org_exchanges.itervalues()}

    m_rxns = interacting_community.merged.reactions
    media_metabolites = {met
           for exch_id, exch_mets in interacting_community.merged.get_exchange_reactions().iteritems()
           for met in exch_mets
           if exch_id in m_rxns and m_rxns[exch_id].lb < 0}

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
                                   n_solutions=30000,
                                   max_uptake=max_uptake, validate=True, abstol=abstol)
        solutions.append(sol)

        if medium_list:
            medium_list_n = float(len(medium_list))
            scores[org_id] = {rxn2met[rxn_id]: count / medium_list_n
                              for rxn_id, count in Counter(chain(*medium_list)).iteritems()}
        else:
            scores[org_id] = {} if sol.status == Status.OPTIMAL else None
        extras['dependencies'][org_id] = medium_list

    return scores, extras


def metabolite_production_score(community, environment, max_uptake=100, abstol=1e-6):
    interacting_community = community.copy(copy_models=True, interacting=True, create_biomass=False)
    environment.apply(interacting_community.merged, inplace=True)

    reactions = interacting_community.merged.reactions
    rxn2met = {ex.organism_reaction: ex.original_metabolite
           for org_exchanges in interacting_community.organisms_exchange_reactions.itervalues()
           for ex in org_exchanges.itervalues()}
    media_metabolites = {met
           for exch_id, exch_mets in interacting_community.merged.get_exchange_reactions().iteritems()
           for met in exch_mets
           if exch_id in reactions and reactions[exch_id].lb < 0}

    solver = solver_instance(interacting_community.merged)

    # Binary constraints that forces biomass production of any model that activates exchanges
    for org_id, exchanges in interacting_community.organisms_exchange_reactions.iteritems():
        org_var = 'y_{}'.format(org_id)
        solver.add_variable(org_var, 0, 1, vartype=VarType.BINARY, update_problem=False)
        for r_id in exchanges:
            lb = -1000 if reactions[r_id].lb is None else reactions[r_id].lb
            ub = 1000 if reactions[r_id].ub is None else reactions[r_id].ub
            solver.add_constraint('c_{}_lb'.format(r_id), {r_id: 1, org_var: -lb}, '>', 0, update_problem=False)
            solver.add_constraint('c_{}_ub'.format(r_id), {r_id: 1, org_var: -ub}, '<', 0, update_problem=False)

    scores = {}
    for org_id, exchange_rxns in community.organisms_exchange_reactions.iteritems():
        # Remove metabolites present in the medium from the list of product candidates
        exchange_rxns = {rxn_id for rxn_id, cnm in exchange_rxns.iteritems()
                         if cnm.extracellular_metabolite not in media_metabolites}

        org_biomass = community.organisms_biomass_reactions[org_id]
        solver.add_constraint('SMETANA_Biomass', {org_biomass: 1}, '>', 1, update_problem=False)
        solver.update()

        org_products = set()
        for i in xrange(30000):
            if not exchange_rxns:
                break

            objective = {r_id: 1.0 for r_id in exchange_rxns}
            solution = solver.solve(objective, minimize=False)

            if i == 0:
                solver.problem.write("{}_mpscore.lp".format(org_id))

            if solution.status != Status.OPTIMAL:
                if i == 0: org_products = None
                break

            i_products = {r_id for r_id in exchange_rxns if solution.values[r_id] > abstol}
            if not i_products:
                break

            org_products = org_products.union(i_products)
            exchange_rxns = exchange_rxns - i_products

        scores[org_id] = {rxn2met[r_id] for r_id in org_products} if org_products else None
        solver.remove_constraint('SMETANA_Biomass')

    return scores, {}


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
