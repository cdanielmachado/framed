from framed.model.cbmodel import CBModel, CBReaction
from framed.model.model import Compartment
from ..model.model import AttrOrderedDict
from warnings import warn
from copy import deepcopy


class Community:

    def __init__(self, community_id):
        self.id = community_id
        self.organisms = AttrOrderedDict()

    def add_organism(self, model, copy=True):
        if model.id in self.organisms:
            warn('Organism {} is already in this community'.format(model.id))
        else:
            if copy:
                model = model.copy()

            self.organisms[model.id] = model

    def remove_organism(self, organism):
        if organism not in self.organisms:
            warn('Organism {} is not in this community'.format(organism))
        else:
            del self.organisms[organism]

    def merge_models(self, extracellular_id, merge_extracellular=True):

        merged_model = CBModel(self.id)

        for organism, model in self.organisms.items():

            if not merge_extracellular:
                pool_compartment = Compartment('pool', 'common pool')
                merged_model.add_compartment(pool_compartment)

            for c_id, comp in model.compartments.items():
                if c_id != extracellular_id or not merge_extracellular:
                    new_comp = deepcopy(comp)
                    new_comp.id = '{}_{}'.format(c_id, organism)
                    new_comp.name = '{} ({})'.format(comp.name, organism)
                    merged_model.add_compartment(new_comp)
                elif c_id not in merged_model.compartments:
                    merged_model.add_compartment(deepcopy(comp))

            for m_id, met in model.metabolites.items():
                if met.compartment != extracellular_id or not merge_extracellular:
                    new_met = deepcopy(met)
                    new_met.id = '{}_{}'.format(m_id, organism)
                    new_met.name = '{} ({})'.format(met.name, organism)
                    new_met.compartment = '{}_{}'.format(met.compartment, organism)
                    merged_model.add_metabolite(new_met)
                elif m_id not in merged_model.metabolites:
                    merged_model.add_metabolite(deepcopy(met))

                if met.compartment == extracellular_id and not merge_extracellular:
                    pool_id = '{}_pool'.format(m_id)
                    if pool_id not in merged_model.metabolites:
                        new_met = deepcopy(met)
                        new_met.id = pool_id
                        new_met.name = '{} (pool)'.format(met.name)
                        new_met.compartment = 'pool'
                        merged_model.add_metabolite(new_met)

                        exch_id = 'EX_{}_pool'.format(m_id)
                        exch_name = '{} (pool exchange)'.format(met.name)
                        new_rxn = CBReaction(exch_id, exch_name, False)
                        new_rxn.stoichiometry[pool_id] = -1.0
                        merged_model.add_reaction(new_rxn)

            for r_id, rxn in model.reactions.items():

                compartments = {model.metabolites[m_id].compartment for m_id in rxn.stoichiometry}
                is_exchange = compartments == {extracellular_id}

                if not is_exchange or not merge_extracellular:
                    new_rxn = deepcopy(rxn)
                    new_rxn.id = '{}_{}'.format(r_id, organism)
                    new_rxn.name = '{} ({})'.format(rxn.name, organism)

                    for m_id, coeff in rxn.stoichiometry.items():
                        if model.metabolites[m_id].compartment != extracellular_id or not merge_extracellular:
                            del new_rxn.stoichiometry[m_id]
                            new_id = '{}_{}'.format(m_id, organism)
                            new_rxn.stoichiometry[new_id] = coeff
                        if (is_exchange and model.metabolites[m_id].compartment == extracellular_id
                                and not merge_extracellular):
                            pool_id = '{}_pool'.format(m_id)
                            new_rxn.stoichiometry[pool_id] = -coeff

                    if is_exchange and not merge_extracellular:
                        new_rxn.reversible = True
                        new_rxn.lb, new_rxn.ub = None, None

                    merged_model.add_reaction(new_rxn)

                elif r_id not in merged_model.reactions:
                    merged_model.add_reaction(deepcopy(rxn))

        return merged_model




