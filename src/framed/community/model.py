from collections import OrderedDict

from framed.model.cbmodel import CBModel, CBReaction
from framed.model.model import Compartment, Metabolite
from ..model.model import AttrOrderedDict
from warnings import warn
from copy import deepcopy


class Community:
    """
    This class implements a microbial community model.

    It serves as a container for multiple organisms, and can be used to merge multiple single-species models (CBModel)
    into a single multi-species model (CBModel) which is compatible with most types of constraint-based methods.
    """

    def __init__(self, community_id, models=None, abundances=None, copy=True):
        """

        Args:
            community_id (str): community identifier
        """
        self.id = community_id
        self.organisms = AttrOrderedDict()
        self.abundances = {}

        if models is not None:
            if abundances is not None:
                assert set(models.keys()) == set(abundances.keys()), 'Abundance organism ids must match model ids'
            else:
                abundances = {model.id: 1.0 for model in models}

            for model in models:
                self.add_organism(model, abundances[model.id], copy)

    def __str__(self):
        if set(self.abundances.values()) != {1.0}:
            lines = ['{} ({})'.format(org_id, self.abundances[org_id]) for org_id in self.organisms]
        else:
            lines = self.organisms.keys()
        return '\n'.join(lines)

    def add_organism(self, model, abundance=1.0, copy=True):
        """ Add an organism to this community.

        Args:
            model (CBModel): model of the organism
            abundance (float): abundance of this organism in the community (default: 1.0)
            copy (bool): create a copy of the given model (default: True)

        """
        if model.id in self.organisms:
            warn('Organism {} is already in this community'.format(model.id))
        else:
            if copy:
                model = model.copy()

            self.organisms[model.id] = model
            self.abundances[model.id] = abundance

    def remove_organism(self, organism):
        """ Remove an organism from this community

        Args:
            organism (str): organism id

        """
        if organism not in self.organisms:
            warn('Organism {} is not in this community'.format(organism))
        else:
            del self.organisms[organism]
            del self.abundances[organism]

    def merge_models(self, extracellular_id, merge_extracellular=True, common_biomass=False):
        """ Merge all organisms into a single multi-species model.

        Args:
            extracellular_id (str): extracellular compartment id
            merge_extracellular (bool): merge all extracellular compartments into a single one (default: True)
            common_biomass (bool): create a common biomass reaction for the community (default: False)

        Returns:
            CBModel: merged multi-species model

        """

        merged_model = CBModel(self.id)

        for org_id, model in self.organisms.items():

            if not merge_extracellular:
                pool_compartment = Compartment('pool', 'common pool')
                merged_model.add_compartment(pool_compartment)

            for c_id, comp in model.compartments.items():
                if c_id != extracellular_id or not merge_extracellular:
                    new_comp = deepcopy(comp)
                    new_comp.id = '{}_{}'.format(c_id, org_id)
                    new_comp.name = '{} ({})'.format(comp.name, org_id)
                    merged_model.add_compartment(new_comp)
                elif c_id not in merged_model.compartments:
                    merged_model.add_compartment(deepcopy(comp))

            for m_id, met in model.metabolites.items():
                if met.compartment != extracellular_id or not merge_extracellular:
                    new_met = deepcopy(met)
                    new_met.id = '{}_{}'.format(m_id, org_id)
                    new_met.name = '{} ({})'.format(met.name, org_id)
                    new_met.compartment = '{}_{}'.format(met.compartment, org_id)
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
                    new_rxn.id = '{}_{}'.format(r_id, org_id)
                    new_rxn.name = '{} ({})'.format(rxn.name, org_id)

                    for m_id, coeff in rxn.stoichiometry.items():
                        if model.metabolites[m_id].compartment != extracellular_id or not merge_extracellular:
                            del new_rxn.stoichiometry[m_id]
                            new_id = '{}_{}'.format(m_id, org_id)
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

            if common_biomass:
                m_id = 'Biomass_{}'.format(org_id)
                name = 'Biomass ({})'.format(org_id)
                comp = extracellular_id if merge_extracellular else 'pool'
                biomass_met = Metabolite(m_id, name, comp)
                merged_model.add_metabolite(biomass_met)

                biomass_rxn = model.biomass_reaction
                new_id = '{}_{}'.format(biomass_rxn, org_id)
                merged_model.reactions[new_id].stoichiometry[m_id] = 1.0
                merged_model.reactions[new_id].objective = 0

        if common_biomass:

            total = float(sum(self.abundances.values()))
            rel_abundance = {org_id: val / total for org_id, val in self.abundances.items()}

            stoichiometry = {'Biomass_{}'.format(org_id): -rel_abundance[org_id] for org_id in self.organisms}
            community_biomass = CBReaction('Biomass_community', 'Community biomass composition', reversible=False,
                                           stoichiometry=stoichiometry, objective=1.0)
            merged_model.add_reaction(community_biomass)

        return merged_model

    def split_fluxes(self, fluxes):
        """ Decompose a flux balance solution of the merged community into organism-specific flux vectors.

        Args:
            fluxes (dict): flux distribution as a single dict

        Returns:
            dict: community flux distribution as a nested dict
        """

        comm_fluxes = OrderedDict()

        for org_id in self.organisms:
            org_fluxes = [(r_id[:-(1+len(org_id))], val) for r_id, val in fluxes.items() if r_id.endswith(org_id)]
            comm_fluxes[org_id] = OrderedDict(org_fluxes)

        return comm_fluxes

