from collections import OrderedDict

from framed.model.cbmodel import CBModel, CBReaction
from framed.model.model import Compartment, Metabolite
from ..model.model import AttrOrderedDict
from warnings import warn
from copy import deepcopy


class CommunityNameMapping(object):
    def __init__(self, original_reaction=None, organism_reaction=None, original_metabolite=None,
                 organism_metabolite=None, extracellular_metabolite=None):
        """
        This class is used to represent mapping between original and merged community model metabolites andreactions
         
        Args:
            original_reaction (str): Name of reaction in original model
            organism_reaction (str): Name of reaction in merged community model 
            original_metabolite (str): Name of metabolite in original model
            organism_metabolite (str): Name of metabolite in merged community model
            extracellular_metabolite (str): Name of "common environment" metabolite in merged community model
        """
        self.original_reaction = original_reaction
        self.organism_reaction = organism_reaction
        self.original_metabolite = original_metabolite
        self.organism_metabolite = organism_metabolite
        self.extracellular_metabolite = extracellular_metabolite

    def __repr__(self):
        return "<orig_m: {}, org_m: {}, ex_m: {}, orig_r: {}, org_r: {}>".format(self.original_metabolite,
                                                                           self.organism_metabolite,
                                                                           self.extracellular_metabolite,
                                                                           self.original_reaction,
                                                                           self.organism_reaction)


class Community(object):
    """
    This class implements a microbial community model.

    It serves as a container for multiple organisms, and can be used to merge multiple single-species models (CBModel)
    into a single multi-species model (CBModel) which is compatible with most types of constraint-based methods.
    """

    """ Merge all organisms into a single multi-species model.

    Args:
        extracellular_compartment_id (str): extracellular compartment id
        merge_extracellular (bool): merge all extracellular compartments into a single one (default: True)
        common_biomass (bool): create a common biomass reaction for the community (default: False)

    Returns:
        CBModel: merged multi-species model

    """
    def __init__(self, community_id, models=None, abundances=None, copy_models=True, extracellular_compartment_id="e",
                 merge_extracellular_compartments=False, create_biomass=True, interacting=True):
        """

        Args:
            community_id (str): community identifier
            models (list): list of models to be merged into single community
            abundances (dict): a dictionary representing each organism initial abundance (keys are organism names)
            copy_models (bool): If true copies for merged models  are created
            extracellular_compartment_id (str): Extracellular compartment id is used when merging extracellular compartments
            merge_extracellular_compartments (bool): Do not create organism specific extracellular compartment
            create_biomass (bool): create biomass reaction with biomass metabolites as reactants
            interacting (bool): If true models will be able to exchange metabolites. Otherwise all produced metabolites will go to sink
        """

        if not interacting and merge_extracellular_compartments:
            raise RuntimeError("Non-interacting models are not supported when merging extracellular compartment")

        self.id = community_id
        self._organisms = AttrOrderedDict(immutable=True)
        self._abundances = AttrOrderedDict(immutable=True)
        self._extracellular_compartment = extracellular_compartment_id # TODO: maybe merge and compartment id arguments should be merged?
        self._merge_extracellular_compartments = merge_extracellular_compartments
        self._create_biomass = create_biomass
        self._merged_model = None
        self._copy_models = copy_models
        self._interacting = interacting
        self._organisms_exchange_reactions = {}
        self._organisms_biomass_reactions = {}

        if models is not None:
            if abundances is not None:
                assert set(m.id for m in models) == set(abundances.keys()), 'Abundance organism ids must match model ids'
            else:
                abundances = {model.id: 1.0 for model in models}

            for model in models:
                self.add_organism(model, abundances[model.id], copy_models)

    @property
    def copy_models(self):
        """
        If true copies for merged models  are created
        
        Returns: bool
        """
        return self._copy_models

    @property
    def create_biomass_reaction(self):
        """
        Create biomass reaction with biomass metabolites as reactants
        
        Returns: bool
        """
        return self._create_biomass

    @create_biomass_reaction.setter
    def create_biomass_reaction(self, value):
        """
        Create biomass reaction with biomass metabolites as reactants
        
        Args:
            value: bool
        """
        self._clear_merged_model()
        self._create_biomass = value

    @property
    def interacting(self):
        """
        If true models will be able to exchange metabolites. Otherwise all produced metabolites will go to sink
        
        Returns: bool
        """
        return self._interacting

    @interacting.setter
    def interacting(self, value):
        """
        If true models will be able to exchange metabolites. Otherwise all produced metabolites will go to sink
        
        Args:
            value: bool
        """
        self._clear_merged_model()
        self._interacting = value

    @property
    def organisms_exchange_reactions(self):
        """
        Returns dictionary containing list of reactions exchanging model metabolites with common environment. 
        Dictionary keys are model ids. Values are dictionaries with keys containing exchange reaction ids and values
        containing various information about these reactions.
        
        Returns: dict
        """
        if not self._merged_model:
            self._merged_model = self._generate_merged_model()

        return self._organisms_exchange_reactions

    @property
    def organisms_biomass_reactions(self):
        """
        Returns dictionary containing reaction exporting biomass to common environment. Keys are model ids, and values
        are reaction ids
        
        Returns: dict
        """
        if not self._merged_model:
            self._merged_model = self._generate_merged_model()

        return self._organisms_biomass_reactions

    @property
    def merge_extracellular_compartments(self):
        """
        Do not create organism specific extracellular compartment
        
        Returns: bool
        """
        return self._merge_extracellular_compartments

    @merge_extracellular_compartments.setter
    def merge_extracellular_compartments(self, value):
        """
        Do not create organism specific extracellular compartment
        
        Args:
            value: bool
        """
        self._clear_merged_model()
        self._merge_extracellular_compartments = value

    @property
    def merged(self):
        """
        Merged models containing every organism as separate compartment
        
        Returns: CBModel
        """
        if not self._merged_model:
            self._merged_model = self._generate_merged_model()

        return self._merged_model

    @property
    def organisms(self):
        """
        Dictionary of organism models which are part of the community. Keys are model ids and values are models
        Returns: dict
        """
        return self._organisms

    @property
    def abundances(self):
        """
        Dictionary of organism abundances. Keys are model ids and values are abundances
        Returns: dict
        """
        return self._abundances

    def __str__(self):
        if set(self._abundances.values()) != {1.0}:
            lines = ['{} ({})'.format(org_id, self._abundances[org_id]) for org_id in self._organisms]
        else:
            lines = self._organisms.keys()
        return '\n'.join(lines)

    def _clear_merged_model(self):
        self._merged_model = None
        self._organisms_exchange_reactions = {}

    def add_organism(self, model, abundance=1.0, copy=True):
        """ Add an organism to this community.

        Args:
            model (CBModel): model of the organism
            abundance (float): abundance of this organism in the community (default: 1.0)
            copy (bool): create a copy of the given model (default: True)

        """
        self._clear_merged_model()

        if model.id in self._organisms:
            warn('Organism {} is already in this community'.format(model.id))
        else:
            if copy:
                model = model.copy()

            self._organisms.__setitem__(model.id, model, force=True)
            self._abundances.__setitem__(model.id, abundance, force=True)

    def remove_organism(self, organism):
        """ Remove an organism from this community

        Args:
            organism (str): organism id

        """
        self._clear_merged_model()

        if organism not in self._organisms:
            warn('Organism {} is not in this community'.format(organism))
        else:
            self._organisms.__delitem__(organism, force=True)
            self._abundances.__delitem__(organism, force=True)

    def _generate_merged_model(self):
        def _id_pattern(object_id, organism_id):
            return "{}_{}".format(object_id, organism_id)
    
        def _name_pattern(object_name, organism_name):
            return "{} ({})".format(object_name, organism_name)

        def _copy_object(obj, org_id, compartment=None):
            new_obj = deepcopy(obj)
            new_obj.id = _id_pattern(obj.id, org_id)
            new_obj.name = _name_pattern(obj.name, org_id)
            if compartment:
                new_obj.compartment = compartment

            return new_obj
        
        models_missing_extracelullar_compartment = [m.id for m in self._organisms.itervalues()
                                                    if self._extracellular_compartment not in m.compartments]
        if models_missing_extracelullar_compartment:
            raise RuntimeError("Extracellular compartment '{}' missing from models: '{}'".format(
                self._extracellular_compartment, "', '".join(models_missing_extracelullar_compartment)))
        
        merged_model = CBModel(self.id)
        merged_model.biomass_reaction = None

        organisms_biomass_metabolites = {}
        
        for org_id, model in self._organisms.items():
            self._organisms_exchange_reactions[org_id] = {}
            self._organisms_biomass_reactions[org_id] = {}


            metabolite_lookup = model.metabolite_reaction_lookup()
            #
            # Create additional extracellular compartment
            #
            if not self._merge_extracellular_compartments:
                pool_compartment = Compartment('pool', 'common pool')
                merged_model.add_compartment(pool_compartment)

            for c_id, comp in model.compartments.items():
                if c_id != self._extracellular_compartment or not self._merge_extracellular_compartments:
                    new_comp = _copy_object(comp, org_id)
                    merged_model.add_compartment(new_comp)
                elif c_id not in merged_model.compartments:
                    merged_model.add_compartment(deepcopy(comp))

            for m_id, met in model.metabolites.items():
                if met.compartment != self._extracellular_compartment or not self._merge_extracellular_compartments:
                    new_met = _copy_object(met, org_id, _id_pattern(met.compartment, org_id))
                    merged_model.add_metabolite(new_met, clear_tmp=False)
                elif m_id not in merged_model.metabolites:
                    merged_model.add_metabolite(deepcopy(met), clear_tmp=False)

                if met.compartment == self._extracellular_compartment and not self._merge_extracellular_compartments:
                    pool_id = _id_pattern(m_id, "pool")
                    if pool_id not in merged_model.metabolites:
                        new_met = _copy_object(met, "pool", "pool")
                        merged_model.add_metabolite(new_met, clear_tmp=False)

                        exch_id = _id_pattern("EX_"+m_id, "pool")
                        exch_name = _name_pattern(met.name, "pool exchange")
                        new_rxn = CBReaction(exch_id, name=exch_name, reversible=True, is_exchange=True)
                        new_rxn.stoichiometry[pool_id] = -1.0
                        merged_model.add_reaction(new_rxn)

            for r_id, rxn in model.reactions.items():
                is_exchange = rxn.is_exchange

                if not is_exchange or not self._merge_extracellular_compartments:
                    new_rxn = _copy_object(rxn, org_id)
                    new_rxn.is_exchange = False

                    for m_id, coeff in rxn.stoichiometry.items():
                        if model.metabolites[m_id].compartment != self._extracellular_compartment or \
                                not self._merge_extracellular_compartments:
                            del new_rxn.stoichiometry[m_id]
                            new_id = _id_pattern(m_id, org_id)
                            new_rxn.stoichiometry[new_id] = coeff

                        if is_exchange and model.metabolites[m_id].compartment == self._extracellular_compartment \
                                and not self._merge_extracellular_compartments:
                            pool_id = _id_pattern(m_id, "pool")
                            new_rxn.stoichiometry[pool_id] = -coeff
                            self._organisms_exchange_reactions[org_id][new_rxn.id] = CommunityNameMapping(
                                organism_reaction=new_rxn.id,
                                original_reaction=r_id,
                                organism_metabolite=new_id,
                                extracellular_metabolite=pool_id,
                                original_metabolite=m_id)

                            if not self.interacting:
                                sink_rxn = CBReaction('Sink_{}'.format(new_id), is_exchange=True, reversible=False)
                                sink_rxn.stoichiometry = {new_id: -1}
                                sink_rxn.lb = 0.0
                                merged_model.add_reaction(sink_rxn)

                    if is_exchange and not self._merge_extracellular_compartments:
                        new_rxn.reversible = True
                        new_rxn.lb = None
                        new_rxn.ub = None if self.interacting else 0.0

                    if rxn.id == model.biomass_reaction:
                        new_rxn.reversible = False

                    if self._create_biomass and rxn.id == model.biomass_reaction:
                        new_rxn.objective = False

                        # Add biomass metabolite to biomass equation
                        m_id = _id_pattern('M_framed_biomass', org_id)
                        name = _name_pattern('Framed biomass', org_id)
                        comp = 'pool' if not self._merge_extracellular_compartments else self._extracellular_compartment
                        biomass_met = Metabolite(m_id, name, comp)
                        merged_model.add_metabolite(biomass_met, clear_tmp=False)
                        new_rxn.stoichiometry[m_id] = 1
                        organisms_biomass_metabolites[org_id] = m_id

                    merged_model.add_reaction(new_rxn)

                else:
                    if is_exchange and self._merge_extracellular_compartments:
                        self._organisms_exchange_reactions[org_id][rxn.id] = CommunityNameMapping(
                            organism_reaction=r_id,
                            original_reaction=r_id,
                            extracellular_metabolite=rxn.stoichiometry.keys()[0],
                            original_metabolite=rxn.stoichiometry.keys()[0],
                            organism_metabolite=None)

                    if r_id in merged_model.reactions:
                        continue

                    new_rxn = deepcopy(rxn)
                    new_rxn.is_exchange = True
                    if rxn.id == model.biomass_reaction and self._create_biomass:
                        new_rxn.reversible = False
                        new_rxn.objective = False

                        m_id = _id_pattern('M_framed_biomass', org_id)
                        name = _name_pattern('Framed biomass', org_id)
                        comp = 'pool' if not self._merge_extracellular_compartments else self._extracellular_compartment
                        biomass_met = Metabolite(m_id, name, comp)
                        merged_model.add_metabolite(biomass_met, clear_tmp=False)
                        new_rxn.stoichiometry[m_id] = 1
                        organisms_biomass_metabolites[org_id] = m_id
                        
                    merged_model.add_reaction(new_rxn)

                if r_id == model.biomass_reaction:
                    self._organisms_biomass_reactions[org_id] = new_rxn.id

        if self._create_biomass:
            biomass_rxn = CBReaction('M_framed_biomass', name="Framed biomass", 
                                     reversible=False, is_exchange=True, objective=1.0)
            for org_biomass in organisms_biomass_metabolites.itervalues():
                biomass_rxn.stoichiometry[org_biomass] = -1

            merged_model.add_reaction(biomass_rxn)
            merged_model.biomass_reaction = biomass_rxn.id

        return merged_model

    def copy(self, copy_models=None, interacting=None, create_biomass=None):
        """
        Copy model object
        Args:
            copy_models (bool): If true copies for merged models  are created
            interacting (bool): If true models will be able to exchange metabolites. Otherwise all produced metabolites will go to sink
            create_biomass (bool): create biomass reaction with biomass metabolites as reactants
        Returns:
            Community
        """
        if copy_models is None:
            copy_models = self._copy_models

        if interacting is None:
            interacting = self._interacting

        if create_biomass is None:
            create_biomass = self._create_biomass

        copy_community = Community(self.id, models=self._organisms.values(), abundances=self._abundances,
                                   copy_models=copy_models, create_biomass=create_biomass,
                                   extracellular_compartment_id=self._extracellular_compartment,
                                   merge_extracellular_compartments=self._merge_extracellular_compartments,
                                   interacting=interacting)

        return copy_community


    def split_fluxes(self, fluxes):
        """ Decompose a flux balance solution of the merged community into organism-specific flux vectors.

        Args:
            fluxes (dict): flux distribution as a single dict

        Returns:
            dict: community flux distribution as a nested dict
        """

        comm_fluxes = OrderedDict()

        for org_id, model in self._organisms.iteritems():
            org_fluxes = [(r_id[:-(1+len(org_id))], val) for r_id, val in fluxes.items() if r_id.endswith(org_id)]
            comm_fluxes[org_id] = OrderedDict(org_fluxes)

        return comm_fluxes