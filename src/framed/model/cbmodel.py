import re
import warnings
from collections import OrderedDict

from .model import Model, Metabolite, Reaction, Compartment, AttrOrderedDict
from .parser import ReactionParser


class Gene:
    """ Base class for modeling genes. """

    def __init__(self, elem_id, name=None):
        """
        Arguments:
            elem_id (str): a valid unique identifier
            name (str): common gene name
        """
        self.id = elem_id
        self.name = name if name is not None else elem_id
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name

    # TODO: 5_program_deepcopy.prof
    def copy(self):
        g = Gene(elem_id=self.id, name=self.name)
        if len(self.metadata):
            g.metadata.update(self.metadata)

        return g


class Protein:
    """ Base class for modeling proteins. 
        
        One protein is composed of a list of genes encoding one or more subunits.
    """

    def __init__(self):
        self.genes = []
        self.metadata = OrderedDict()

    def __str__(self):
        return self.to_string()

    def to_string(self, sort=False):

        if len(self.genes) > 1:
            if sort:
                return '(' + ' and '.join(sorted(self.genes)) + ')'
            else:
                return '(' + ' and '.join(self.genes) + ')'
        else:
            return str(self.genes[0])

    # TODO: 5_program_deepcopy.prof
    def copy(self):
        p = Protein()
        p.genes = self.genes
        if len(self.metadata):
            p.metadata.update(self.metadata)

        return p


class GPRAssociation:
    """ Base class for modeling Gene-Protein-Reaction associations. 

        Each GPR association is composed by a list of proteins that can catalyze a reaction.
        Each protein is encoded by one or several genes.
    """

    def __init__(self):
        self.proteins = []
        self.metadata = OrderedDict()

    def __str__(self):
        return self.to_string()

    def to_string(self, sort=False):
        proteins_str = [protein.to_string(sort) for protein in self.proteins]

        if sort:
            proteins_str.sort()

        gpr_str = ' or '.join(proteins_str)

        if len(self.proteins) > 1:
            return '(' + gpr_str + ')'
        else:
            return gpr_str

    def get_genes(self):
        """ Return the set of all genes associated with the respective reaction. """

        genes = set()
        for protein in self.proteins:
            genes |= set(protein.genes)
        return genes

    # TODO: 5_program_deepcopy.prof
    def copy(self):
        gpr = GPRAssociation()
        gpr.proteins = [p.copy() for p in self.proteins]
        if len(self.metadata):
            gpr.metadata.update(self.metadata)

        return gpr


class CBReaction(Reaction):

    def __init__(self, elem_id, name=None, reversible=True, stoichiometry=None, regulators=None,
                 lb=None, ub=None, objective=0, gpr_association=None, is_exchange=None, is_sink=None):
        Reaction.__init__(self, elem_id, name=name, reversible=reversible, stoichiometry=stoichiometry,
                          regulators=regulators, is_exchange=is_exchange, is_sink=is_sink)

        if lb is None and not reversible:
            lb = 0

        self.lb = lb
        self.ub = ub
        self.objective = objective
        self.gpr = gpr_association
        self._bool_function = None

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['_bool_function']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._bool_function = None

    def set_lower_bound(self, value):
        self.lb = value
        self.reversible = bool(value is None or value < 0)

    def set_upper_bound(self, value):
        self.ub = value

    def set_flux_bounds(self, lb, ub):
        self.lb, self.ub = lb, ub
        self.reversible = bool(lb is None or lb < 0)

    def set_gpr_association(self, gpr_association):
        self.gpr = gpr_association

    def set_objective(self, value):
        self.objective = value

    def get_associated_genes(self):
        if self.gpr:
            return self.gpr.get_genes()
        else:
            return []

    def evaluate_gpr(self, active_genes):
        """ Boolean evaluation of the GPR association for a given set of active genes.

        Arguments:
            active_genes (list): list of active genes

        Returns:
            bool: is the reaction active
        """

        if self._bool_function is None:
            self._gpr_to_function()

        return self._bool_function(active_genes)

    def _gpr_to_function(self):
        rule = str(self.gpr) if self.gpr else None
        if not rule:
            rule = 'True'
        else:
            rule = ' ' + rule.replace('(', '( ').replace(')', ' )') + ' '
            for gene in self.get_associated_genes():
                rule = rule.replace(' ' + gene + ' ', ' x[\'' + gene + '\'] ')
        self._bool_function = eval('lambda x: ' + rule)

    def to_string(self, metabolite_names=None):
        """ Print a reaction to a text based representation.

        Arguments:
            metabolite_names (dict): replace metabolite id's with names (optional)

        Returns:
            str: reaction string
        """

        res = Reaction.to_string(self, metabolite_names)

        if self.lb is not None and (self.reversible or self.lb != 0.0) or self.ub is not None:
            res += ' [{}, {}]'.format(self.lb if self.lb is not None else '',
                                      self.ub if self.ub is not None else '')
        if self.objective:
            res += ' @{}'.format(self.objective)

        return res

    # TODO: 5_program_deepcopy.prof
    def copy(self):
        r = CBReaction(elem_id=self.id, name=self.name, reversible=self.reversible,
                        stoichiometry=self.stoichiometry, regulators=self.regulators,
                        is_exchange=self.is_exchange, is_sink=self.is_sink, lb=self.lb, ub=self.ub,
                        objective=self.objective, gpr_association=self.gpr.copy() if self.gpr else self.gpr)
        r._bool_function = self._bool_function
        if len(self.metadata):
            r.metadata.update(self.metadata)

        return r


class CBModel(Model):

    def __init__(self, model_id):
        """
        Arguments:
            model_id (string): a valid unique identifier
        """
        Model.__init__(self, model_id)
        self.genes = AttrOrderedDict()
        self._biomass_reaction = None
        self._biomass_reaction_set = False

    def _clear_temp(self):
        Model._clear_temp(self)
        self._g_r_lookup = None

    @property
    def biomass_reaction(self):
        if not self._biomass_reaction_set:
            self.detect_biomass_reaction()

        return self._biomass_reaction

    @biomass_reaction.setter
    def biomass_reaction(self, r_id):
        if r_id:
            assert r_id in self.reactions, "{} is not a valid reaction id in this model".format(r_id)

        self._biomass_reaction = r_id
        self._biomass_reaction_set = True

    def get_flux_bounds(self, r_id):
        """ Get flux bounds for reaction

        Arguments:
            r_id (str): reaction id

        Returns:
            float: lower bound
            float: upper bound
        """

        reaction = self.reactions[r_id]
        return reaction.lb, reaction.ub

    def set_flux_bounds(self, r_id, lb, ub):
        """ Define flux bounds for one reaction

        Arguments:
            r_id (str): reaction id
            lb (float): lower bound (use None to represent negative infinity)
            ub (float): upper bound (use None to represent positive infinity)
        """
        if r_id in self.reactions:
            self.reactions[r_id].set_flux_bounds(lb, ub)
        else:
            raise KeyError("Reaction '{}' not found".format(r_id))

    def set_lower_bound(self, r_id, lb):
        """ Define lower bound for one reaction

        Arguments:
            r_id (str): reaction id
            lb (float): lower bound (use None to represent negative infinity)
        """
        if r_id in self.reactions:
            self.reactions[r_id].lb = lb
        else:
            raise KeyError("Reaction '{}' not found".format(r_id))

    def set_upper_bound(self, r_id, ub):
        """ Define upper bound for one reaction

        Arguments:
            r_id (str): reaction id
            ub (float): upper bound (use None to represent positive infinity)
        """
        if r_id in self.reactions:
            self.reactions[r_id].ub = ub
        else:
            raise KeyError("Reaction '{}' not found".format(r_id))

    def set_objective(self, coefficients):
        """ Define objective coefficients for a list of reactions

        Arguments:
            coefficients (dict): dictionary of reactions and coefficients

        """
        for r_id, coeff, in coefficients.items():
            self.set_reaction_objective(r_id, coeff)

    def set_reaction_objective(self, r_id, coeff=0):
        """ Define objective coefficient for a single reaction

        Arguments:
            r_id (str): reaction id
            coeff (float): reaction objective (default: 0)
        """
        if r_id in self.reactions:
            self.reactions[r_id].objective = coeff
        else:
            raise KeyError("Reaction '{}' not found".format(r_id))

    def add_reaction(self, reaction, clear_tmp=True):
        """ Add a single reaction to the model.
        If a reaction with the same id exists, it will be replaced.

        Arguments:
            reaction (CBReaction): reaction
        """

        if not isinstance(reaction, CBReaction):
            cbreaction = CBReaction(
                reaction.id,
                name=reaction.name,
                reversible=reaction.reversible,
                stoichiometry=reaction.stoichiometry,
                regulators=reaction.regulators,
                is_exchange=reaction.is_exchange
            )

            cbreaction.metadata = reaction.metadata
            Model.add_reaction(self, cbreaction, clear_tmp=clear_tmp)
        else:
            Model.add_reaction(self, reaction, clear_tmp=clear_tmp)

    def detect_biomass_reaction(self, update=False):
        """ Detects biomass reaction in the model (searches by objective coefficient)

        Returns:
            str: first reaction that matches (or else None)
        """

        if self._biomass_reaction is None or update:
            matches = [r_id for r_id, rxn in self.reactions.items() if rxn.objective]

            if matches:
                self._biomass_reaction = matches[0]
                if len(matches) > 1:
                    w = 'Multiple biomass reactions detected (first selected): {}'.format(" ".join(matches))
                    warnings.warn(w, UserWarning)

                if not re.search("biomass|growth", self._biomass_reaction, re.IGNORECASE):
                    w = "Suspicious biomass reaction '{}' name".format(self._biomass_reaction)
                    warnings.warn(w, UserWarning)
            else:
                w = 'No biomass reaction detected'
                warnings.warn(w, UserWarning)

        return self._biomass_reaction

    def add_gene(self, gene):
        """ Add a gene metabolite to the model.
        If a gene with the same id exists, it will be replaced.

        Arguments:
            gene (Gene): gene
        """
        self.genes[gene.id] = gene

    def remove_gene(self, gene_id):
        """ Remove a gene from the model.

        Arguments:
            str : Gene id
        """
        self.remove_genes([gene_id])

    def remove_genes(self, gene_list):
        """ Remove a set of genes from the model.

        Arguments:
            list : Gene ids
        """

        #TODO: remove genes from GPR associations as well

        for gene_id in gene_list:
            if gene_id in self.genes:
                del self.genes[gene_id]
            else:
                warnings.warn("No such gene '{}'".format(gene_id), RuntimeWarning)

    def set_gpr_association(self, r_id, gpr, add_genes=True):
        """ Set GPR association for a given reaction:

        Arguments:
            r_id (str): reaction id
            gpr (GPRAssociation): GPR association
            add_genes (bool): check if associated genes need to be added to the model
        """

        if r_id in self.reactions:
            self.reactions[r_id].gpr = gpr

            if add_genes and gpr is not None:
                for gene_id in gpr.get_genes():
                    if gene_id not in self.genes:
                        self.add_gene(Gene(gene_id))
        else:
            warnings.warn("No such reaction '{}'".format(r_id), RuntimeWarning)

    def evaluate_gprs(self, active_genes):
        """ Boolean evaluation of the GPR associations for a given set of active genes.

        Arguments:
            active_genes (list): list of active genes

        Returns:
            list: list of active reactions
        """
        genes_state = {gene: gene in active_genes for gene in self.genes}
        return [r_id for r_id, rxn in self.reactions.items() if rxn.evaluate_gpr(genes_state)]

    def add_ratio_constraint(self, r_id_num, r_id_den, ratio):
        """ Add a flux ratio constraint to the model.

        Arguments:
            r_id_num (str): id of the numerator
            r_id_den (str): id of the denominator
            ratio (float): ratio value

        Returns:
            str : identifier of the pseudo-metabolite
        """

        if r_id_num in self.reactions and r_id_den in self.reactions:
            m_id = 'ratio_{}_{}'.format(r_id_num, r_id_den)
            if 'pseudo' not in self.compartments:
                self.add_compartment(Compartment('pseudo'))
            self.add_metabolite(Metabolite(m_id, m_id, compartment='pseudo'))
            self.reactions[r_id_num].stoichiometry[m_id] = 1
            self.reactions[r_id_den].stoichiometry[m_id] = -ratio
            return m_id
        else:
            raise KeyError('Invalid reactions in ratio {}/{}'.format(r_id_num, r_id_den))

    def remove_ratio_constraint(self, r_id_num, r_id_den):
        """ Remove a flux ratio constraint from the model.

        Arguments:
            r_id_num (str): id of the numerator
            r_id_den (str): id of the denominator

        """

        if r_id_num in self.reactions and r_id_den in self.reactions:
            m_id = 'ratio_{}_{}'.format(r_id_num, r_id_den)
            if m_id in self.metabolites:
                self.remove_metabolite(m_id)
            else:
                warnings.warn('No ratio constraint for {}/{}'.format(r_id_num, r_id_den), RuntimeWarning)
        else:
            raise KeyError('Invalid reactions in ratio {}/{}'.format(r_id_num, r_id_den))

    def add_reaction_from_str(self, reaction_str, default_compartment=None, clear_tmp=True):
        """ Parse a reaction from a string and add it to the model.

        Arguments:
            reaction_str (str): string representation a the reaction
            default_compartment (str): default compartment id (optional)

        Notes:
            If the metabolites specified in the reaction are not yet in the model, they will be automatically added.
            You can specify the compartment for new metabolites using the optional argument. However, if you want to
            use multiple compartments you will have to change them manually afterwards.
        """

        if not self._parser:
            self._parser = ReactionParser()

        r_id, reversible, stoichiometry, lb, ub, obj_coeff = \
            self._parser.parse_reaction(reaction_str, kind='cb')

        for m_id in stoichiometry:
            if m_id not in self.metabolites:
                self.add_metabolite(Metabolite(m_id, m_id, compartment=default_compartment), clear_tmp=clear_tmp)

        reaction = CBReaction(r_id, r_id, reversible, stoichiometry, None, lb, ub, obj_coeff)
        self.add_reaction(reaction, clear_tmp=clear_tmp)

        return r_id

    def get_objective(self):
        return {r_id: rxn.objective for r_id, rxn in self.reactions.items() if rxn.objective}

    def gene_to_reaction_lookup(self):
        """ Build a dictionary from genes to associated reactions.

        Returns:
            dict: gene to reaction mapping

        """
        if not self._g_r_lookup:
            self._g_r_lookup = OrderedDict([(g_id, []) for g_id in self.genes])

            for r_id, rxn in self.reactions.items():
                genes = rxn.get_associated_genes()
                for g_id in genes:
                    self._g_r_lookup[g_id].append(r_id)

        return self._g_r_lookup

    def get_reactions_by_gene(self, g_id):
        """ Get a list of reactions associated with a given gene.

        Args:
            g_id (str): gene id

        Returns:
            list: reactions catalyzed by any proteins (or subunits) encoded by this gene
        """
        g_r_lookup = self.gene_to_reaction_lookup()
        return g_r_lookup[g_id]

    def print_objective(self):
        coeffs = ['{:+g} {}'.format(rxn.objective, r_id) for r_id, rxn in self.reactions.items() if rxn.objective]
        return ' '.join(coeffs)


