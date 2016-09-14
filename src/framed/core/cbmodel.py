from collections import OrderedDict

from framed.core.model import Model, Metabolite


class Gene:
    """ Base class for modeling genes. """

    def __init__(self, elem_id, name=None):
        """
        Arguments:
            elem_id : String -- a valid unique identifier
            name : String -- common gene name
        """
        self.id = elem_id
        self.name = name if name is not None else elem_id
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name


class Protein:

    def __init__(self):
        self.genes = []
        self.metadata = OrderedDict()

    def __str__(self):
        if len(self.genes) > 1:
            return '(' + ' and '.join(self.genes) + ')'
        else:
            return str(self.genes[0])


class GPRAssociation:

    def __init__(self):
        self.proteins = []
        self.metadata = OrderedDict()

    def __str__(self):
        gpr_str = ' or '.join(map(str, self.proteins))
        if len(self.proteins) > 1:
            return '(' + gpr_str + ')'
        else:
            return gpr_str

    def get_genes(self):
        genes = set()
        for protein in self.proteins:
            genes |= set(protein.genes)
        return genes


class CBModel(Model):

    def __init__(self, model_id):
        """
        Arguments:
            model_id : String -- a valid unique identifier
        """
        Model.__init__(self, model_id)
        self.bounds = OrderedDict()
        self.objective = OrderedDict()
        self.genes = OrderedDict()
        self.gpr_associations = OrderedDict()
        self.rule_functions = OrderedDict()
        self.biomass_reaction = None

    def set_flux_bounds(self, r_id, lb, ub):
        """ Define flux bounds for one reaction

        Arguments:
            r_id : str -- reaction id
            lb : float -- lower bound (use None to represent negative infinity)
            ub : float -- upper bound (use None to represent positive infinity)
        """
        if r_id in self.reactions:
            self.bounds[r_id] = (lb, ub)
        else:
            print 'No such reaction', r_id

    def set_lower_bound(self, r_id, lb):
        """ Define lower bound for one reaction

        Arguments:
            r_id : str -- reaction id
            lb : float -- lower bound (use None to represent negative infinity)
        """
        if r_id in self.reactions:
            _, ub = self.bounds[r_id]
            self.bounds[r_id] = lb, ub
        else:
            print 'No such reaction', r_id

    def set_upper_bound(self, r_id, ub):
        """ Define upper bound for one reaction

        Arguments:
            r_id : str -- reaction id
            ub : float -- upper bound (use None to represent positive infinity)
        """
        if r_id in self.reactions:
            lb, _ = self.bounds[r_id]
            self.bounds[r_id] = lb, ub
        else:
            print 'No such reaction', r_id

    def set_objective(self, coefficients):
        """ Define objective coefficients for a list of reactions

        """
        for r_id, coeff, in coefficients.items():
            self.set_reaction_objective(r_id, coeff)

    def set_reaction_objective(self, r_id, coeff=0):
        """ Define objective coefficient for a single reaction

        Arguments:
            r_id : str -- reaction id
            coeff : float -- reaction objective (default: 0)
        """
        if r_id in self.reactions:
            self.objective[r_id] = coeff
        else:
            print 'No such reaction', r_id

    def add_reaction(self, reaction, lb=None, ub=None, coeff=0):
        """ Add a single reaction to the model.
        If a reaction with the same id exists, it will be replaced.

        Arguments:
            reaction : Reaction
            lb : float -- lower bound (default: None)
            ub : float -- upper bound (default: None)
            coeff : float -- objective coefficient (default: 0)
        """
        Model.add_reaction(self, reaction)

        if lb == None and not reaction.reversible:
            lb = 0

        self.bounds[reaction.id] = (lb, ub)
        self.objective[reaction.id] = coeff
        self.set_gpr_association(reaction.id, None)

    def remove_reactions(self, id_list):
        """ Remove a list of reactions from the model.
        Also removes all the edges connected to the reactions.

        Arguments:
            id_list : list of str -- reaction ids
        """
        for r_id in id_list:
            if r_id in self.reactions:
                del self.bounds[r_id]
                del self.objective[r_id]
                del self.gpr_associations[r_id]
                del self.rule_functions[r_id]
        Model.remove_reactions(self, id_list)

    def print_reaction(self, r_id, reaction_names=False, metabolite_names=False):
        """ Print a reaction to a text based representation.

        Arguments:
            r_id : str -- reaction id

        Returns:
            str -- reaction string
        """
        res = Model.print_reaction(self, r_id, reaction_names, metabolite_names)
        lb, ub = self.bounds[r_id]
        rev = self.reactions[r_id].reversible
        if lb != None and (rev or lb != 0.0) or ub != None:
            res += ' [{}, {}]'.format(lb if lb != None else '',
                                      ub if ub != None else '')
        coeff = self.objective[r_id]
        if coeff:
            res += ' @{}'.format(coeff)

        return res

    def detect_biomass_reaction(self):
        """ Detects biomass reaction in the model (searches by objective coefficient)

        Returns:
            str -- first reaction that matches (or else None)
        """

        if not self.biomass_reaction:
            matches = [r_id for r_id, coeff in self.objective.items() if coeff]

            if matches:
                self.biomass_reaction = matches[0]
                if len(matches) == 1:
                    print 'Biomass reaction detected:', self.biomass_reaction
                else:
                    print 'Multiple biomass reactions detected (first selected):', " ".join(matches)
            else:
                print 'No biomass reaction detected.'

        return self.biomass_reaction

    def add_gene(self, gene):
        """ Add a gene metabolite to the model.
        If a gene with the same id exists, it will be replaced.

        Arguments:
            gene : Gene
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
            TODO: When switching to a GPR class representation, make sure to clean up rules as well

        Arguments:
            list of str : Gene ids
        """

        for gene_id in gene_list:
            if gene_id in self.genes:
                del self.genes[gene_id]
            else:
                print 'No such gene', gene_id

    def set_gpr_association(self, r_id, gpr):

        if r_id in self.reactions:
            self.gpr_associations[r_id] = gpr
            self.rule_functions[r_id] = self._rule_to_function(gpr)
        else:
            print 'No such reaction', r_id

    def eval_GPR(self, active_genes):
        """ Evaluate the GPR associations.

        Arguments:
            active_genes : list (of str) -- set of active genes

        Returns:
            list (of str) -- set of active reactions
        """
        genes_state = {gene: gene in active_genes for gene in self.genes}
        return [r_id for r_id, f in self.rule_functions.items() if f(genes_state)]

    def _rule_to_function(self, gpr):
        rule = str(gpr) if gpr else None
        if not rule:
            rule = 'True'
        else:
            rule = ' ' + rule.replace('(', '( ').replace(')', ' )') + ' '
            for gene in self.genes:
                rule = rule.replace(' ' + gene + ' ', ' x[\'' + gene + '\'] ')
        return eval('lambda x: ' + rule)

    def add_ratio_constraint(self, r_id_num, r_id_den, ratio):
        """ Add a flux ratio constraint to the model.

        Arguments:
            r_id_num : str -- id of the numerator
            r_id_num : str -- id of the denominator
            ratio : float -- ratio value

        Returns:
            str : identifier of the pseudo-metabolite
        """

        if r_id_num in self.reactions and r_id_den in self.reactions:
            m_id = 'ratio_{}_{}'.format(r_id_num, r_id_den)
            self.add_metabolite(Metabolite(m_id))
            self.set_stoichiometry(m_id, r_id_num, 1)
            self.set_stoichiometry(m_id, r_id_den, -ratio)
            return m_id
        else:
            print 'Invalid reactions.'

    def remove_ratio_constraint(self, r_id_num, r_id_den):
        """ Remove a flux ratio constraint from the model.

        Arguments:
            r_id_num : str -- id of the numerator
            r_id_num : str -- id of the denominator

        """

        if r_id_num in self.reactions and r_id_den in self.reactions:
            m_id = 'ratio_{}_{}'.format(r_id_num, r_id_den)
            if m_id in self.metabolites:
                self.remove_metabolite(m_id)
            else:
                print 'No ratio constraint for {}/{}'.format(r_id_num, r_id_den)
        else:
            print 'Invalid reactions.'