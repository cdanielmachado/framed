from collections import OrderedDict

from framed.core.model import Model


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
        self.reaction_genes = OrderedDict()
        self.rules = OrderedDict()
        self.rule_functions = OrderedDict()
        self.biomass_reaction = None

    def set_multiple_bounds(self, bounds):
        """ Define flux bounds for a set of reactions

        """
        for r_id, (lb, ub) in bounds:
            self.set_flux_bounds(r_id, lb, ub)

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
        for r_id, coeff, in coefficients:
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
        self.set_rule(reaction.id, '')
        self.reaction_genes[reaction.id] = set()



    def remove_reactions(self, id_list):
        """ Remove a list of reactions from the model.
        Also removes all the edges connected to the reactions.

        Arguments:
            id_list : list of str -- reaction ids
        """
        Model.remove_reactions(self, id_list)
        for r_id in id_list:
            del self.bounds[r_id]
            del self.objective[r_id]
            del self.rules[r_id]
            del self.rule_functions[r_id]
            del self.reaction_genes[r_id]

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


    def add_genes(self, genes):
        """ Add a list of genes to the model.

        Arguments:
            genes : list of Gene
        """
        for gene in genes:
            self.add_gene(gene)

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
            del self.genes[gene_id]
        else:
            print 'No such gene', gene_id


    def set_rules(self, rules):
        """ Define GPR association rules for a set of reactions

        Arguments:
            rules : list (of (str, str)) -- reaction id, rule
        """
        for r_id, rule in rules:
            self.set_rule(r_id, rule)

    def set_rule(self, r_id, rule):
        """ Define GPR association rule for one reaction

        Arguments:
            r_id : str -- reaction id
            rule : str -- GPR association rule
        """
        if r_id in self.reactions:
            self.rules[r_id] = rule
            self.rule_functions[r_id] = self._rule_to_function(rule)
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


    def _rule_to_function(self, rule):
        if not rule:
            rule = 'True'
        else:
            rule = ' ' + rule.replace('(', '( ').replace(')', ' )') + ' '
            for gene in self.genes:
                rule = rule.replace(' ' + gene + ' ', ' x[\'' + gene + '\'] ')
        return eval('lambda x: ' + rule)


class Gene:
    """ Base class for modeling genes. """

    def __init__(self, elem_id, name=None):
        """
        Arguments:
            elem_id : String -- a valid unique identifier
            name : String -- common gene name
        """
        self.id = elem_id
        self.name = name

    def __str__(self):
        return self.name if self.name else self.id