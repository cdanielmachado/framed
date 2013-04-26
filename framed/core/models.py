""" This module defines base classes for metabolic modeling.

TODO: Add pretty printing
TODO: Add support for compartments
TODO: Add self consistency check (e.g: no disconnected components)
TODO: Add gpr parsing
TODO: Add explicit (graph-based) gene-reaction associations
"""


from collections import OrderedDict

class Metabolite:
    """ Base class for metabolite specific details.
    
    Arguments:
        elem_id : String -- a valid unique identifier
        name : String -- common metabolite name
    """
    
    def __init__(self, elem_id, name=None):
        self.id = elem_id
        self.name = name


class Reaction:
    """ Base class for reaction specific details.
    
    Arguments:
        elem_id : String -- a valid unique identifier
        name : String -- common reaction name
    """
    
    def __init__(self, elem_id, name=None, reversible=True):
        self.id = elem_id
        self.name = name
        self.reversible = reversible


class Gene:
    """ Base class for gene specific details.
    
    Arguments:
        elem_id : String -- a valid unique identifier
        name : String -- common gene name
    """
    
    def __init__(self, elem_id, name=None):
        self.id = elem_id
        self.name = name


class StoichiometricModel:
    """ Base class for all metabolic models implemented as a bipartite network.
    Defines the list of metabolites, reactions, and stoichiometry.
    """
    
    def __init__(self, model_id):
        self.id = model_id
        self.metabolites = OrderedDict()
        self.reactions = OrderedDict()
        self.stoichiometry = OrderedDict()
        
    def add_metabolites(self, metabolites):
        for metabolite in metabolites:
            self.add_metabolite(metabolite)

    def add_metabolite(self, metabolite):
        self.metabolites[metabolite.id] = metabolite

    def add_reactions(self, reactions):
        for reaction in reactions:
            self.add_reaction(reaction)

    def add_reaction(self, reaction):
        self.reactions[reaction.id] = reaction
        
    def add_stoichiometry(self, stoichiometry):
        for m_id, r_id, coeff in stoichiometry:
            if m_id in self.metabolites and r_id in self.reactions:
                self.stoichiometry[(m_id, r_id)] = coeff
    
    def remove_metabolites(self, id_list):
        for m_id in id_list:
            del self.metabolites[m_id]
        for (m2_id, r_id) in self.stoichiometry:
            if m2_id in id_list:
                del self.stoichiometry[(m2_id, r_id)]
    
    def remove_metabolite(self, m_id):
        self.remove_metabolites([m_id])

    def remove_reactions(self, id_list):
        for r_id in id_list:
            del self.reactions[r_id]
        for (m_id, r2_id) in self.stoichiometry:
            if r2_id in id_list:
                del self.stoichiometry[(m_id, r2_id)]
    
    def remove_reaction(self, r_id):
        self.remove_reactions([r_id])
    
    def full_matrix(self):
        return [[self.stoichiometry[(m_id, r_id)] if (m_id, r_id) in self.stoichiometry else 0
                 for r_id in self.reactions]
                for m_id in self.metabolites]
    
class ConstraintBasedModel(StoichiometricModel):
    """ Base class for constraint-based models.
    Extends StoichiometricModel with flux bounds
    """
    
    def __init__(self, model_id):
        StoichiometricModel.__init__(self, model_id)
        self.bounds = OrderedDict()
    
    def add_bounds(self, bounds_list):
        for r_id, lb, ub in bounds_list:
            self.add_flux_bounds(r_id, lb, ub)
    
    def add_flux_bounds(self, reaction_id, lb=None, ub=None):
        if reaction_id in self.reactions:
            self.bounds[reaction_id] = (lb, ub)
            
    def add_reaction(self, reaction, lb=None, ub=None):
        StoichiometricModel.add_reaction(self, reaction)
        self.bounds[reaction.id] = (lb, ub)
    
    def remove_reactions(self, id_list):
        StoichiometricModel.remove_reactions(self, id_list)
        for r_id in id_list:
            del self.bounds[r_id]



class GPRConstrainedModel(ConstraintBasedModel):
    """ Base class for constraint-based models with GPR associations.
    Extends ConstraintBasedModel with genes and rules
    """
    
    def __init__(self, model_id):
        ConstraintBasedModel.__init__(self, model_id)
        self.genes = OrderedDict()
        self.rules = OrderedDict()

    def add_genes(self, genes):
        for gene in genes:
            self.add_gene(gene)

    def add_gene(self, gene):
        self.genes[gene.id] = gene
    
    def add_rules(self, rules):
        for r_id, rule in rules:
            self.add_rule(r_id, rule)
    
    def add_rule(self, reaction_id, rule):
        self.rules[reaction_id] = rule

    def add_reaction(self, reaction, lb=None, ub=None, rule=None):
        ConstraintBasedModel.add_reaction(self, reaction, lb, ub)
        self.rules[reaction.id] = rule
    
    def remove_reactions(self, id_list):
        ConstraintBasedModel.remove_reactions(self, id_list)
        for r_id in id_list:
            del self.rules[r_id]
