""" This module defines base classes for metabolic modeling.

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
    
    def __init__(self, elem_id, name=None, compartment=None):
        self.id = elem_id
        self.name = name
        self.compartment = compartment


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

class Compartment:
    """ Base class for compartments.
    
    Arguments:
        elem_id : String -- a valid unique identifier
        name : String -- compartment name
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
        self.compartments = OrderedDict()
        self.stoichiometry = OrderedDict()
        
    def add_metabolites(self, metabolites):
        for metabolite in metabolites:
            self.add_metabolite(metabolite)

    def add_metabolite(self, metabolite):
        if metabolite.compartment in self.compartments or not metabolite.compartment:
            self.metabolites[metabolite.id] = metabolite

    def add_reactions(self, reactions):
        for reaction in reactions:
            self.add_reaction(reaction)

    def add_reaction(self, reaction):
        self.reactions[reaction.id] = reaction
    
    def add_compartments(self, compartments):
        for compartment in compartments:
            self.add_compartment(compartment)

    def add_compartment(self, compartment):
        self.compartments[compartment.id] = compartment
        
        
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
        
        
    def remove_compartment(self, c_id, delete_metabolites=True):
        if c_id in self.compartments:
            del self.compartments[c_id]
            
            if delete_metabolites:
                self.remove_metabolites([m_id for m_id, metabolite in self.metabolites 
                                         if metabolite.compartment == c_id]) 
        

    def metabolite_reaction_lookup_table(self):
        table = {m_id : dict() for m_id in self.metabolites} 
    
        for (m_id, r_id), coeff in self.stoichiometry.items():
            table[m_id][r_id] = coeff
        
        return table
 
   
    def reaction_metabolite_lookup_table(self):
        table = {r_id : dict() for r_id in self.reactions} 
    
        for (m_id, r_id), coeff in self.stoichiometry.items():
            table[r_id][m_id] = coeff
        
        return table
 
    
    def stoichiometric_matrix(self):
        return [[self.stoichiometry[(m_id, r_id)] if (m_id, r_id) in self.stoichiometry else 0
                 for r_id in self.reactions]
                for m_id in self.metabolites]        
    
    
    def __repr__(self):
        return '\n'.join([self.print_reaction(r_id) for r_id in self.reactions])
    
    def print_reaction(self, r_id):
        res = r_id + ': '
        res += ' + '.join([m_id if coeff == -1.0 else str(-coeff) + ' ' + m_id
                         for (m_id, r_id2), coeff in self.stoichiometry.items()
                         if r_id2 == r_id and coeff < 0])
        res += ' <-> ' if self.reactions[r_id].reversible else ' --> '
        res += ' + '.join([m_id if coeff == 1.0 else str(coeff) + ' ' + m_id
                         for (m_id, r_id2), coeff in self.stoichiometry.items()
                         if r_id2 == r_id and coeff > 0])
        return res   
            
             
class ConstraintBasedModel(StoichiometricModel):
    """ Base class for constraint-based models.
    Extends StoichiometricModel with flux bounds
    """
    
    def __init__(self, model_id):
        StoichiometricModel.__init__(self, model_id)
        self.bounds = OrderedDict()
    
    def set_bounds(self, bounds_list):
        for r_id, lb, ub in bounds_list:
            self.set_flux_bounds(r_id, lb, ub)
    
    def set_flux_bounds(self, reaction_id, lb, ub):
        if reaction_id in self.reactions:
            self.bounds[reaction_id] = (lb, ub)
            
    def set_lower_bound(self, reaction_id, lb):
        if reaction_id in self.reactions:
            self.bounds[reaction_id][0] = lb
                
    def set_upper_bound(self, reaction_id, ub):
        if reaction_id in self.reactions:
            self.bounds[reaction_id][1] = ub

    def add_reaction(self, reaction, lb=None, ub=None):
        StoichiometricModel.add_reaction(self, reaction)
        self.bounds[reaction.id] = (lb, ub)
    
    def remove_reactions(self, id_list):
        StoichiometricModel.remove_reactions(self, id_list)
        for r_id in id_list:
            del self.bounds[r_id]

    def print_reaction(self, r_id):
        res = StoichiometricModel.print_reaction(self, r_id)
        lb, ub = self.bounds[r_id]
        rev = self.reactions[r_id].reversible
        if lb != None and (rev or lb != 0.0) or ub != None:
            res += ' [{}, {}]'.format(lb if lb != None else '',
                                      ub if ub != None else '')
        return res 

    def detect_biomass_reaction(self):
        matches = [r_id for r_id in self.reactions if 'biomass' in r_id.lower()]
        return matches[0] if matches else None

class GPRConstrainedModel(ConstraintBasedModel):
    """ Base class for constraint-based models with GPR associations.
    Extends ConstraintBasedModel with genes and rules
    """
    
    def __init__(self, model_id):
        ConstraintBasedModel.__init__(self, model_id)
        self.genes = OrderedDict()
        self.rules = OrderedDict()
        self.rule_functions = OrderedDict()

    def add_genes(self, genes):
        for gene in genes:
            self.add_gene(gene)

    def add_gene(self, gene):
        self.genes[gene.id] = gene

    def add_reaction(self, reaction, lb=None, ub=None, rule=None):
        ConstraintBasedModel.add_reaction(self, reaction, lb, ub)
        self.set_rule(reaction.id, rule)
    
    def remove_reactions(self, id_list):
        ConstraintBasedModel.remove_reactions(self, id_list)
        for r_id in id_list:
            del self.rules[r_id]
            del self.rule_functions[r_id]
    
    def set_rules(self, rules):
        for r_id, rule in rules:
            self.set_rule(r_id, rule)
    
    def set_rule(self, r_id, rule):
        if r_id in self.reactions:
            self.rules[r_id] = rule
            self.rule_functions[r_id] = self._rule_to_function(rule)
    
    def _rule_to_function(self, rule):
        if not rule:
            rule = 'True'
        else:
            for gene in self.genes:
                rule = rule.replace(gene, 'x[\'' + gene + '\']')
        return eval('lambda x: ' + rule)
    
    def eval_GPR(self, active_genes):
        genes_state = {gene: gene in active_genes for gene in self.genes}
        return [r_id for r_id, f in self.rule_functions.items() if f(genes_state)]
