""" This module implements methods for reading and writing SBML files.

TODO: Add support for sbml-fbc package.

@author: Daniel Machado

   Copyright 2013 Novo Nordisk Foundation Center for Biosustainability,
   Technical University of Denmark.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
   
"""
from ..core.model import Model, Metabolite, Reaction, Compartment
from ..core.odemodel import ODEModel
from ..core.cbmodel import CBModel, Gene
from ..core.fixes import fix_bigg_model

from collections import OrderedDict
from libsbml import SBMLReader, SBMLWriter, SBMLDocument, XMLNode, AssignmentRule, parseL3FormulaWithModel

CB_MODEL = 'cb'
ODE_MODEL = 'ode'

BIGG_MODEL = 'bigg'

LB_TAG = 'LOWER_BOUND'
UB_TAG = 'UPPER_BOUND'
OBJ_TAG = 'OBJECTIVE_COEFFICIENT'
GPR_TAG = 'GENE_ASSOCIATION:'

DEFAULT_SBML_LEVEL = 3
DEFAULT_SBML_VERSION = 1


def load_sbml_model(filename, kind=None):
    """ Loads a metabolic model from a file.
    
    Arguments:
        filename : String -- SBML file path
        kind : {None (default), CB_MODEL, ODE_MODEL} -- define kind of model to load (optional)
    
    Returns:
        Model -- Simple model or respective subclass
    """
    reader = SBMLReader()
    document = reader.readSBML(filename)
    sbml_model = document.getModel()

    if sbml_model is None:
        raise IOError('Failed to load model.')

    if kind and kind.lower() == CB_MODEL:
        model = _load_cbmodel(sbml_model)
    elif kind and kind.lower() == ODE_MODEL:
        model = _load_odemodel(sbml_model)
    else:
        model = _load_stoichiometric_model(sbml_model)

    return model


def load_cbmodel(filename, flavor=None):
    model = load_sbml_model(filename, kind=CB_MODEL)

    if flavor and flavor.lower() == BIGG_MODEL:
        fix_bigg_model(model)

    return model


def load_odemodel(filename):
    return load_sbml_model(filename, ODE_MODEL)


def _load_stoichiometric_model(sbml_model):
    model = Model(sbml_model.getId())
    model.add_compartments(_load_compartments(sbml_model))
    model.add_metabolites(_load_metabolites(sbml_model))
    model.add_reactions(_load_reactions(sbml_model))
    return model


def _load_compartments(sbml_model):
    return [_load_compartment(compartment) for compartment in sbml_model.getListOfCompartments()]


def _load_compartment(compartment):
    return Compartment(compartment.getId(), compartment.getName(), compartment.getSize())


def _load_metabolites(sbml_model):
    return [_load_metabolite(species) for species in sbml_model.getListOfSpecies()]


def _load_metabolite(species):
    return Metabolite(species.getId(), species.getName(), species.getCompartment())


def _load_reactions(sbml_model):
    return [_load_reaction(reaction) for reaction in sbml_model.getListOfReactions()]


def _load_reaction(reaction):

    stoichiometry = OrderedDict()

    for reactant in reaction.getListOfReactants():
        m_id = reactant.getSpecies()
        coeff = -reactant.getStoichiometry()
        if m_id not in stoichiometry:
            stoichiometry[m_id] = coeff
        else:
            stoichiometry[m_id] += coeff

    for product in reaction.getListOfProducts():
        m_id = product.getSpecies()
        coeff = product.getStoichiometry()
        if m_id not in stoichiometry:
            stoichiometry[m_id] = coeff
        else:
            stoichiometry[m_id] += coeff
        if stoichiometry[m_id] == 0.0:
            del stoichiometry[m_id]

    modifiers = [modifier.getSpecies() for modifier in reaction.getListOfModifiers()]

    return Reaction(reaction.getId(), reaction.getName(), reaction.getReversible(), stoichiometry, modifiers)


def _load_cbmodel(sbml_model):
    model = CBModel(sbml_model.getId())
    model.add_compartments(_load_compartments(sbml_model))
    model.add_metabolites(_load_metabolites(sbml_model))
    model.add_reactions(_load_reactions(sbml_model))
    model.set_multiple_bounds(_load_bounds(sbml_model))
    model.set_objective(_load_objective_coefficients(sbml_model))
    genes, rules, reaction_genes = _load_gpr(sbml_model)
    model.add_genes(genes)
    model.set_rules(rules)
    for r_id, gene_set in reaction_genes.items():
        model.reaction_genes[r_id] = gene_set
    return model


def _load_bounds(sbml_model):
    return [(reaction.getId(), (_get_cb_parameter(reaction, LB_TAG),
                                _get_cb_parameter(reaction, UB_TAG)))
            for reaction in sbml_model.getListOfReactions()]


def _load_objective_coefficients(sbml_model):
    return [(reaction.getId(), _get_cb_parameter(reaction, OBJ_TAG, default_value=0))
            for reaction in sbml_model.getListOfReactions()]


def _get_cb_parameter(reaction, tag, default_value=None):
    param_value = default_value
    kinetic_law = reaction.getKineticLaw()
    if kinetic_law:
        parameter = kinetic_law.getParameter(tag)
        if parameter:
            param_value = parameter.getValue()
    return param_value


def _load_gpr(sbml_model):
    genes = set()
    rules = []
    reaction_genes = {}
    for reaction in sbml_model.getListOfReactions():
        rule = _extract_rule(reaction)
        new_genes = rule.replace('(', '').replace(')', '').replace(' and ', ' ').replace(' or ', ' ').split()
        genes = genes | set(new_genes)
        rules.append((reaction.getId(), rule))
        reaction_genes[reaction.getId()] = set(new_genes)
    genes = [Gene(gene) for gene in sorted(genes)]
    return genes, rules, reaction_genes


def _extract_rule(reaction):
    notes = reaction.getNotesString()
    if GPR_TAG in notes:
        rule = notes.partition(GPR_TAG)[2].partition('<')[0].strip()
    else:
        rule = ''
    return rule


def _load_odemodel(sbml_model):
    model = ODEModel(sbml_model.getId())
    model.add_compartments(_load_compartments(sbml_model))
    model.add_metabolites(_load_metabolites(sbml_model))
    model.add_reactions(_load_reactions(sbml_model))
    model.set_concentrations(_load_concentrations(sbml_model))
    model.set_constant_parameters(_load_constant_parameters(sbml_model))
    model.set_variable_parameters(_load_variable_parameters(sbml_model))
    model.set_local_parameters(_load_local_parameters(sbml_model))
    model.set_ratelaws(_load_ratelaws(sbml_model))
    model.set_assignment_rules(_load_assignment_rules(sbml_model))

    return model


def _load_concentrations(sbml_model):
    return [(species.getId(), species.getInitialConcentration())
            for species in sbml_model.getListOfSpecies()]


def _load_constant_parameters(sbml_model):
    return [(parameter.getId(), parameter.getValue())
            for parameter in sbml_model.getListOfParameters() if parameter.getConstant()]

def _load_variable_parameters(sbml_model):
    return [(parameter.getId(), parameter.getValue())
            for parameter in sbml_model.getListOfParameters() if parameter.getConstant() == False]

def _load_local_parameters(sbml_model):
    params = OrderedDict()
    for reaction in sbml_model.getListOfReactions():
        params[reaction.getId()] = [(parameter.getId(), parameter.getValue())
                                    for parameter in reaction.getKineticLaw().getListOfParameters()]
    return params


def _load_ratelaws(sbml_model):
    return [(reaction.getId(), reaction.getKineticLaw().getFormula())
            for reaction in sbml_model.getListOfReactions()]

def _load_assignment_rules(sbml_model):
    return [(rule.getVariable(), rule.getFormula()) for rule in sbml_model.getListOfRules()
            if isinstance(rule, AssignmentRule)]


def save_sbml_model(model, filename):
    """ Save a model to an SBML file.
    
    Arguments:
        model : Model (or any subclass) -- model
        filename : String -- SBML file path
    """

    document = SBMLDocument(DEFAULT_SBML_LEVEL, DEFAULT_SBML_VERSION)
    sbml_model = document.createModel(model.id)
    _save_compartments(model, sbml_model)
    _save_metabolites(model, sbml_model)
    _save_reactions(model, sbml_model)
    if isinstance(model, CBModel):
        _save_cb_parameters(model, sbml_model)
        _save_gpr(model, sbml_model)
    if isinstance(model, ODEModel):
        _save_concentrations(model, sbml_model)
        _save_global_parameters(model, sbml_model)
        _save_kineticlaws(model, sbml_model)
        _save_assignment_rules(model, sbml_model)
    writer = SBMLWriter()
    writer.writeSBML(document, filename)


def save_cbmodel(model, filename, flavor=None):

    if flavor and flavor.lower() == BIGG_MODEL:
        old_bounds = model.bounds.copy()
        for r_id, (lb, ub) in model.bounds.items():
            lb = -1000 if lb is None else lb
            ub = 1000 if ub is None else ub
            model.set_flux_bounds(r_id, lb, ub)
        save_sbml_model(model, filename)
        model.bounds = old_bounds
    else:
         save_sbml_model(model, filename)       


def _save_compartments(model, sbml_model):
    for compartment in model.compartments.values():
        sbml_compartment = sbml_model.createCompartment()
        sbml_compartment.setId(compartment.id)
        sbml_compartment.setName(compartment.name)
        sbml_compartment.setSize(compartment.size)


def _save_metabolites(model, sbml_model):
    for metabolite in model.metabolites.values():
        species = sbml_model.createSpecies()
        species.setId(metabolite.id)
        species.setName(metabolite.name)
        species.setCompartment(metabolite.compartment)


def _save_reactions(model, sbml_model):
    for reaction in model.reactions.values():
        sbml_reaction = sbml_model.createReaction()
        sbml_reaction.setId(reaction.id)
        sbml_reaction.setName(reaction.name)
        sbml_reaction.setReversible(reaction.reversible)
        for m_id, coeff in reaction.stoichiometry.items():
            if coeff < 0:
                speciesReference = sbml_reaction.createReactant()
                speciesReference.setSpecies(m_id)
                speciesReference.setStoichiometry(-coeff)
            elif coeff > 0:
                speciesReference = sbml_reaction.createProduct()
                speciesReference.setSpecies(m_id)
                speciesReference.setStoichiometry(coeff)
        for m_id in reaction.modifiers:
            speciesReference = sbml_reaction.createModifier()
            speciesReference.setSpecies(m_id)


def _save_cb_parameters(model, sbml_model):
    for r_id in model.reactions:
        lb, ub = model.bounds[r_id]
        coeff = model.objective[r_id]
        sbml_reaction = sbml_model.getReaction(r_id)
        kineticLaw = sbml_reaction.createKineticLaw()
        kineticLaw.setFormula('0')
        if lb is not None:
            lbParameter = kineticLaw.createParameter()
            lbParameter.setId(LB_TAG)
            lbParameter.setValue(lb)
        if ub is not None:
            ubParameter = kineticLaw.createParameter()
            ubParameter.setId(UB_TAG)
            ubParameter.setValue(ub)
        objParameter = kineticLaw.createParameter()
        objParameter.setId(OBJ_TAG)
        objParameter.setValue(coeff)

                
def _save_gpr(model, sbml_model):
    for r_id in model.reactions:
        sbml_reaction = sbml_model.getReaction(r_id)
        #sbml_reaction.appendNotes(GPR_TAG + ' ' + model.rules[r_id])
        note = XMLNode.convertStringToXMLNode('<html><p>' + GPR_TAG + ' ' + model.rules[r_id] + '</p></html>')
        note.getNamespaces().add('http://www.w3.org/1999/xhtml')
        sbml_reaction.setNotes(note)


def _save_concentrations(model, sbml_model):
    for m_id, value in model.concentrations.items():
        species = sbml_model.getSpecies(m_id)
        species.setInitialConcentration(value)

def _save_global_parameters(model, sbml_model):
    for p_id, value in model.global_parameters.items():
        parameter = sbml_model.createParameter()
        parameter.setId(p_id)
        parameter.setValue(value)

def _save_kineticlaws(model, sbml_model):
    for r_id, ratelaw in model.ratelaws.items():
        sbml_reaction = sbml_model.getReaction(r_id)
        kineticLaw = sbml_reaction.createKineticLaw()
        #kineticLaw.setFormula(ratelaw)
        kineticLaw.setMath(parseL3FormulaWithModel(ratelaw, sbml_model)) #avoids conversion of Pi to pi
        for p_id, value in model.local_parameters[r_id].items():
            parameter = kineticLaw.createParameter()
            parameter.setId(p_id)
            parameter.setValue(value)

def _save_assignment_rules(model, sbml_model):
    for p_id, formula in model.assignment_rules.items():
        rule = sbml_model.createAssignmentRule()
        rule.setVariable(p_id)
        rule.setFormula(formula)
        sbml_model.getParameter(p_id).setConstant(False)
