""" This module implements methods for reading and writing SBML files.

Author: Daniel Machado
   
"""

from ..model.model import Model, Metabolite, Reaction, Compartment
from ..model.odemodel import ODEModel
from ..model.cbmodel import CBModel, Gene, Protein, GPRAssociation
from ..model.fixes import fix_cb_model

import os
from collections import OrderedDict
from sympy.parsing.sympy_parser import parse_expr
from sympy import to_dnf, Or, And
from sympy.logic.boolalg import is_dnf
from libsbml import SBMLReader, SBMLWriter, SBMLDocument, XMLNode, AssignmentRule, parseL3FormulaWithModel, FbcExtension

import cgi
import warnings
import re

DEFAULT_SBML_LEVEL = 3
DEFAULT_SBML_VERSION = 1

CB_MODEL = 'cb'
ODE_MODEL = 'ode'

LB_TAG = 'LOWER_BOUND'
UB_TAG = 'UPPER_BOUND'
OBJ_TAG = 'OBJECTIVE_COEFFICIENT'
GPR_TAG = 'GENE_ASSOCIATION'

DEFAULT_LOWER_BOUND_ID = 'cobra_default_lb'
DEFAULT_UPPER_BOUND_ID = 'cobra_default_ub'
DEFAULT_ZERO_BOUND_ID = 'cobra_0_bound'

DEFAULT_LOWER_BOUND = -1000
DEFAULT_UPPER_BOUND = 1000

ACTIVATOR_TAG = 'SBO:0000459'
INHIBITOR_TAG = 'SBO:0000020'

non_alphanum = re.compile('\W+')
re_type = type(non_alphanum)


class Flavor:
    """ Enumeration of available model flavors. """

    COBRA = 'cobra'  # UCSD models in the old cobra toolbox format
    COBRA_OTHER = 'cobra:other'  # other models using the old cobra toolbox format
    SEED = 'seed'  # modelSEED format
    BIGG = 'bigg'  # BiGG database format (uses sbml-fbc2)
    FBC2 = 'fbc2'  # other models in sbml-fbc2 format


def load_sbml_model(filename, kind=None, flavor=None, exchange_detection_mode=None):
    """ Loads a metabolic model from a file.
    
    Arguments:
        filename (str): SBML file path
        kind (str): define kind of model to load ('cb' or 'ode', optional)
        flavor (str): adapt to different modeling conventions (optional, see Notes)
        exchange_detection_mode (str): detect exchange reactions (optional, see Notes)
    
    Returns:
        Model: Simple model or respective subclass

    Notes:
        Currently supported flavors:
            * 'cobra': UCSD models in the old cobra toolbox format
            * 'cobra:other': other models using the old cobra toolbox format
            * 'seed': modelSEED format
            * 'bigg': BiGG database format (uses sbml-fbc2)
            * 'fbc2': other models using sbml-fbc2

        Supported exchange detection modes:
            * 'unbalanced': Exchange reactions is the one that have either only reactants or products
            * 'boundary': Exchange reaction is the one that have single boundary metabolite on one side
            * <regular expression pattern>: Regular expression which is executed against reaction ID

        Note that some flavors (cobra, bigg) have their own exchange detection mode.

    """
    if not os.path.exists(filename):
        raise IOError("Model file was not found")

    reader = SBMLReader()
    document = reader.readSBML(filename)
    sbml_model = document.getModel()

    if sbml_model is None:
        document.printErrors()
        raise IOError('Failed to load model.')

    if kind and kind.lower() == CB_MODEL:
        model = _load_cbmodel(sbml_model, flavor, exchange_detection_mode=exchange_detection_mode)
    elif kind and kind.lower() == ODE_MODEL:
        model = _load_odemodel(sbml_model)
    else:
        model = _load_stoichiometric_model(sbml_model)

    _load_metadata(sbml_model, model)

    return model


def load_cbmodel(filename, flavor=None, exchange_detection_mode=None):
    """
    Args:
        filename (str): SBML file path
        flavor (str): adapt to different modeling conventions (optional, see Notes)
        exchange_detection_mode (str): detect exchange reactions (optional, see Notes)

    Returns:
        Model: Simple model or respective subclass

    Notes:
        Currently supported flavors:
            * 'cobra': UCSD models in the old cobra toolbox format
            * 'cobra:other': other models using the old cobra toolbox format
            * 'seed': modelSEED format
            * 'bigg': BiGG database format (uses sbml-fbc2)
            * 'fbc2': other models using sbml-fbc2

        Supported exchange detection modes:
            * 'unbalanced': Exchange reactions is the one that have either only reactants or products
            * 'boundary': Exchange reaction is the one that have single boundary metabolite on one side
            * <regular expression pattern>: Regular expression which is executed against reaction ID

        Note that some flavors (cobra, bigg) have their own exchange detection mode.
            
    Returns:
        CBModel: constraint-based model
    """

    model = load_sbml_model(filename, kind=CB_MODEL, flavor=flavor, exchange_detection_mode=exchange_detection_mode)

    fix_cb_model(model, flavor=flavor)

    return model


def load_odemodel(filename):
    return load_sbml_model(filename, ODE_MODEL)


def _load_stoichiometric_model(sbml_model):
    model = Model(sbml_model.getId())
    _load_compartments(sbml_model, model)
    _load_metabolites(sbml_model, model)
    _load_reactions(sbml_model, model)
    return model


def _load_compartments(sbml_model, model):
    for compartment in sbml_model.getListOfCompartments():
        model.add_compartment(_load_compartment(compartment))


def _load_compartment(compartment):
    comp = Compartment(compartment.getId(), compartment.getName(), compartment.getSize())
    _load_metadata(compartment, comp)
    return comp


def _load_metabolites(sbml_model, model, flavor=None):
    for species in sbml_model.getListOfSpecies():
        model.add_metabolite(_load_metabolite(species, flavor), clear_tmp=False)


def _load_metabolite(species, flavor=None):
    metabolite = Metabolite(species.getId(), species.getName(), species.getCompartment(),
                            species.getBoundaryCondition(), species.getConstant())

    if flavor in {Flavor.BIGG, Flavor.FBC2}:
        fbc_species = species.getPlugin('fbc')
        if fbc_species.isSetChemicalFormula():
            formula = fbc_species.getChemicalFormula()
            metabolite.metadata['FORMULA'] = formula

        if fbc_species.isSetCharge():
            charge = fbc_species.getCharge()
            metabolite.metadata['CHARGE'] = str(charge)

    _load_metadata(species, metabolite)
    return metabolite


def _load_reactions(sbml_model, model, exchange_detection_mode=None):
    for reaction in sbml_model.getListOfReactions():
        r = _load_reaction(reaction, sbml_model=sbml_model, exchange_detection_mode=exchange_detection_mode)
        model.add_reaction(r, clear_tmp=False)


def _load_reaction(reaction, sbml_model, exchange_detection_mode=None):
    """
    Args:
        reaction: <SBMLReaction> object 
        exchange_detection_mode: Argument describing how to detect exchange reaction (possible values
            'unbalanced' - Exchange reactions is the one that have either only reactants or products
            'boundary' - Exchange reaction is the one that have single boundary metabolite on one side
            Regex object - Regular expression which is executed against reaction ID
            None - All reactions are NOT exchange reactions

    Returns:

    """

    stoichiometry = OrderedDict()
    modifiers = OrderedDict()

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

    for modifier in reaction.getListOfModifiers():
        m_id = modifier.getSpecies()
        kind = '?'
        sboterm = modifier.getSBOTermID()
        if sboterm == ACTIVATOR_TAG:
            kind = '+'
        if sboterm == INHIBITOR_TAG:
            kind = '-'
        modifiers[m_id] = kind

    is_exchange = False
    if exchange_detection_mode == "unbalanced":
        sign = None
        is_exchange = True
        for m_id, c in stoichiometry.iteritems():
            if sign is None:
                sign = c > 0
            else:
                if sign != c > 0:
                    is_exchange = False
    elif exchange_detection_mode == "boundary":
        products = {m_id for m_id, c in stoichiometry.iteritems() if c > 0}
        reactants = {m_id for m_id, c in stoichiometry.iteritems() if c < 0}
        boundary_products = {m_id for m_id in products if sbml_model.getSpecies(m_id).getBoundaryCondition()}
        is_exchange = (boundary_products and not (products - boundary_products))
        if not is_exchange:
            boundary_reactants = {m_id for m_id in products if sbml_model.getSpecies(m_id).getBoundaryCondition()}
            is_exchange = (boundary_reactants and not (reactants - boundary_reactants))
    elif exchange_detection_mode is None:
        pass
    else:
        is_exchange = exchange_detection_mode.match(reaction.getId()) is not None

    rxn = Reaction(reaction.getId(), name=reaction.getName(), reversible=reaction.getReversible(),
                   stoichiometry=stoichiometry, regulators=modifiers, is_exchange=is_exchange)
    _load_metadata(reaction, rxn)
    return rxn


def _load_cbmodel(sbml_model, flavor, exchange_detection_mode=None):
    if exchange_detection_mode and exchange_detection_mode not in {None, 'unbalanced', 'boundary'}:
        try:
            exchange_detection_mode = re.compile(exchange_detection_mode)
        except:
            raise RuntimeError("Exchange detection mode must be: 'unbalanced', 'boundary', or a valid regular expression.")

    if exchange_detection_mode is None:
        if flavor in {Flavor.COBRA, Flavor.BIGG}:
            exchange_detection_mode = re.compile('^R_EX')
        elif flavor in {Flavor.COBRA_OTHER, Flavor.SEED}:
            exchange_detection_mode = 'boundary'
        elif flavor in {Flavor.FBC2}:
            exchange_detection_mode = 'unbalanced'

    model = CBModel(sbml_model.getId())
    _load_compartments(sbml_model, model)
    _load_metabolites(sbml_model, model, flavor)
    _load_reactions(sbml_model, model, exchange_detection_mode=exchange_detection_mode)
    if flavor in {None, Flavor.COBRA, Flavor.COBRA_OTHER, Flavor.SEED}:
        _load_cobra_bounds(sbml_model, model)
        _load_cobra_objective(sbml_model, model)
        _load_cobra_gpr(sbml_model, model)
    elif flavor in {Flavor.BIGG, Flavor.FBC2}:
        _load_fbc2_bounds(sbml_model, model)
        _load_fbc2_objective(sbml_model, model)
        _load_fbc2_gpr(sbml_model, model)
    else:
        raise TypeError("Unsupported SBML flavor: {}".format(flavor))

    if exchange_detection_mode and len(model.get_exchange_reactions()) == 0:
        warnings.warn("Exchange reactions were not detected")

    return model


def _load_cobra_bounds(sbml_model, model):
    for reaction in sbml_model.getListOfReactions():
        default_lb = None if reaction.getReversible() else 0
        lb = _get_cb_parameter(reaction, LB_TAG, default_lb)
        ub = _get_cb_parameter(reaction, UB_TAG)
        model.set_flux_bounds(reaction.getId(), lb, ub)


def _load_cobra_objective(sbml_model, model):
    objective = OrderedDict()
    for reaction in sbml_model.getListOfReactions():
        coeff = _get_cb_parameter(reaction, OBJ_TAG, default_value=0)
        if coeff:
            objective[reaction.getId()] = coeff
    model.set_objective(objective)


def _get_cb_parameter(reaction, tag, default_value=None):
    param_value = default_value
    kinetic_law = reaction.getKineticLaw()
    if kinetic_law:
        parameter = kinetic_law.getParameter(tag)
        if parameter:
            param_value = parameter.getValue()
    return param_value


def _load_cobra_gpr(sbml_model, model):
    genes = set()
    gprs = OrderedDict()

    for reaction in sbml_model.getListOfReactions():
        rule = model.reactions[reaction.getId()].metadata.pop(GPR_TAG, None)
        if rule:
            gpr = parse_gpr_rule(rule, prefix='G_')
            for protein in gpr.proteins:
                genes |= set(protein.genes)
            gprs[reaction.getId()] = gpr
        else:
            gprs[reaction.getId()] = None

    for gene in sorted(genes):
        model.add_gene(Gene(gene, gene[2:]))

    for r_id, gpr in gprs.items():
        model.set_gpr_association(r_id, gpr, add_genes=False)


def sanitize_id(identifier):
    return non_alphanum.sub('_', identifier)


def parse_gpr_rule(rule, prefix=None):

    if not rule:
        return None

    rule = rule.replace('(', '( ').replace(')', ' )')

    def replacement(token):
        if token.lower() == 'and':
            return '&'
        elif token.lower() == 'or':
            return '|'
        elif token == '(' or token == ')':
            return token
        elif prefix is not None and not token.startswith(prefix):
            return prefix + sanitize_id(token)
        else:
            return sanitize_id(token)

    rule = ' '.join(map(replacement, rule.split()))

    expr = parse_expr(rule)

    if not is_dnf(expr):
        expr = to_dnf(expr)

    gpr = GPRAssociation()

    if type(expr) is Or:
        for sub_expr in expr.args:
            protein = Protein()
            if type(sub_expr) is And:
                protein.genes = [str(gene) for gene in sub_expr.args]
            else:
                protein.genes = [str(sub_expr)]
            gpr.proteins.append(protein)
    elif type(expr) is And:
        protein = Protein()
        protein.genes = [str(gene) for gene in expr.args]
        gpr.proteins = [protein]
    else:
        protein = Protein()
        protein.genes = [str(expr)]
        gpr.proteins = [protein]


    return gpr


def _load_fbc2_bounds(sbml_model, model):
    params = {param.getId(): param.getValue() for param in sbml_model.getListOfParameters()}

    for reaction in sbml_model.getListOfReactions():
        fbc_rxn = reaction.getPlugin('fbc')
        lb = fbc_rxn.getLowerFluxBound()
        ub = fbc_rxn.getUpperFluxBound()
        model.set_flux_bounds(reaction.getId(), params[lb], params[ub])


def _load_fbc2_objective(sbml_model, model):
    fbcmodel = sbml_model.getPlugin('fbc')
    active_obj = fbcmodel.getActiveObjective()
    objective = OrderedDict()
    for rxn_obj in active_obj.getListOfFluxObjectives():
        r_id = rxn_obj.getReaction()
        coeff = rxn_obj.getCoefficient()
        if coeff:
            objective[r_id] = coeff
    model.set_objective(objective)


def _load_fbc2_gpr(sbml_model, model):
    fbcmodel = sbml_model.getPlugin('fbc')

    for gene in fbcmodel.getListOfGeneProducts():
        model.add_gene(Gene(gene.getId(), gene.getName()))

    for reaction in sbml_model.getListOfReactions():
        fbcrxn = reaction.getPlugin('fbc')
        gpr_assoc = fbcrxn.getGeneProductAssociation()
        if gpr_assoc:
            gpr = _parse_fbc_association(gpr_assoc.getAssociation(), reaction.id)
            model.set_gpr_association(reaction.getId(), gpr, add_genes=False)
        else:
            model.set_gpr_association(reaction.getId(), None)


def _parse_fbc_association(gpr_assoc, reaction_id):

    gpr = GPRAssociation()

    parsing_error = False
    if gpr_assoc.isFbcOr():
        for item in gpr_assoc.getListOfAssociations():
            protein = Protein()
            if item.isFbcAnd():
                for subitem in item.getListOfAssociations():
                    if subitem.isGeneProductRef():
                        protein.genes.append(subitem.getGeneProduct())
                    else:
                        w = "Gene association for reaction '{}' is not DNF".format(reaction_id)
                        warnings.warn(w, SyntaxWarning)
                        parsing_error = True
            elif item.isGeneProductRef:
                protein.genes.append(item.getGeneProduct())
            else:
                w = "Gene association for reaction '{}' is not DNF".format(reaction_id)
                warnings.warn(w, SyntaxWarning)
                parsing_error = True
            gpr.proteins.append(protein)

    elif gpr_assoc.isFbcAnd():
        protein = Protein()
        for item in gpr_assoc.getListOfAssociations():
            if item.isGeneProductRef():
                protein.genes.append(item.getGeneProduct())
            else:
                w = "Gene association for reaction '{}' is not DNF".format(reaction_id)
                warnings.warn(w, SyntaxWarning)
                parsing_error = True
        gpr.proteins = [protein]
    elif gpr_assoc.isGeneProductRef():
        protein = Protein()
        protein.genes = [gpr_assoc.getGeneProduct()]
        gpr.proteins = [protein]
    else:
        w = "Gene association for reaction '{}' is not DNF".format(reaction_id)
        warnings.warn(w, SyntaxWarning)
        parsing_error = True

    if not parsing_error:
        return gpr


def _load_odemodel(sbml_model):
    model = ODEModel(sbml_model.getId())
    _load_compartments(sbml_model, model)
    _load_metabolites(sbml_model, model)
    _load_reactions(sbml_model, model)
    _load_concentrations(sbml_model, model)
    _load_global_parameters(sbml_model, model)
    _load_local_parameters(sbml_model, model)
    _load_ratelaws(sbml_model, model)
    _load_assignment_rules(sbml_model, model)

    return model


def _load_concentrations(sbml_model, model):
    for species in sbml_model.getListOfSpecies():
        model.set_concentration(species.getId(), species.getInitialConcentration())


def _load_global_parameters(sbml_model, model):
    for parameter in sbml_model.getListOfParameters():
            model.set_global_parameter(parameter.getId(), parameter.getValue(), parameter.getConstant())


def _load_local_parameters(sbml_model, model):
    for reaction in sbml_model.getListOfReactions():
        for parameter in reaction.getKineticLaw().getListOfParameters():
            model.set_local_parameter(reaction.getId(), parameter.getId(), parameter.getValue())


def _load_ratelaws(sbml_model, model):
    for reaction in sbml_model.getListOfReactions():
        model.set_ratelaw(reaction.getId(), reaction.getKineticLaw().getFormula())


def _load_assignment_rules(sbml_model, model):
    for rule in sbml_model.getListOfRules():
        if isinstance(rule, AssignmentRule):
            model.set_assignment_rule(rule.getVariable(), rule.getFormula())


def save_sbml_model(model, filename, flavor=None):
    """ Save a model to an SBML file.
    
    Arguments:
        model (Model): model
        filename (str): file path
        flavor (str): adapt to different modeling conventions (optional, currently available: 'cobra', 'fbc2')
    """

    document = SBMLDocument(DEFAULT_SBML_LEVEL, DEFAULT_SBML_VERSION)
    sbml_model = document.createModel(model.id)

    if flavor in {Flavor.BIGG, Flavor.FBC2}:
        document.enablePackage(FbcExtension.getXmlnsL3V1V2(), 'fbc', True)
        fbc_model = sbml_model.getPlugin('fbc')
        fbc_model.setStrict(True)
        document.setPackageRequired('fbc', False)
    _save_compartments(model, sbml_model)
    _save_metabolites(model, sbml_model, flavor)
    _save_reactions(model, sbml_model)
    if isinstance(model, CBModel):
        _save_cb_parameters(model, sbml_model, flavor)
        _save_gpr_associations(model, sbml_model, flavor)
    if isinstance(model, ODEModel):
        _save_concentrations(model, sbml_model)
        _save_global_parameters(model, sbml_model)
        _save_kineticlaws(model, sbml_model)
        _save_assignment_rules(model, sbml_model)
    _save_metadata(model, sbml_model)
    writer = SBMLWriter()
    writer.writeSBML(document, filename)


def save_cbmodel(model, filename, flavor=Flavor.COBRA):
    save_sbml_model(model, filename, flavor)


def _save_compartments(model, sbml_model):
    for compartment in model.compartments.values():
        sbml_compartment = sbml_model.createCompartment()
        sbml_compartment.setId(compartment.id)
        sbml_compartment.setName(compartment.name)
        sbml_compartment.setSize(compartment.size)
        sbml_compartment.setConstant(True)
        _save_metadata(compartment, sbml_compartment)


def _save_metabolites(model, sbml_model, flavor):
    for metabolite in model.metabolites.values():
        species = sbml_model.createSpecies()
        species.setId(metabolite.id)
        species.setName(metabolite.name)
        species.setCompartment(metabolite.compartment)
        species.setBoundaryCondition(metabolite.boundary)
        species.setConstant(metabolite.constant)
        species.setHasOnlySubstanceUnits(True)

        # if flavor in {Flavor.BIGG, Flavor.FBC2}:
        #     fbc_species = species.getPlugin('fbc')
        #
        #     if 'FORMULA' in metabolite.metadata:
        #         fbc_species.setChemicalFormula(metabolite.metadata['FORMULA'])
        #     if 'CHARGE' in metabolite.metadata:
        #         try:
        #             charge = int(metabolite.metadata['CHARGE'])
        #             fbc_species.setCharge(charge)
        #         except ValueError:
        #             pass

        _save_metadata(metabolite, species)


def _save_reactions(model, sbml_model):
    for reaction in model.reactions.values():
        sbml_reaction = sbml_model.createReaction()
        sbml_reaction.setId(reaction.id)
        sbml_reaction.setName(reaction.name)
        sbml_reaction.setReversible(reaction.reversible)
        sbml_reaction.setFast(False)
        _save_metadata(reaction, sbml_reaction)

        for m_id, coeff in reaction.stoichiometry.items():
            if coeff < 0:
                speciesReference = sbml_reaction.createReactant()
                speciesReference.setSpecies(m_id)
                speciesReference.setStoichiometry(-coeff)
                speciesReference.setConstant(True)
            elif coeff > 0:
                speciesReference = sbml_reaction.createProduct()
                speciesReference.setSpecies(m_id)
                speciesReference.setStoichiometry(coeff)
                speciesReference.setConstant(True)
        for m_id, kind in reaction.regulators.items():
            speciesReference = sbml_reaction.createModifier()
            speciesReference.setSpecies(m_id)
            if kind == '+':
                speciesReference.setSBOTerm(ACTIVATOR_TAG)
            if kind == '-':
                speciesReference.setSBOTerm(INHIBITOR_TAG)


def _save_cb_parameters(model, sbml_model, flavor):

    if flavor == Flavor.COBRA:
        _save_cobra_parameters(model, sbml_model, set_default_bounds=True)
    elif flavor in {Flavor.BIGG, Flavor.FBC2}:
        _save_fbc_fluxbounds(model, sbml_model)
        _save_fbc_objective(model, sbml_model)
    else:
        _save_cobra_parameters(model, sbml_model)


def _save_gpr_associations(model, sbml_model, flavor):
    if flavor in {Flavor.BIGG, Flavor.FBC2}:
        _save_fbc_gprs(model, sbml_model)
    else:
        _save_cobra_gprs(model, sbml_model)


def _save_cobra_parameters(model, sbml_model, set_default_bounds=False):
    for r_id, reaction in model.reactions.items():
        sbml_reaction = sbml_model.getReaction(r_id)
        kineticLaw = sbml_reaction.createKineticLaw()
        kineticLaw.setFormula('0')
        lb, ub = reaction.lb, reaction.ub
        if set_default_bounds:
            lb = DEFAULT_LOWER_BOUND if lb is None else lb
            ub = DEFAULT_UPPER_BOUND if ub is None else ub
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
        objParameter.setValue(reaction.objective)


def _save_cobra_gprs(model, sbml_model):
    for r_id, reaction in model.reactions.items():
        if reaction.gpr:
            reaction.metadata[GPR_TAG] = str(reaction.gpr)
            sbml_reaction = sbml_model.getReaction(r_id)
            _save_metadata(reaction, sbml_reaction)


def _save_fbc_fluxbounds(model, sbml_model):

    default_lb = sbml_model.createParameter()
    default_lb.setId(DEFAULT_LOWER_BOUND_ID)
    default_lb.setValue(DEFAULT_LOWER_BOUND)
    default_lb.setConstant(True)

    default_ub = sbml_model.createParameter()
    default_ub.setId(DEFAULT_UPPER_BOUND_ID)
    default_ub.setValue(DEFAULT_UPPER_BOUND)
    default_ub.setConstant(True)

    zero_bound = sbml_model.createParameter()
    zero_bound.setId(DEFAULT_ZERO_BOUND_ID)
    zero_bound.setValue(0)
    zero_bound.setConstant(True)

    for r_id, reaction in model.reactions.items():
        fbcrxn = sbml_model.getReaction(r_id).getPlugin('fbc')

        if reaction.lb is None or reaction.lb <= DEFAULT_LOWER_BOUND:
            fbcrxn.setLowerFluxBound(DEFAULT_LOWER_BOUND_ID)
        elif reaction.lb == 0:
            fbcrxn.setLowerFluxBound(DEFAULT_ZERO_BOUND_ID)
        else:
            lb_id = '{}_lower_bound'.format(r_id)
            lb_param = sbml_model.createParameter()
            lb_param.setId(lb_id)
            lb_param.setValue(reaction.lb)
            lb_param.setConstant(True)
            fbcrxn.setLowerFluxBound(lb_id)

        if reaction.ub is None or reaction.ub >= DEFAULT_UPPER_BOUND:
            fbcrxn.setUpperFluxBound(DEFAULT_UPPER_BOUND_ID)
        elif reaction.ub == 0:
            fbcrxn.setUpperFluxBound(DEFAULT_ZERO_BOUND_ID)
        else:
            ub_id = '{}_upper_bound'.format(r_id)
            ub_param = sbml_model.createParameter()
            ub_param.setId(ub_id)
            ub_param.setValue(reaction.ub)
            ub_param.setConstant(True)
            fbcrxn.setUpperFluxBound(ub_id)


def _save_fbc_objective(model, sbml_model):
    fbcmodel = sbml_model.getPlugin('fbc')
    obj = fbcmodel.createObjective()
    obj.setId('objective')
    fbcmodel.setActiveObjectiveId('objective')
    obj.setType('maximize')
    for r_id, reaction in model.reactions.items():
        if reaction.objective:
            r_obj = obj.createFluxObjective()
            r_obj.setReaction(r_id)
            r_obj.setCoefficient(reaction.objective)


def _save_fbc_gprs(model, sbml_model):
    fbcmodel = sbml_model.getPlugin('fbc')
    for gene in model.genes.values():
        gene_prod = fbcmodel.createGeneProduct()
        gene_prod.setId(gene.id)
        gene_prod.setName(gene.name)
        gene_prod.setLabel(gene.name)

    for r_id, reaction in model.reactions.items():
        if reaction.gpr:
            fbcrxn = sbml_model.getReaction(r_id).getPlugin('fbc')
            gpr_assoc = fbcrxn.createGeneProductAssociation()

            if len(reaction.gpr.proteins) > 1:
                gpr_assoc = gpr_assoc.createOr()

            for protein in reaction.gpr.proteins:
                if len(protein.genes) > 1:
                    protein_assoc = gpr_assoc.createAnd()
                else:
                    protein_assoc = gpr_assoc

                for gene in protein.genes:
                    gene_ref = protein_assoc.createGeneProductRef()
                    gene_ref.setGeneProduct(gene)


def _save_concentrations(model, sbml_model):
    for m_id, value in model.concentrations.items():
        species = sbml_model.getSpecies(m_id)
        species.setInitialConcentration(value)


def _save_global_parameters(model, sbml_model):
    for p_id, value in model.constant_params.items():
        parameter = sbml_model.createParameter()
        parameter.setId(p_id)
        parameter.setValue(value)
        parameter.setConstant(True)
    for p_id, value in model.variable_params.items():
        parameter = sbml_model.createParameter()
        parameter.setId(p_id)
        parameter.setValue(value)
        parameter.setConstant(False)


def _save_kineticlaws(model, sbml_model):
    for r_id, ratelaw in model.ratelaws.items():
        sbml_reaction = sbml_model.getReaction(r_id)
        kineticLaw = sbml_reaction.createKineticLaw()
        #kineticLaw.setFormula(ratelaw)
        kineticLaw.setMath(parseL3FormulaWithModel(ratelaw, sbml_model)) #avoids conversion of Pi to pi
        for p_id, value in model.local_params[r_id].items():
            parameter = kineticLaw.createParameter()
            parameter.setId(p_id)
            parameter.setValue(value)


def _save_assignment_rules(model, sbml_model):
    for p_id, formula in model.assignment_rules.items():
        rule = sbml_model.createAssignmentRule()
        rule.setVariable(p_id)
        rule.setFormula(formula)
        sbml_model.getParameter(p_id).setConstant(False)


def _save_metadata(elem, sbml_elem):
    if elem.metadata:
        try:
            notes = ['<p>{}: {}</p>'.format(key, cgi.escape(value))
                     for key, value in elem.metadata.items()]
            note_string = '<html>' + ''.join(notes) + '</html>'
            note_xml = XMLNode.convertStringToXMLNode(note_string)
            note_xml.getNamespaces().add('http://www.w3.org/1999/xhtml')
            sbml_elem.setNotes(note_xml)
        except AttributeError:
            warnings.warn("Unable to save metadata for object {}:".format(sbml_elem.getId()), RuntimeWarning)


def _load_metadata(sbml_elem, elem):
    notes = sbml_elem.getNotes()

    if notes:
        _recursive_node_parser(notes, elem.metadata)


def _recursive_node_parser(node, cache):
    node_data = node.getCharacters()
    if ':' in node_data:
        key, value = node_data.split(':', 1)
        cache[key.strip()] = value.strip()

    for i in range(node.getNumChildren()):
        _recursive_node_parser(node.getChild(i), cache)