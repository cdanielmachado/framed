"""
Module for model transformation operations.

Author: Daniel Machado
   
"""

from .model import Reaction, Compartment, Metabolite
from .cbmodel import CBModel, CBReaction, GPRAssociation
from ..cobra.variability import blocked_reactions


def simplify(model, reactions=None, clean_compartments=True, inplace=True):
    """ Removes all blocked reactions in a constraint based model
    
    Arguments:
        model (CBModel): model
        reactions (list): List of reactions which will be checked for being blocked (default: None - check all reactions)
        clean_compartments (bool): remove empty compartments (default: True)
        inplace (bool): change model in place (default), otherwise create a copy first
        
    Returns:
        CBModel: simplified model (if not in place)
    """

    if not inplace:
        model = model.copy()

    del_reactions = blocked_reactions(model, reactions=reactions)
    model.remove_reactions(del_reactions)

    del_metabolites = disconnected_metabolites(model)
    model.remove_metabolites(del_metabolites, safe_delete=False)

    del_genes = disconnected_genes(model)
    model.remove_genes(del_genes)

    if clean_compartments:
        empty = empty_compartments(model)
        model.remove_compartments(empty, delete_metabolites=False)

    if not inplace:
        return model


def make_irreversible(model, inplace=True, reactions=None):
    """ Splits all reversible reactions into forward and backward directions.
    
    Arguments:
        model : Model (or CBmodel)
        inplace (bool): change model inplace (default), otherwise create a copy first
        reactions (list) : split only reactions in this list (optional)

    Returns:
        dict: mapping of old reaction ids to splitted reaction ids
    """

    if not inplace:
        model = model.copy()

    if reactions is None:
        reactions = model.reactions.keys()

    mapping = dict()

    for r_id, reaction in model.reactions.items():
        if reaction.reversible and r_id in reactions:
            fwd_id = reaction.id + '_f'
            bwd_id = reaction.id + '_b'
            mapping[r_id] = (fwd_id, bwd_id)
            bwd_stoichiometry = [(m_id, -coeff) for m_id, coeff in reaction.stoichiometry.items()]

            if isinstance(model, CBModel):
                lb, ub = reaction.lb, reaction.ub
                lb_fwd = max(0, lb) if lb is not None else 0
                ub_fwd = max(0, ub) if ub is not None else None
                lb_bwd = max(-ub, 0) if ub is not None else 0
                ub_bwd = max(-lb, 0) if lb is not None else None
                obj = reaction.objective
                obj_fwd = obj if obj >= 0 else 0
                obj_bwd = -obj if obj < 0 else 0
                r_fwd = CBReaction(fwd_id, reaction.name, False, reaction.stoichiometry, reaction.regulators,
                                   lb_fwd, ub_fwd, obj_fwd, reaction.gpr)
                r_bwd = CBReaction(bwd_id, reaction.name, False, bwd_stoichiometry, reaction.regulators,
                                   lb_bwd, ub_bwd, obj_bwd, reaction.gpr)
            else:
                r_fwd = Reaction(fwd_id, reaction.name, False, reaction.stoichiometry, reaction.regulators)
                r_bwd = Reaction(bwd_id, reaction.name, False, bwd_stoichiometry, reaction.regulators)

            model.add_reaction(r_fwd)
            model.add_reaction(r_bwd)
            model.remove_reaction(r_id)

    return mapping


def disconnected_metabolites(model):
    m_r_table = model.metabolite_reaction_lookup()
    return [m_id for m_id, edges in m_r_table.items() if not edges]


def disconnected_genes(model):
    disconnected = set(model.genes)
    for reaction in model.reactions.values():
        disconnected -= set(reaction.get_associated_genes())
    return disconnected


def empty_compartments(model):
    used = [met.compartment for met in model.metabolites.values()]
    empty = set(model.compartments) - set(used)
    return empty


def split_isozymes(model):
    mapping = dict()

    for r_id, reaction in model.reactions.items():

        if reaction.gpr is not None and len(reaction.gpr.proteins) > 1:
            mapping[r_id] = []
            for i, protein in enumerate(reaction.gpr.proteins):
                r_id_new = '{}_iso{}'.format(reaction.id, i + 1)
                mapping[r_id].append(r_id_new)
                gpr_new = GPRAssociation()
                gpr_new.proteins.append(protein)
                reaction_new = CBReaction(r_id_new, reaction.name,
                                          reversible=reaction.reversible,
                                          stoichiometry=reaction.stoichiometry,
                                          regulators=reaction.regulators,
                                          lb=reaction.lb, ub=reaction.ub,
                                          objective=reaction.objective,
                                          gpr_association=gpr_new)
                model.add_reaction(reaction_new)
            model.remove_reaction(r_id)

    return mapping


def genes_to_species(model, gene_prefix='G_', usage_prefix='u_', pseudo_genes=None):

    if pseudo_genes is None:
        pseudo_genes = {'G_s0001', 'G_S0001', 'G_s_0001', 'G_S_0001'}

    new_reactions = []
    compartment = Compartment('genes', 'gene pool')
    model.add_compartment(compartment)

    for gene in model.genes.values():
        if gene.id in pseudo_genes:
            continue
        model.add_metabolite(Metabolite(gene.id, gene.id, 'genes'))
        r_id = usage_prefix + gene.id[len(gene_prefix):]
        reaction = CBReaction(r_id, r_id, False, {gene.id: 1})
        model.add_reaction(reaction)
        new_reactions.append(r_id)

    for r_id, reaction in model.reactions.items():

        if reaction.gpr is not None:
            if len(reaction.gpr.proteins) > 1:
                print 'error: isozymes not split:', r_id
                return
            elif len(reaction.gpr.proteins) == 1:
                for g_id in reaction.gpr.proteins[0].genes:
                    if g_id not in pseudo_genes:
                        reaction.stoichiometry[g_id] = -1

    return new_reactions


def merge_fluxes(fluxes, mapping_rev, mapping_iso, net=True):

    fluxes = fluxes.copy()

    for r_id, r_ids in mapping_iso.items():
        fluxes[r_id] = sum([fluxes[r_id2] for r_id2 in r_ids])
        for r_id2 in r_ids:
            del fluxes[r_id2]

    for r_id, (fwd_id, bwd_id) in mapping_rev.items():
        if net:
            fluxes[r_id] = fluxes[fwd_id] - fluxes[bwd_id]
        else:
            fluxes[r_id] = fluxes[fwd_id] + fluxes[bwd_id]
        del fluxes[fwd_id]
        del fluxes[bwd_id]

    return fluxes


def convert_constraints(constraints, mapping_rev, mapping_iso):
    constraints = constraints.copy()

    for r_id, (fwd_id, bwd_id) in mapping_rev.items():
        if r_id in constraints:
            x = constraints[r_id]
            lb, ub = x if isinstance(x, tuple) else (x, x)
            lb_fwd = max(0, lb) if lb is not None else 0
            ub_fwd = max(0, ub) if ub is not None else None
            lb_bwd = max(-ub, 0) if ub is not None else 0
            ub_bwd = max(-lb, 0) if lb is not None else None
            constraints[fwd_id] = (lb_fwd, ub_fwd)
            constraints[bwd_id] = (lb_bwd, ub_bwd)
            del constraints[r_id]

    for r_id, r_ids in mapping_iso.items():
        if r_id in constraints:
            x = constraints[r_id]
            ub = x[1] if isinstance(x, tuple) else x
            for r_id2 in r_ids:
                constraints[r_id2] = (0, ub)
            del constraints[r_id]

    return constraints


def gpr_transform(model, inplace=True, gene_prefix='G_', usage_prefix='u_', pseudo_genes=None):
    """ Transformation method that integrates GPR associations directly into the stoichiometric matrix.

    Notes:
        This method extends the stoichiometric matrix where genes become pseudo-metabolites (rows)
        and enzyme usage variables become pseudo-reactions columns.

        See http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005140 for details.

    Args:
        model (CBModel): original model
        inplace (bool): change model in place (default), otherwise create a copy first
        gene_prefix (str): prefix used for gene ids (default: 'G_')
        usage_prefix (str): prefix to create enzyme usage variables (default: 'u_')
        pseudo_genes (list): ignore pseudo/fake genes in the model (default: 's0001')

    Returns:
        CBModel: only when inplace=False
        list: list of newly created enzyme usage variables
    """

    if not inplace:
        model = model.copy()

    mapping_rev = make_irreversible(model)
    mapping_iso = split_isozymes(model)
    u_reactions = genes_to_species(model, gene_prefix=gene_prefix, usage_prefix=usage_prefix, pseudo_genes=pseudo_genes)
    model.convert_fluxes = lambda x, net=True: merge_fluxes(x, mapping_rev, mapping_iso, net)
    model.convert_constraints = lambda x: convert_constraints(x, mapping_rev, mapping_iso)
    if inplace:
        return u_reactions
    else:
        return model, u_reactions

