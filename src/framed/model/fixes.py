"""
Fixes to clean up common models from different sources/groups.

Author: Daniel Machado
   
"""


def fix_cb_model(model, flavor=None):
    """ Apply fixes to known issues in a constraint baised model.

    Args:
        model (CBModel): model
        flavor (str): model *flavor* (optional)

    Notes:
        Currently supported *flavors*: 'cobra' (old style format from the Palsson lab)

    """

    if flavor == 'cobra':
        fix_cobra_model(model)
    else:
        default_fixes(model)


def default_fixes(model):
    """ Apply default fixes to a constraint-based model

    Args:
        model: CBModel

    """
    remove_boundary_metabolites(model)
    fix_reversibility(model)
    clean_bounds(model)
    fix_sink_reactions(model)


def fix_cobra_model(model, remove_boundary=True, set_reversibilty=True, use_infinity=True, clean_ids=True, fix_sinks=True):

    if remove_boundary:
        remove_boundary_metabolites(model, tag='_b')
    if set_reversibilty:
        fix_reversibility(model)
    if use_infinity:
        clean_bounds(model)
    if clean_ids:
        clean_bigg_ids(model)
    if fix_sinks:
        fix_sink_reactions(model)


def fix_sink_reactions(model):
    exchange_compartments = {}
    for r in model.reactions.itervalues():
        if not r.is_exchange: continue

        for m_id in r.stoichiometry:
            met = model.metabolites[m_id]
            if met.boundary: continue
            if met.compartment not in exchange_compartments:
                exchange_compartments[met.compartment] = []
            exchange_compartments[met.compartment].append(r.id)

    if exchange_compartments:
        extracellular = max(exchange_compartments.iteritems(), key=lambda x: len(x[1]))[0]
        for compartment, reactions in exchange_compartments.iteritems():
            if compartment == extracellular:
                continue

            for r_id in reactions:
                reaction = model.reactions[r_id]
                reaction.is_sink = True
                reaction.is_exchange = False


def remove_boundary_metabolites(model, tag=None):
    """ Remove remove boundary metabolites. """

    if tag:
        boundary = filter(lambda m_id: m_id.endswith(tag), model.metabolites)
    else:
        boundary = [m_id for m_id, met in model.metabolites.items() if met.boundary]

    model.remove_metabolites(boundary)


def fix_reversibility(model):
    """ Make reaction reversibility consistent with the bounds. """

    for reaction in model.reactions.values():
        reaction.reversible = (reaction.lb is None or reaction.lb < 0)


def clean_bounds(model, threshold=1000):
    """ Remove artificially large bounds (unbounded = no bounds). """

    for reaction in model.reactions.values():
        reaction.lb = reaction.lb if reaction.lb > -threshold else None
        reaction.ub = reaction.ub if reaction.ub < threshold else None


def clean_bigg_ids(model):
    model._clear_temp()

    clean = lambda x: x.replace('_LPAREN_', '_').replace('_RPAREN_', '_').replace('_DASH_', '__').rstrip('_')

    def key_replace(ord_dict, key, new_key):
        item = ord_dict[key]
        del ord_dict[key]
        ord_dict[new_key] = item

    for m_id, metabolite in model.metabolites.items():
        metabolite.id = clean(m_id)
        key_replace(model.metabolites, m_id, metabolite.id)

    for r_id, reaction in model.reactions.items():
        reaction.id = clean(r_id)
        key_replace(model.reactions, r_id, reaction.id)

        for m_id in reaction.stoichiometry.keys():
            key_replace(reaction.stoichiometry, m_id, clean(m_id))
