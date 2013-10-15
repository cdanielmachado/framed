__author__ = 'kaizhuang'

from ..core.models import Metabolite

def add_ratio_constraint(model, r_id_num, r_id_den, ratio):
    """ Add a flux ratio constraint to a model.

    Arguments:
        model : StoichiometricModel
        r_id_num : str -- id of the numerator
        r_id_num : str -- id of the denominator
        ratio : float -- ratio value

    Returns:
        m_id : id of the pseudo-metabolite used to enable ratio constraint
    """

    if r_id_num in model.reactions and r_id_den in model.reactions:
        m_id = 'ratio_{}_{}'.format(r_id_num, r_id_den)
        model.add_metabolite(Metabolite(m_id))
        model.add_stoichiometry([(m_id, r_id_num, 1), (m_id, r_id_den, -ratio)])

    return m_id


def remove_ratio_constraint(model, m_id):
    """ Remove a flux ratio constraint from a model

    Arguments:
        model : StoichiometricModel
        m_id : id of the pseudo-metabolite used to establish the ratio constraint
    """
    model.remove_metabolites(m_id)


def list_pseudo_metabolites(model):
    """ Show all ratio constraints in the model.
        Specifically, this function prints out the ids of the pseudo-metabolites used to enable the ratio constraints.

    Arguments:
        model : StoichiometricModel
    """
    for m_id in model.metabolites.keys():
        if 'ratio_' in m_id:
            print m_id

