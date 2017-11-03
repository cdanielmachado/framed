""" This module implements gene-wise reformulations of some commonly used constraint-based methods.

Notes:

    Essentially, all methods use gene-level objective functions (rather than reaction-level). For each gene (e.g: b0001)
    a new variable (u_b0001) is introduced that represents the amount of flux carried by the respective enzyme.
    Constraints can also be specified at gene level as a dictionary with the format {'b0001': (lb, ub)}.

    For instance, to simulate the enzyme/flux redistribution of E. coli after knockout of pfkA (b3916) but not pfkB
    (i.e. reaction PFK is still active) using gene-wise MOMA as simulation method (pFBA, lMOMA and ROOM are also available):

    ::

        sol = gene_MOMA(model, constraints={'u_b3916': (0,0)})

        print sol.show_values(pattern='R_')  # look at reaction fluxes
        print sol.show_values(pattern='u_')  # look at enzyme usage values


    See http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005140 for details.

Author: Daniel Machado

"""

from ..model.transformation import gpr_transform
from .simulation import pFBA, lMOMA, MOMA, ROOM
from functools import wraps


def gene_wise(method):
    """ Generic decorator to transform reaction-based simulation methods to gene-based methods. """

    @wraps(method)
    def func_wrapper(model, transformed=False, constraints=None, **kwargs):
        if transformed:
            reactions = [r_id for r_id in model.reactions if r_id.startswith('u_')]
        else:
            model, reactions = gpr_transform(model, inplace=False)

        if constraints:
            constraints = model.convert_constraints(constraints)

        if kwargs.get('reference') is not None:
            reference = kwargs['reference']
            if not set(reference.keys()).issubset(model.reactions):
                raise RuntimeError('Reference fluxes must be calculated for extended model.')

        sol = method(model, constraints=constraints, reactions=reactions, **kwargs)

        sol.extended = sol.values
        sol.values = model.convert_fluxes(sol.values)

        return sol
    return func_wrapper


@gene_wise
def gene_pFBA(model, objective=None, minimize=False, constraints=None, reactions=None, solver=None):
    """ Perform a gene-wise version of parsimonious FBA.

    This method minimizes the total enzyme usage by re-formulating the pFBA objective function at gene level.
    It is also possible to introduce constraints directly at gene level as well.

    Notes:
        See http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005140 for details.

    Arguments:
        model (CBModel): a constraint-based model
        transformed (bool): if the model is already expanded with stoichiometric integration of GPRs (default: False)
        constraints (dict): additional constraints (at gene or reaction level) (optional)
        **kwargs: any other arguments supported by pFBA (see documentation for details)

    Returns:
        Solution: solution

    """
    return pFBA(model, objective=objective, minimize=minimize, constraints=constraints, reactions=reactions, solver=solver)


@gene_wise
def gene_MOMA(model, reference=None, constraints=None, reactions=None, solver=None):
    """ Perform a gene-wise version of MOMA.

    This method minimizes the metabolic adjustment at gene level (rather than reaction level).
    It is also possible to introduce constraints directly at gene level as well.

    Notes:
        See http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005140 for details.

    Arguments:
        model (CBModel): a constraint-based model
        transformed (bool): if the model is already expanded with stoichiometric integration of GPRs (default: False)
        constraints (dict): additional constraints (at gene or reaction level) (optional)
        **kwargs: any other arguments supported by pFBA (see documentation for details)

    Returns:
        Solution: solution
    """

    return MOMA(model, reference=reference, constraints=constraints, reactions=reactions, solver=solver)


@gene_wise
def gene_lMOMA(model, reference=None, constraints=None, reactions=None, solver=None):
    """ Perform a gene-wise version of linear MOMA.

    This method minimizes the metabolic adjustment at gene level (rather than reaction level).
    It is also possible to introduce constraints directly at gene level as well.

    Notes:
        See http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005140 for details.

    Arguments:
        model (CBModel): a constraint-based model
        transformed (bool): if the model is already expanded with stoichiometric integration of GPRs (default: False)
        constraints (dict): additional constraints (at gene or reaction level) (optional)
        **kwargs: any other arguments supported by pFBA (see documentation for details)

    Returns:
        Solution: solution
    """

    return lMOMA(model, reference=reference, constraints=constraints, reactions=reactions, solver=solver)


@gene_wise
def gene_ROOM(model, reference=None, constraints=None, reactions=None, solver=None, delta=0.03, epsilon=0.001):
    """ Perform a gene-wise version of ROOM.

    This method minimizes the regulatory on/off adjustment at gene level (rather than reaction level).
    It is also possible to introduce constraints directly at gene level as well.

    Notes:
        See http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005140 for details.

    Arguments:
        model (CBModel): a constraint-based model
        transformed (bool): if the model is already expanded with stoichiometric integration of GPRs (default: False)
        constraints (dict): additional constraints (at gene or reaction level) (optional)
        **kwargs: any other arguments supported by pFBA (see documentation for details)

    Returns:
        Solution: solution
    """

    return ROOM(model, reference=reference, constraints=constraints, reactions=reactions,
                solver=solver, delta=delta, epsilon=epsilon)
