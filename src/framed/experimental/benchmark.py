from .fluxutils import fit_fluxes_to_model, flux_distance
from ..cobra.deletion import gene_deletion


def evaluate(method, model, dataset, condition, **kwargs):

    measured_fluxes = kwargs.get('measured_fluxes', None)
    fit_fluxes = kwargs.get('fit_fluxes', False)
    fit_growth = kwargs.get('fit_growth', False)
    fit_quadratic = kwargs.get('fit_quadratic', False)
    fit_model = kwargs.get('fit_model', model)
    error_quadratic = kwargs.get('error_quadratic', False)
    error_normalize = kwargs.get('error_normalize', True)

    if fit_growth:
        growth = dataset.growth_rate[condition]
        biomass = model.biomass_reaction

    if measured_fluxes and isinstance(measured_fluxes, dict):
        measured_fluxes = measured_fluxes[condition] # allows condition-specific measurements

    if measured_fluxes:
        input_fluxes = dataset.get_fluxomics(measured_fluxes, condition).to_dict()

        if fit_fluxes:
            fit_constraints = {biomass: growth} if fit_growth else None
            input_fluxes = fit_fluxes_to_model(fit_model, input_fluxes, fit_constraints, quadratic=fit_quadratic)
    else:
        input_fluxes = {}

    constraints = {r_id: val for r_id, val in input_fluxes.items()}

    if fit_growth:
        constraints[biomass] = growth

    sol, sim_fluxes = method(model, constraints, dataset, condition, **kwargs)

    exp_fluxes = dataset.get_fluxomics(conditions=condition).to_dict()
    error = flux_distance(exp_fluxes, sim_fluxes, normalize=error_normalize, quadratic=error_quadratic)

    return error, sim_fluxes, sol


def benchmark(method, model, dataset, conditions=None, **kwargs):

    if conditions is None:
        conditions = dataset.conditions

    results = {}
    for condition in conditions:
        results[condition] = evaluate(method, model, dataset, condition, **kwargs)

    return results


def run_method(method, model, constraints, dataset, condition, **kwargs):

    genes = dataset.get_gene_deletions(condition)
    reference = None

    if method in ['MOMA', 'lMOMA', 'ROOM']:
        ref_condition = kwargs.get('reference_condition', None)
        if ref_condition:
            run_pFBA = lambda *args, **kwargs: run_method('pFBA', *args, **kwargs)
            _, _, sol = evaluate(run_pFBA, model, dataset, ref_condition, **kwargs)
            reference = sol.values
        else:
            raise NameError("Must specify reference condition (named argument: 'reference_condition') for MOMA/lMOMA/ROOM simulations.")

    sol = gene_deletion(model, genes, method=method, reference=reference, constraints=constraints, compute_silent_deletions=True)

    return sol, sol.values

