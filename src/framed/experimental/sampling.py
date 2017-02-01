from collections import OrderedDict
from random import gauss, random

from framed.cobra.simulation import pFBA
from framed.cobra.variability import FVA
from framed.solvers import solver_instance
from framed.solvers.solver import Status


def lp_sampler(model, n_samples=1000, weights=None, constraints=None, select_probability=0.01,
               futile_cycle_threshold=1e2, variation_threshold=1e-4, merge_keys=False, verbose=True):

    if not weights:
        variability = FVA(model, constraints=constraints)

        weights = {r_id: 1.0/(ub - lb) for r_id, (lb, ub) in variability.items()
                   if ub is not None and lb is not None
                       and variation_threshold < (ub - lb) < futile_cycle_threshold}

    samples = []
    solver = solver_instance(model)

    for i in range(n_samples):
        objective = {r_id: gauss(0, 1)*W for r_id, W in weights.items()
                     if random() < select_probability}

        sol = pFBA(model, objective=objective, constraints=constraints, solver=solver)
        if sol.status == Status.OPTIMAL:
            samples.append(sol.values)

    if verbose:
        print 'Sampling success rate: {} (of {})'.format(len(samples), n_samples)

    if merge_keys:
        merged = OrderedDict()
        for r_id in model.reactions:
            merged[r_id] = [sample[r_id] for sample in samples]
        samples = merged

    return samples
