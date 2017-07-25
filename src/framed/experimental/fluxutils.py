from collections import OrderedDict
from framed.cobra.simulation import MOMA, lMOMA
from framed.solvers.solver import Status
from numpy import array, sqrt, sum, abs


def fit_fluxes_to_model(model, fluxes, constraints=None, quadratic=False):

    if quadratic:
        sol = MOMA(model, reference=fluxes, constraints=constraints)
    else:
        sol = lMOMA(model, reference=fluxes, constraints=constraints)

    if sol.status == Status.OPTIMAL:
        fitted = OrderedDict([(r_id, sol.values[r_id]) for r_id in fluxes])
    else:
        fitted = None

    return fitted


def flux_distance(original, other, normalize=False, quadratic=False):
    x = array(original.values())
    y = array([other[r_id] for r_id in original])

    if quadratic:
        dist = sqrt(sum((x-y)**2))
        size = sqrt(sum(x**2))
    else:
        dist = sum(abs(x-y))
        size = sum(abs(x))

    if normalize:
        return dist / size
    else:
        return dist


def compare_fluxes(original, other, tolerance=1e-6, sort_values=False):

    common = sorted(set(original.keys()) & set(other.keys()))
    diff = [(r_id, abs(original[r_id] - other[r_id])) for r_id in common]

    if sort_values:
        diff.sort(key=lambda x: x[1], reverse=True)

    for r_id, val in diff:
        if val > tolerance:
            print '{: <16} {: < 10.3g} {: < 10.3g}'.format(r_id, original[r_id], other[r_id])


def compute_turnover(model, v):
    m_r_table = model.metabolite_reaction_lookup_table()
    t = {m_id: 0.5*sum([abs(coeff * v[r_id]) for r_id, coeff in neighbours.items()])
         for m_id, neighbours in m_r_table.items()}
    return t


def merge_fluxes(model, mapping, v):
    v_new = OrderedDict()

    for r_id in model.reactions:
        if r_id in mapping:
            fwd_id, bwd_id = mapping[r_id]
            v_new[r_id] = v[fwd_id] - v[bwd_id]
        else:
            v_new[r_id] = v[r_id]

    return v_new


