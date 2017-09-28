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


def compare_fluxes(original, other, tolerance=1e-6, abstol=1e-9, sort=False, pattern=None):

    only_left = sorted(set(original.keys()) - set(other.keys()))
    only_right = sorted(set(other.keys()) - set(original.keys()))
    common = sorted(set(original.keys()) & set(other.keys()))

    difference = [(r_id, abs(original[r_id] - other[r_id])) for r_id in common]
    flux_left = [(r_id, original[r_id]) for r_id in only_left]
    flux_right = [(r_id, other[r_id]) for r_id in only_right]

    if pattern is not None:
        difference = filter(lambda (a, b): pattern in a, difference)
        flux_left = filter(lambda (a, b): pattern in a, flux_left)
        flux_right = filter(lambda (a, b): pattern in a, flux_right)

    if sort:
        difference.sort(key=lambda x: x[1], reverse=True)
        flux_left.sort(key=lambda x: x[1], reverse=True)
        flux_right.sort(key=lambda x: x[1], reverse=True)

    for r_id, val in difference:
        if val > tolerance:
            x1 = original[r_id] if abs(original[r_id]) > abstol else 0
            x2 = other[r_id] if abs(other[r_id]) > abstol else 0
            print '{: <16} {: < 10.3g} {: < 10.3g}'.format(r_id, x1, x2)

    for r_id, val in flux_left:
        if abs(val) > tolerance:
            x = original[r_id] if abs(original[r_id]) > abstol else 0
            print '{: <16} {: < 10.3g}   --'.format(r_id, x)

    for r_id, val in flux_right:
        if abs(val) > tolerance:
            x = other[r_id] if abs(other[r_id]) > abstol else 0
            print '{: <16}   --       {: < 10.3g}'.format(r_id, x)


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


