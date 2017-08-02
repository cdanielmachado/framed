"""
Package implementing interfaces to linear programming solvers.

Author: Daniel Machado, Marta Matos

"""

solvers = dict()

try:
    from .gurobi_interface import GurobiSolver

    solvers['gurobi'] = GurobiSolver
except ImportError:
    pass


try:
    from .cplex_interface import CplexSolver

    solvers['cplex'] = CplexSolver
except ImportError:
    pass


default_solver = None


def get_default_solver():

    global default_solver

    if default_solver:
        return default_solver

    solver_order = ['gurobi', 'cplex']

    for solver in solver_order:
        if solver in solvers.keys():
            default_solver = solver
            break

    if not default_solver:
        raise EnvironmentError("No solver available")

    return default_solver


def set_default_solver(solvername):
    """ Sets default solver.

    Arguments:
        solvername : str -- solver name (currently available: gurobi, cplex, glpk)
    """

    global default_solver

    if solvername.lower() in solvers.keys():
        default_solver = solvername.lower()
    else:
        raise EnvironmentError("Solver '{}' not available.".format(solvername))


def solver_instance(model=None):
    """ Returns a new instance of the currently selected solver.

    Arguments:
        model : CBModel (optional) -- immediatly instantiate problem with given model

    Returns:
        Solver
    """

    solver = get_default_solver()

    if solver:
        return solvers[solver](model)
