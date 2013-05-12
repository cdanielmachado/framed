""" 
Package implementing interfaces to common solvers.

"""
from warnings import warn

solvers = dict()

try:
    from .pulp_wrapper import PuLPSolver
    solvers['pulp'] = PuLPSolver
except:
    pass

try:
    from .gurobi_wrapper import GurobiSolver
    solvers['gurobi'] = GurobiSolver
except:
    pass

default_solver = 'gurobi' 

def set_default_solver(solvername):
    if solvername.lower() in solvers.keys():
        default_solver = solvername.lower()
    else:
        warn('Solver ' + solvername + ' not available.')

def solver_instance():
    return solvers[default_solver]()