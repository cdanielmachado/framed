""" 
Package implementing interfaces to common solvers.

@author: Daniel Machado

   Copyright 2013 Novo Nordisk Foundation Center for Biosustainability,
   Technical University of Denmark.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

"""

from warnings import warn
from .solver import Parameter

solvers = dict()

try:
    from .glpk_wrapper import GlpkSolver

    solvers['glpk'] = GlpkSolver
except:
    pass

try:
    from .glpk_wrapper_lazy import GlpkSolverLazy

    solvers['glpk_lazy'] = GlpkSolverLazy
except:
    pass


try:
    from .gurobi_wrapper import GurobiSolver

    solvers['gurobi'] = GurobiSolver
except:
    pass


try:
    from .cplex_wrapper import CplexSolver

    solvers['cplex'] = CplexSolver
except:
    pass

default_solver = 'gurobi'


def set_default_solver(solvername):
    """ Sets default solver.

    Arguments:
        solvername : str -- solver name (currently available: gurobi, pulp)
    """

    global default_solver

    if solvername.lower() in solvers.keys():
        default_solver = solvername.lower()
    else:
        warn('Solver ' + solvername + ' not available.')


def solver_instance(model=None):
    """ Returns a new instance of the currently selected solver.

    Arguments:
        model : CBModel (optional) -- immediatly instantiate problem with given model

    Returns:
        Solver
    """

    return solvers[default_solver](model)
