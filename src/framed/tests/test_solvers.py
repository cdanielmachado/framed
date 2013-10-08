'''
Unit testing solver interfaces.

@author: Nikolaus Sonnenschein
'''
import unittest
import pickle

from framed.io_utils.sbml import load_sbml_model, CONSTRAINT_BASED
from framed.solvers.gurobi_wrapper import GurobiSolver
from framed.solvers import solver_instance
from framed.analysis.simulation import FBA
from framed.core.fixes import fix_bigg_model

SMALL_TEST_MODEL = '../../../examples/models/ecoli_core_model.xml'


class SolverPickleTest(unittest.TestCase):
    """docstring for SolverPickleTest"""

    def setUp(self):
        self.model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)

    def test_gurobi_solver_pickle(self):
        self.solver = solver_instance()
        self.solver.build_problem(self.model)
        solver_from_pickle = pickle.loads(pickle.dumps(self.solver))
        self.assertEqual(FBA(self.model, solver=self.solver).status, FBA(self.model, solver=solver_from_pickle).status)


def suite():
    tests = [SolverPickleTest]
    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())

# model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
# fix_bigg_model(model)
# solver=solver_instance()
# solver.build_problem(model)



# import multiprocessing

# struct = {'model':model, 'solver':solver}

# def test_fba(kwargs):
# 	return FBA(struct['model'], solver=struct['solver']).fobj

# def f(x):
# 	return x*x

# pool = multiprocessing.Pool(4)
# print pool.map(f, [1,2,3])

# print pool.map(test_fba, [struct, struct])

