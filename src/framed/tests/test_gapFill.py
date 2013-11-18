'''
Unit testing solver interfaces.

@author: Marta Matos
'''
import unittest
import re
from framed.io_utils.sbml import load_sbml_model, CONSTRAINT_BASED
from framed.solvers import *
from framed.analysis.simulation import FBA
from framed.core.fixes import fix_bigg_model
from framed.analysis.gapFind import *
from framed.analysis.gapFill import GapFill
from framed.analysis.variability import *

class GapFillTest(unittest.TestCase):

    #@unittest.skip("demonstrating skipping")
    def test_gapFind_gurobi(self):
        
        solver = GurobiSolver()
        model = '../../../examples/models/toy_model_holes.xml'
        reactionsDB = '../../../examples/models/toy_model_DB'
        
        self.model = load_sbml_model(model, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)
        
        framed_result = GapFill(self.model, reactionsDB, solver)
    
    #@unittest.skip("demonstrating skipping")
    def test_gapFind_glpk(self):

        solver = GlpkSolver()
        model = '../../../examples/models/toy_model_holes.xml'
        reactionsDB = '../../../examples/models/toy_model_DB'

        self.model = load_sbml_model(model, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)

        framed_result = GapFill(self.model, reactionsDB, solver)


def suite():
    tests = [GapFillTest]
    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
