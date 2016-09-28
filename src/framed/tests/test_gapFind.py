"""
Unit testing for gap finding.

@author: Marta Matos
"""

import re
import unittest

from framed.core.fixes import fix_cobra_model
from framed.io_utils.plaintext import read_model_from_file
from framed.io_utils.sbml import load_cbmodel
from framed.solvers.glpk_wrapper_lazy import GlpkSolverLazy
from framed.reconstruction.gapFind import GapFind
from framed.solvers.gurobi_wrapper import GurobiSolver

class GapFindTest(unittest.TestCase):

    """ The tests in this class find gaps in 7 networks and
    the results are compared to the ones obtained with the
    same algorithm implemented by the COBRA toolbox.
    The comparison is done separately for Glpk and Gurobi,
    as each solver produces different results.
    COBRA results are stored in files.
    """

    def compare_to_COBRA(self, framed_result, COBRA_mets_result_filename, COBRA_rxns_result_filename):
        """ Compares framed results to COBRA results
        """
        mets_file = open(COBRA_mets_result_filename, 'r')
        rxns_file = open(COBRA_rxns_result_filename, 'r')

        COBRA_mets = mets_file.read()
        COBRA_rxns = rxns_file.read()

        COBRA_mets = COBRA_mets.split()
        COBRA_rxns = COBRA_rxns.split()

        for i in range(0, len(COBRA_mets)):
            COBRA_mets[i] = 'M_' + re.sub('[\[]', '_', COBRA_mets[i])
            COBRA_mets[i] = re.sub('[\]]', '', COBRA_mets[i])

        for i in range(0, len(COBRA_rxns)):
            COBRA_rxns[i] = 'R_' + re.sub('[(]', '_', COBRA_rxns[i])
            COBRA_rxns[i] = re.sub('[\)]', '', COBRA_rxns[i])

        framed_result[1].sort()

        mets_file.close()
        rxns_file.close()

        return (framed_result[0] == COBRA_mets and framed_result[1] == COBRA_rxns)

#===============================================================================
#     def test_gapFind_glpk(self):
#         """ Finds gaps in 7 models using Glpk
#         """
#         for i in range(1,8):
#             solver = GlpkSolver()
# 
#             model = '../../../examples/reconstruction/gapFind/ecoli_core_model' + str(i) + '.xml'
#             self.model = load_sbml_model(model, kind=CONSTRAINT_BASED)
#             fix_cobra_model(self.model)
# 
#             framed_result = GapFind(self.model, solver)
#             COBRA_mets_result_filename = '../../../examples/reconstruction/gapFind/gapFind_results/ecoli0' + str(i) + '_mets_glpk.txt'
#             COBRA_rxns_result_filename = '../../../examples/reconstruction/gapFind/gapFind_results/ecoli0' + str(i) + '_rxns_glpk.txt'
# 
#             result = self.compare_to_COBRA(framed_result, COBRA_mets_result_filename, COBRA_rxns_result_filename)
#             self.assertTrue(result)
#===============================================================================

    def test_gapFind_glpk_lazy(self):
        """ Finds gaps in 7 models using Glpk
        """
        for i in range(1,8):
            solver = GlpkSolverLazy(tol_int="default")
 
            model = '../../../examples/reconstruction/gapFind/ecoli_core_model' + str(i) + '.xml'
            self.model = load_cbmodel(model, flavor='cobra')

            framed_result = GapFind(self.model, solver)
            COBRA_mets_result_filename = '../../../examples/reconstruction/gapFind/gapFind_results/ecoli0' + str(i) + '_mets_glpk.txt'
            COBRA_rxns_result_filename = '../../../examples/reconstruction/gapFind/gapFind_results/ecoli0' + str(i) + '_rxns_glpk.txt'
 
            result = self.compare_to_COBRA(framed_result, COBRA_mets_result_filename, COBRA_rxns_result_filename)
            self.assertTrue(result)
            
            
    def test_gapFind_gurobi(self):
        """ Finds gaps in 7 models using Gurobi
        """
        for i in range(1,8):
            if i == 6:
                continue
            solver = GurobiSolver()

            model = '../../../examples/reconstruction/gapFind/ecoli_core_model' + str(i) + '.xml'
            self.model = load_cbmodel(model, flavor='cobra')

            framed_result = GapFind(self.model, solver)
            COBRA_mets_result_filename = '../../../examples/reconstruction/gapFind/gapFind_results/ecoli0' + str(i) + '_mets_gurobi.txt'
            COBRA_rxns_result_filename = '../../../examples/reconstruction/gapFind/gapFind_results/ecoli0' + str(i) + '_rxns_gurobi.txt'

            result = self.compare_to_COBRA(framed_result, COBRA_mets_result_filename, COBRA_rxns_result_filename)
            self.assertTrue(result)            

class GapFind_exp(unittest.TestCase):
    """ Class created just to look at the results from 
    the particular method that was tested
    """
    
    def test_gapFind_exp(self):
        solver = GlpkSolverLazy()
        model = '../../../examples/reconstruction/toy_model'
        self.model = read_model_from_file(model, kind='cb')
        fix_cobra_model(self.model)
        print self.model
        framed_result = GapFind(self.model, solver, root_gaps_only=True)
        print framed_result

def suite():
    tests = [GapFindTest]
    #GapFindTest
    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
