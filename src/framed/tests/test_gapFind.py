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
from framed.analysis.gapFind_raven import GapFindRaven
from framed.analysis.variability import *

class GapFindTest(unittest.TestCase):

    #def setUp(self):
    #    self.model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
    #    fix_bigg_model(self.model)
    
    def compare_to_matlab(self, framed_result, matlab_mets_result_filename, matlab_rxns_result_filename):
        mets_file = open(matlab_mets_result_filename, 'r')
        rxns_file = open(matlab_rxns_result_filename, 'r')
        
        matlab_mets = mets_file.read()
        matlab_rxns = rxns_file.read()
        
        matlab_mets = matlab_mets.split()
        matlab_rxns = matlab_rxns.split()

        for i in range(0,len(matlab_mets)):
            matlab_mets[i] = 'M_' + re.sub('[\[]', '_', matlab_mets[i])
            matlab_mets[i] = re.sub('[\]]', '', matlab_mets[i])

        for i in range(0,len(matlab_rxns)):
            matlab_rxns[i] = 'R_' + re.sub('[(]', '_', matlab_rxns[i])
            matlab_rxns[i] = re.sub('[\)]', '', matlab_rxns[i])
        
        #print 'framed'
        #print framed_result[0]
        #print 'matlab'
        #print matlab_mets
        
        #print 'framed'
        framed_result[1].sort()
        #print framed_result[1]
        #print 'matlab'
        #print matlab_rxns

        mets_file.close()
        rxns_file.close()
        
        return (framed_result[0] == matlab_mets and framed_result[1] == matlab_rxns)

    @unittest.skip("demonstrating skipping")
    def test_gapFind_gurobi(self):
        
        for i in range(1,8):
            if i == 6:
                continue
            solver = GurobiSolver()
            print 'iter ' + str(i)
            model = '../../../examples/models/ecoli_core_model' + str(i) + '.xml'
            self.model = load_sbml_model(model, kind=CONSTRAINT_BASED)
            fix_bigg_model(self.model)

            framed_result = GapFind(self.model, solver)
            matlab_mets_result_filename = '../../../examples/models/gapFind_results/ecoli0' + str(i) + '_mets_gurobi.txt'
            matlab_rxns_result_filename = '../../../examples/models/gapFind_results/ecoli0' + str(i) + '_rxns_gurobi.txt'
            
            result = self.compare_to_matlab(framed_result, matlab_mets_result_filename, matlab_rxns_result_filename)
            print result
            self.assertTrue(result)
    
    #@unittest.skip("demonstrating skipping")
    def test_gapFind_glpk(self):

        for i in range(1,2):
            solver = GlpkSolver()
            print 'iter ' + str(i)
            model = '../../../examples/models/ecoli_core_model' + str(i) + '.xml'
            self.model = load_sbml_model(model, kind=CONSTRAINT_BASED)
            fix_bigg_model(self.model)
            
            framed_result = GapFind(self.model, solver)
            matlab_mets_result_filename = '../../../examples/models/gapFind_results/ecoli0' + str(i) + '_mets_glpk.txt'
            matlab_rxns_result_filename = '../../../examples/models/gapFind_results/ecoli0' + str(i) + '_rxns_glpk.txt'
            
            result = self.compare_to_matlab(framed_result, matlab_mets_result_filename, matlab_rxns_result_filename)
            print result
            self.assertTrue(result)
            

def suite():
    tests = [GapFindTest]
    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
