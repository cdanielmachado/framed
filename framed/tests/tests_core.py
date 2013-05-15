'''
Unit testing module.
'''
import unittest

from framed.io_utils.sbml import load_sbml_model, save_sbml_model, CONSTRAINT_BASED, GPR_CONSTRAINED
from framed.core.fixes import fix_bigg_model
from framed.analysis.simulation import FBA, MOMA
from framed.analysis.variability import FVA
from framed.io_utils.plaintext import read_model_from_file, write_model_to_file
from framed.analysis.deletion import gene_deletion
from framed.analysis.essentiality import essential_genes

SMALL_TEST_MODEL = '../../misc/ecoli_core_model.xml'
TEST_MODEL_COPY = '../../misc/model_copy.xml'
PLAIN_TEXT_COPY = '../../misc/model_copy.txt'

GROWTH_RATE = 0.8739

DOUBLE_GENE_KO = ['b3731', 's0001']
DOUBLE_KO_GROWTH_RATE = 0.108
DOUBLE_KO_SUCC_EX = 3.8188

MOMA_GENE_KO = ['b0721']
MOMA_GROWTH_RATE = 0.5745
MOMA_SUCC_EX = 4.467

ESSENTIAL_GENES = ['b0720', 'b1136', 'b1779', 'b2415', 'b2416', 'b2779', 'b2926']


class SBMLTest(unittest.TestCase):
    """ Test SBML import and export. """
        
    def testRun(self):
        model = load_sbml_model(SMALL_TEST_MODEL, kind=GPR_CONSTRAINED)
        save_sbml_model(model, TEST_MODEL_COPY)
        model_copy = load_sbml_model(TEST_MODEL_COPY, kind=GPR_CONSTRAINED)
        self.assertEqual(model.id, model_copy.id)
        self.assertListEqual(model.metabolites.keys(), model_copy.metabolites.keys())
        self.assertListEqual(model.reactions.keys(), model_copy.reactions.keys())
        self.assertDictEqual(model.stoichiometry, model_copy.stoichiometry)
        self.assertDictEqual(model.bounds, model_copy.bounds)
        self.assertListEqual(model.genes.keys(), model_copy.genes.keys())
        self.assertDictEqual(model.rules, model_copy.rules)

class PlainTextIOTest(unittest.TestCase):
    """ Test plain text import and export. """
        
    def testRun(self):
        model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
        fix_bigg_model(model)
        write_model_to_file(model, PLAIN_TEXT_COPY)
        model_copy = read_model_from_file(PLAIN_TEXT_COPY, kind=CONSTRAINT_BASED)
        self.assertListEqual(sorted(model.metabolites.keys()),
                             sorted(model_copy.metabolites.keys()))
        self.assertListEqual(model.reactions.keys(), model_copy.reactions.keys())
        self.assertDictEqual(dict(model.stoichiometry),
                             dict(model_copy.stoichiometry))
        self.assertDictEqual(model.bounds, model_copy.bounds)
        
class FBATest(unittest.TestCase):
    """ Test FBA simulation. """
    
    def testRun(self):
        model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
        fix_bigg_model(model)
        solution = FBA(model)
        self.assertTrue(solution.status)
        self.assertAlmostEqual(solution.fobj, GROWTH_RATE, places=2)

class FBATest2(unittest.TestCase):
    """ Test FBA simulation from plain text model. """
    
    def testRun(self):
        model = read_model_from_file(PLAIN_TEXT_COPY, kind=CONSTRAINT_BASED)
        solution = FBA(model)
        self.assertTrue(solution.status)
        self.assertAlmostEqual(solution.fobj, GROWTH_RATE, places=2)
        
class FVATest(unittest.TestCase):
    """ Test flux variability analysis """
    
    def testRun(self):
        model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
        fix_bigg_model(model)
        variability = FVA(model)        
        self.assertTrue(all([lb <= ub if lb is not None and ub is not None else True
                             for lb, ub in variability.values()]))


class GeneDeletionFBATest(unittest.TestCase):
    """ Test gene deletion with FBA. """
    
    def testRun(self):
        model = load_sbml_model(SMALL_TEST_MODEL, kind=GPR_CONSTRAINED)
        fix_bigg_model(model)
        solution = gene_deletion(model, DOUBLE_GENE_KO)
        self.assertTrue(solution.status)
        biomass = model.reactions.keys().index(model.detect_biomass_reaction())
        succ_ex = model.reactions.keys().index('R_EX_succ_e')
        self.assertAlmostEqual(solution.values[biomass], DOUBLE_KO_GROWTH_RATE, 3)
        self.assertAlmostEqual(solution.values[succ_ex], DOUBLE_KO_SUCC_EX, 3)

class GeneDeletionMOMATest(unittest.TestCase):
    """ Test gene deletion with MOMA. """
    
    def testRun(self):
        model = load_sbml_model(SMALL_TEST_MODEL, kind=GPR_CONSTRAINED)
        fix_bigg_model(model)
        solution = gene_deletion(model, MOMA_GENE_KO, 'MOMA')
        self.assertTrue(solution.status)
        biomass = model.reactions.keys().index(model.detect_biomass_reaction())
        succ_ex = model.reactions.keys().index('R_EX_succ_e')
        self.assertAlmostEqual(solution.values[biomass], MOMA_GROWTH_RATE, 3)
        self.assertAlmostEqual(solution.values[succ_ex], MOMA_SUCC_EX, 3)
                
class GeneEssentialityTest(unittest.TestCase):
    """ Test gene deletion with MOMA. """
    
    def testRun(self):
        model = load_sbml_model(SMALL_TEST_MODEL, kind=GPR_CONSTRAINED)
        fix_bigg_model(model)
        essential = essential_genes(model)
        self.assertListEqual(essential, ESSENTIAL_GENES)
        
        
                
def suite():
    #tests = [SBMLTest, PlainTextIOTest, FBATest, FVATest, FBATest2, GeneDeletionFBATest, GeneDeletionMOMATest]
    tests = [GeneEssentialityTest]
    
    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite
        

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())