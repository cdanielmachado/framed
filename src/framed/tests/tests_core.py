"""
Unit testing module for core features.

@author: Daniel Machado
"""
import unittest

from framed.io_utils.sbml import load_cbmodel, load_odemodel, save_sbml_model
from framed.analysis.simulation import FBA, pFBA, qpFBA
from framed.analysis.variability import FVA
from framed.io_utils.plaintext import read_model_from_file, write_model_to_file
from framed.analysis.deletion import gene_deletion
from framed.analysis.essentiality import essential_genes
from framed.solvers.solver import Status
from framed.core.transformation import make_irreversible, simplify, add_ratio_constraint


SMALL_TEST_MODEL = '../../../examples/models/ecoli_core_model.xml'
LARGE_TEST_MODEL = '../../../examples/models/Ec_iAF1260_flux1.xml'
TEST_MODEL_COPY = '../../../examples/models/cbmodel_copy.xml'
PLAIN_TEXT_COPY = '../../../examples/models/cbmodel_copy.txt'
KINETIC_MODEL = '../../../examples/models/BIOMD0000000051.xml'
KINETIC_MODEL_COPY = '../../../examples/models/odemodel_copy.xml'

GROWTH_RATE = 0.8739

DOUBLE_GENE_KO = ['b3731', 's0001']
DOUBLE_KO_GROWTH_RATE = 0.108
DOUBLE_KO_SUCC_EX = 3.8188

MOMA_GENE_KO = ['b0721']
MOMA_GROWTH_RATE = 0.5745
MOMA_SUCC_EX = 4.467

LMOMA_GENE_KO = ['b0721']
LMOMA_GROWTH_RATE = 0.5066
LMOMA_SUCC_EX = 5.311

ROOM_GENE_KO = ['b0721']
ROOM_GROWTH_RATE = 0.3120
ROOM_SUCC_EX = 2.932

ESSENTIAL_GENES = ['b0720', 'b1136', 'b1779', 'b2415', 'b2416', 'b2779', 'b2926']


class SBMLTest(unittest.TestCase):
    """ Test SBML import and export. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL)
        save_sbml_model(model, TEST_MODEL_COPY)
        model_copy = load_cbmodel(TEST_MODEL_COPY)
        self.assertEqual(model.id, model_copy.id)
        self.assertListEqual(model.compartments.keys(), model_copy.compartments.keys())
        self.assertListEqual(model.metabolites.keys(), model_copy.metabolites.keys())
        self.assertListEqual(model.reactions.keys(), model_copy.reactions.keys())
        for r1, r2 in zip(model.reactions.values(), model_copy.reactions.values()):
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.reversible, r2.reversible)
            self.assertDictEqual(r1.stoichiometry, r2.stoichiometry)
            self.assertDictEqual(model.bounds, model_copy.bounds)
        self.assertListEqual(model.genes.keys(), model_copy.genes.keys())
        self.assertDictEqual(model.rules, model_copy.rules)


class PlainTextIOTest(unittest.TestCase):
    """ Test plain text import and export. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL)
        write_model_to_file(model, PLAIN_TEXT_COPY)
        model_copy = read_model_from_file(PLAIN_TEXT_COPY, kind='cb')
        self.assertListEqual(sorted(model.metabolites.keys()),
                             sorted(model_copy.metabolites.keys()))
        self.assertListEqual(model.reactions.keys(), model_copy.reactions.keys())
        self.assertDictEqual(model.bounds, model_copy.bounds)


class FBATest(unittest.TestCase):
    """ Test FBA simulation. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        solution = FBA(model, get_shadow_prices=True, get_reduced_costs=True)
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.fobj, GROWTH_RATE, places=2)

class FBAwithRatioTest(unittest.TestCase):
    """ Test FBA with ratio constraints simulation. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        r_id1 = 'R_PGI'
        r_id2 = 'R_G6PDH2r'
        ratio = 2.0
        add_ratio_constraint(model, r_id1, r_id2, ratio)
        solution = FBA(model, get_shadow_prices=True, get_reduced_costs=True)
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertEqual(solution.values[r_id1] / solution.values[r_id2], ratio)

class pFBATest(unittest.TestCase):
    """ Test pFBA simulation. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        solution1 = pFBA(model)
        solution2 = FBA(model)
        self.assertEqual(solution1.status, Status.OPTIMAL)
        self.assertEqual(solution2.status, Status.OPTIMAL)
        growth1 = solution1.values[model.detect_biomass_reaction()]
        growth2 = solution2.values[model.detect_biomass_reaction()]
        self.assertAlmostEqual(growth1, growth2, places=4)
        norm1 = sum([abs(solution1.values[r_id]) for r_id in model.reactions])
        norm2 = sum([abs(solution2.values[r_id]) for r_id in model.reactions])
        self.assertLessEqual(norm1, norm2 + 1e-6)

class qpFBATest(unittest.TestCase):
    """ Test qpFBA simulation. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        solution1 = qpFBA(model)
        solution2 = FBA(model)
        self.assertEqual(solution1.status, Status.OPTIMAL)
        self.assertEqual(solution2.status, Status.OPTIMAL)
        growth1 = solution1.values[model.detect_biomass_reaction()]
        growth2 = solution2.values[model.detect_biomass_reaction()]
        self.assertAlmostEqual(growth1, growth2, places=4)
        norm1 = sum([abs(solution1.values[r_id]) for r_id in model.reactions])
        norm2 = sum([abs(solution2.values[r_id]) for r_id in model.reactions])
        self.assertLessEqual(norm1, norm2 + 1e-5)

class FBAFromPlainTextTest(unittest.TestCase):
    """ Test FBA simulation from plain text model. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        write_model_to_file(model, PLAIN_TEXT_COPY)
        model_copy = read_model_from_file(PLAIN_TEXT_COPY, kind='cb')
        solution = FBA(model_copy)
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.fobj, GROWTH_RATE, places=2)


class IrreversibleModelFBATest(unittest.TestCase):
    """ Test FBA simulation after reversible decomposition. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        make_irreversible(model)
        self.assertTrue(all([not reaction.reversible for reaction in model.reactions.values()]))
        solution = FBA(model)
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.fobj, GROWTH_RATE, places=2)


class SimplifiedModelFBATest(unittest.TestCase):
    """ Test FBA simulation after model simplification. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        simplify(model)
        solution = FBA(model)
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.fobj, GROWTH_RATE, places=2)


class TransformationCommutativityTest(unittest.TestCase):
    """ Test commutativity between transformations. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        simplify(model)
        make_irreversible(model)
        simplify(model) #remove directionally blocked reactions

        model2 = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        make_irreversible(model2)
        simplify(model2)

        self.assertEqual(model.id, model2.id)
        self.assertListEqual(model.metabolites.keys(), model2.metabolites.keys())
        self.assertListEqual(model.reactions.keys(), model2.reactions.keys())
        self.assertDictEqual(model.bounds, model2.bounds)
        self.assertListEqual(model.genes.keys(), model2.genes.keys())
        self.assertDictEqual(model.rules, model2.rules)


class FVATest(unittest.TestCase):
    """ Test flux variability analysis """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        variability = FVA(model)
        self.assertTrue(all([lb <= ub if lb is not None and ub is not None else True
                             for lb, ub in variability.values()]))


class GeneDeletionFBATest(unittest.TestCase):
    """ Test gene deletion with FBA. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        solution = gene_deletion(model, DOUBLE_GENE_KO)
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.values[model.detect_biomass_reaction()], DOUBLE_KO_GROWTH_RATE, 3)
        self.assertAlmostEqual(solution.values['R_EX_succ_e'], DOUBLE_KO_SUCC_EX, 3)


class GeneDeletionMOMATest(unittest.TestCase):
    """ Test gene deletion with MOMA. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        solution = gene_deletion(model, MOMA_GENE_KO, 'MOMA')
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.values[model.detect_biomass_reaction()], MOMA_GROWTH_RATE, 3)
        self.assertAlmostEqual(solution.values['R_EX_succ_e'], MOMA_SUCC_EX, 3)


class GeneDeletionLMOMATest(unittest.TestCase):
    """ Test gene deletion with MOMA. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        solution = gene_deletion(model, LMOMA_GENE_KO, 'lMOMA')
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.values[model.detect_biomass_reaction()], LMOMA_GROWTH_RATE, 3)
        self.assertAlmostEqual(solution.values['R_EX_succ_e'], LMOMA_SUCC_EX, 3)

class GeneDeletionROOMTest(unittest.TestCase):
    """ Test gene deletion with ROOM. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        solution = gene_deletion(model, ROOM_GENE_KO, 'ROOM')
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.values[model.detect_biomass_reaction()], ROOM_GROWTH_RATE, 3)
        self.assertAlmostEqual(solution.values['R_EX_succ_e'], ROOM_SUCC_EX, 3)


class GeneEssentialityTest(unittest.TestCase):
    """ Test gene essentiality. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='bigg')
        essential = essential_genes(model)
        self.assertListEqual(essential, ESSENTIAL_GENES)



class SBMLTestODE(unittest.TestCase):
    """ Test SBML import and export. """

    def testRun(self):
        model = load_odemodel(KINETIC_MODEL)
        save_sbml_model(model, KINETIC_MODEL_COPY)
        model_copy = load_odemodel(KINETIC_MODEL_COPY)
        self.assertEqual(model.id, model_copy.id)
        self.assertListEqual(model.compartments.keys(), model_copy.compartments.keys())
        for c1, c2 in zip(model.compartments.values(), model_copy.compartments.values()):
            self.assertEqual(c1.size, c2.size)
        self.assertListEqual(model.metabolites.keys(), model_copy.metabolites.keys())
        self.assertListEqual(model.reactions.keys(), model_copy.reactions.keys())
        for r1, r2 in zip(model.reactions.values(), model_copy.reactions.values()):
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.reversible, r2.reversible)
            self.assertDictEqual(r1.stoichiometry, r2.stoichiometry)
            self.assertListEqual(r1.modifiers, r2.modifiers)
        self.assertDictEqual(model.ratelaws, model_copy.ratelaws)
        self.assertDictEqual(model.global_parameters, model_copy.global_parameters)
        for p1, p2 in zip(model.local_parameters.values(), model_copy.local_parameters.values()):
            self.assertDictEqual(p1, p2)




def suite():
    tests = [SBMLTest, PlainTextIOTest, FBATest, pFBATest, qpFBATest, FBAFromPlainTextTest, FVATest, IrreversibleModelFBATest,
             SimplifiedModelFBATest, TransformationCommutativityTest, GeneDeletionFBATest, GeneDeletionMOMATest,
             GeneEssentialityTest, GeneDeletionLMOMATest, GeneDeletionROOMTest, SBMLTestODE]
#    tests = [SBMLTestODE]

    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())