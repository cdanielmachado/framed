"""
Unit testing module for core features.

@author: Daniel Machado
"""

import unittest

from framed.io.sbml import load_cbmodel
from framed.cobra.simulation import FBA, pFBA
from framed.cobra.variability import FVA
from framed.io.plaintext import read_model_from_file, write_model_to_file
from framed.cobra.deletion import gene_deletion
from framed.cobra.essentiality import essential_genes
from framed.solvers.solver import Status
from framed.model.transformation import make_irreversible, simplify
from framed.solvers import set_default_solver

#set_default_solver('cplex')

SMALL_TEST_MODEL = '../../../examples/models/ecoli_core_model.xml'
LARGE_TEST_MODEL = '../../../examples/models/Ec_iAF1260_flux1.xml'
TEST_MODEL_COPY = '../../../examples/models/cbmodel_copy.xml'
PLAIN_TEXT_COPY = '../../../examples/models/cbmodel_copy.txt'
KINETIC_MODEL = '../../../examples/models/BIOMD0000000051.xml'
KINETIC_MODEL_COPY = '../../../examples/models/odemodel_copy.xml'

GROWTH_RATE = 0.8739

DOUBLE_GENE_KO = ['G_b3731', 'G_s0001']
DOUBLE_KO_GROWTH_RATE = 0.108
DOUBLE_KO_SUCC_EX = 3.8188

MOMA_GENE_KO = ['G_b0721']
MOMA_GROWTH_RATE = 0.5745
MOMA_SUCC_EX = 4.467

LMOMA_GENE_KO = ['G_b0721']
LMOMA_GROWTH_RATE = 0.5066
LMOMA_SUCC_EX = 5.311

ROOM_GENE_KO = ['G_b0721']
ROOM_GROWTH_RATE = 0.373
ROOM_SUCC_EX = 5.799

ESSENTIAL_GENES = ['G_b0720', 'G_b1136', 'G_b1779', 'G_b2415', 'G_b2416', 'G_b2779', 'G_b2926']


class FBATest(unittest.TestCase):
    """ Test FBA simulation. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        solution = FBA(model, get_shadow_prices=True, get_reduced_costs=True)
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.fobj, GROWTH_RATE, places=2)

class FBAwithRatioTest(unittest.TestCase):
    """ Test FBA with ratio constraints simulation. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        r_id1 = 'R_PGI'
        r_id2 = 'R_G6PDH2r'
        ratio = 2.0
        model.add_ratio_constraint(r_id1, r_id2, ratio)
        solution = FBA(model, get_shadow_prices=True, get_reduced_costs=True)
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertEqual(solution.values[r_id1] / solution.values[r_id2], ratio)

class pFBATest(unittest.TestCase):
    """ Test pFBA simulation. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        solution1 = pFBA(model)
        solution2 = FBA(model)
        self.assertEqual(solution1.status, Status.OPTIMAL)
        self.assertEqual(solution2.status, Status.OPTIMAL)
        growth1 = solution1.values[model.biomass_reaction]
        growth2 = solution2.values[model.biomass_reaction]
        self.assertAlmostEqual(growth1, growth2, places=4)
        norm1 = sum([abs(solution1.values[r_id]) for r_id in model.reactions])
        norm2 = sum([abs(solution2.values[r_id]) for r_id in model.reactions])
        self.assertLessEqual(norm1, norm2 + 1e-6)

class FBAFromPlainTextTest(unittest.TestCase):
    """ Test FBA simulation from plain text model. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        write_model_to_file(model, PLAIN_TEXT_COPY)
        model_copy = read_model_from_file(PLAIN_TEXT_COPY, kind='cb')
        solution = FBA(model_copy)
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.fobj, GROWTH_RATE, places=2)


class IrreversibleModelFBATest(unittest.TestCase):
    """ Test FBA simulation after reversible decomposition. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        make_irreversible(model)
        self.assertTrue(all([not reaction.reversible for reaction in model.reactions.values()]))
        solution = FBA(model)
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.fobj, GROWTH_RATE, places=2)


class SimplifiedModelFBATest(unittest.TestCase):
    """ Test FBA simulation after model simplification. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        simplify(model)
        solution = FBA(model)
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.fobj, GROWTH_RATE, places=2)


class TransformationCommutativityTest(unittest.TestCase):
    """ Test commutativity between transformations. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        simplify(model)
        make_irreversible(model)
        simplify(model) #remove directionally blocked reactions

        model2 = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        make_irreversible(model2)
        simplify(model2)

        self.assertEqual(model.id, model2.id)
        self.assertListEqual(model.metabolites.keys(), model2.metabolites.keys())
        self.assertListEqual(model.reactions.keys(), model2.reactions.keys())
        self.assertListEqual(model.genes.keys(), model2.genes.keys())
        for r1, r2 in zip(model.reactions.values(), model2.reactions.values()):
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.reversible, r2.reversible)
            self.assertDictEqual(r1.stoichiometry, r2.stoichiometry)
            self.assertEqual(r1.lb, r2.lb)
            self.assertEqual(r1.ub, r2.ub)
            self.assertEqual(str(r1.gpr), str(r2.gpr))
        self.assertSetEqual(set(model.genes.keys()), set(model2.genes.keys()))

class FVATest(unittest.TestCase):
    """ Test flux variability analysis """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        variability = FVA(model)
        self.assertTrue(all([lb <= ub if lb is not None and ub is not None else True
                             for lb, ub in variability.values()]))


class GeneDeletionpFBATest(unittest.TestCase):
    """ Test gene deletion with pFBA. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        solution = gene_deletion(model, DOUBLE_GENE_KO, 'pFBA')
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.values[model.biomass_reaction], DOUBLE_KO_GROWTH_RATE, 2)
        self.assertAlmostEqual(solution.values['R_EX_succ_e'], DOUBLE_KO_SUCC_EX, 2)


class GeneDeletionMOMATest(unittest.TestCase):
    """ Test gene deletion with MOMA. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        solution = gene_deletion(model, MOMA_GENE_KO, 'MOMA')
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.values[model.biomass_reaction], MOMA_GROWTH_RATE, 2)
        self.assertAlmostEqual(solution.values['R_EX_succ_e'], MOMA_SUCC_EX, 2)


class GeneDeletionLMOMATest(unittest.TestCase):
    """ Test gene deletion with MOMA. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        solution = gene_deletion(model, LMOMA_GENE_KO, 'lMOMA')
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.values[model.biomass_reaction], LMOMA_GROWTH_RATE, 2)
        self.assertAlmostEqual(solution.values['R_EX_succ_e'], LMOMA_SUCC_EX, 2)

class GeneDeletionROOMTest(unittest.TestCase):
    """ Test gene deletion with ROOM. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        solution = gene_deletion(model, ROOM_GENE_KO, 'ROOM')
        self.assertEqual(solution.status, Status.OPTIMAL)
        self.assertAlmostEqual(solution.values[model.biomass_reaction], ROOM_GROWTH_RATE, 2)
        self.assertAlmostEqual(solution.values['R_EX_succ_e'], ROOM_SUCC_EX, 2)


class GeneEssentialityTest(unittest.TestCase):
    """ Test gene essentiality. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        essential = essential_genes(model)
        self.assertListEqual(essential, ESSENTIAL_GENES)


def suite():
    tests = [FBATest, pFBATest, FBAFromPlainTextTest, FVATest, IrreversibleModelFBATest,
             SimplifiedModelFBATest, TransformationCommutativityTest, GeneDeletionpFBATest, GeneDeletionMOMATest,
             GeneEssentialityTest, GeneDeletionLMOMATest, GeneDeletionROOMTest]

    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())