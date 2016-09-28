"""
Unit testing module for dFBA and MdFBA.
"""
__author__ = 'kaizhuang'

import unittest

from framed.io_utils.sbml import load_cbmodel
from framed.bioreactor.base import *
from framed.analysis.simulation import FBA

SMALL_TEST_MODEL = '../../../examples/models/ecoli_core_model.xml'
ec_core_model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')


class OrganismTest(unittest.TestCase):
    def setUp(self):
        self.ec1 = Organism(ec_core_model)

    def testInitialization(self):
        self.assertNotEqual(ec_core_model, self.ec1.model)
        self.assertEqual(ec_core_model.id, self.ec1.model.id)
        self.assertListEqual(ec_core_model.metabolites.keys(), self.ec1.model.metabolites.keys())
        self.assertListEqual(ec_core_model.reactions.keys(), self.ec1.model.reactions.keys())
        self.assertListEqual(ec_core_model.stoichiometric_matrix(), self.ec1.model.stoichiometric_matrix())
        self.assertDictEqual(ec_core_model.bounds, self.ec1.model.bounds)
        self.assertEqual(self.ec1.fba_objective, {'R_Biomass_Ecoli_core_w_GAM': 1})

    def testBoundChanges(self):
        ec_core_model.bounds['R_EX_glc_e'] = (-15, 0)
        self.assertNotEqual(ec_core_model.bounds['R_EX_glc_e'], self.ec1.model.bounds['R_EX_glc_e'])
        self.ec1.model.bounds['R_EX_glc_e'] = (-15, 0)
        self.assertEqual(ec_core_model.bounds['R_EX_glc_e'], self.ec1.model.bounds['R_EX_glc_e'])

    def testUpdate(self):
        self.assertRaises(NotImplementedError, self.ec1.update)
        self.ec1.update = updateOrganism
        self.assertTrue(self.ec1.update(self.ec1) == self.ec1)

    def testFBA(self):
        correct_solution = FBA(ec_core_model)

        solution1 = FBA(self.ec1.model)
        self.assertTrue(solution1.status)
        self.assertEqual(correct_solution.fobj, solution1.fobj)

        solver = solver_instance(self.ec1.model)
        obj = {self.ec1.model.detect_biomass_reaction(): 1}
        solution2 = solver.solve_lp(obj, minimize=False)
        self.assertTrue(solution2.status)
        self.assertEqual(correct_solution.fobj, solution2.fobj)

    def tearDown(self):
        del self.ec1


class EnvironmentTest(unittest.TestCase):
    def setUp(self):
        self.o1 = Organism(ec_core_model)
        self.o2 = Organism(ec_core_model)
        self.env = Environment()

    def testInitialization(self):
        assert type(self.env) == Environment

    def test_add_organisms(self):
        self.env.add_organism(self.o1)
        assert self.env.organisms == [self.o1]
        self.env.add_organism(self.o2)
        assert self.env.organisms == [self.o1, self.o2]
        self.env.add_organisms([self.o1, self.o2])
        assert self.env.organisms == [self.o1, self.o2, self.o1, self.o2]

    def test_add_metabolites(self):
        self.env.add_metabolite('EX_glc(e)')
        assert self.env.metabolites == ['EX_glc(e)']
        self.env.add_metabolites(['EX_ac(e)', 'EX_o2(e)'])
        assert self.env.metabolites == ['EX_glc(e)', 'EX_ac(e)', 'EX_o2(e)']

    def tearDown(self):
        del self.o1
        del self.o2
        del self.env


class BioreactorTest(unittest.TestCase):
    def setUp(self):
        self.o1 = Organism(ec_core_model)
        self.o2 = Organism(ec_core_model)
        self.br = Bioreactor([self.o1, self.o2], ['R_EX_glc_e', 'R_EX_ac_e', 'R_EX_o2_e'])
        self.br2 = Bioreactor()

    def testInitialization(self):
        assert self.br.organisms == [self.o1, self.o2]
        assert self.br.metabolites == ['R_EX_glc_e', 'R_EX_ac_e', 'R_EX_o2_e']
        assert len(self.br.organisms) == 2

    def test_set_organisms(self):
        self.br2.set_organisms([self.o1, self.o2])
        assert self.br2.organisms == self.br.organisms

    def test_set_metabolites(self):
        self.br2.set_metabolites(['R_EX_glc_e', 'R_EX_ac_e', 'R_EX_o2_e'])
        assert self.br2.metabolites == self.br.metabolites

    def test_set_initial_conditions(self):
        self.br.set_initial_conditions([1], [0.1, 0.1], [10, 1, 0])
        self.assertEqual(self.br.initial_conditions, [1, 0.1, 0.1, 10, 1, 0])
        self.assertRaises(AssertionError, self.br.set_initial_conditions, 1, [0.1, 0.1], [10, 1, 0])

    def tearDown(self):
        del self.br
        del self.o1
        del self.o2


def updateOrganism(self):
    """
    fixture function for testing Organism.update()
    """
    return self


def suite():
    tests = [OrganismTest, EnvironmentTest, BioreactorTest]

    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
