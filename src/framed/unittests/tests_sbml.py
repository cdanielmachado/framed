import unittest

from framed import load_cbmodel, load_odemodel, save_sbml_model
from framed import read_model_from_file, write_model_to_file


SMALL_TEST_MODEL = '../../../examples/models/ecoli_core_model.xml'
TEST_MODEL_FBC = '../../../examples/models/ecoli_core_bigg2.0.xml'
TEST_MODEL_COPY = '../../../examples/models/cbmodel_copy.xml'
TEST_MODEL_COPY2 = '../../../examples/models/cbmodel_copy2.xml'
TEST_MODEL_COPY3 = '../../../examples/models/cbmodel_copy3.xml'
PLAIN_TEXT_COPY = '../../../examples/models/cbmodel_copy.txt'
KINETIC_MODEL = '../../../examples/models/BIOMD0000000051.xml'
KINETIC_MODEL_COPY = '../../../examples/models/odemodel_copy.xml'


class COBRA2COBRATest(unittest.TestCase):
    """ Test SBML import and export. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        save_sbml_model(model, TEST_MODEL_COPY, flavor='cobra')
        model_copy = load_cbmodel(TEST_MODEL_COPY, flavor='cobra')
        self.assertEqual(model.id, model_copy.id)
        self.assertListEqual(model.compartments.keys(), model_copy.compartments.keys())
        self.assertListEqual(model.metabolites.keys(), model_copy.metabolites.keys())
        self.assertListEqual(model.reactions.keys(), model_copy.reactions.keys())
        for r1, r2 in zip(model.reactions.values(), model_copy.reactions.values()):
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.reversible, r2.reversible)
            self.assertDictEqual(r1.stoichiometry, r2.stoichiometry)
            self.assertEqual(r1.lb, r2.lb)
            self.assertEqual(r1.ub, r2.ub)
            self.assertEqual(str(r1.gpr), str(r2.gpr))
        self.assertSetEqual(set(model.genes.keys()), set(model_copy.genes.keys()))


class FBC2FBCTest(unittest.TestCase):
    """ Test SBML import and export. """

    def testRun(self):
        model = load_cbmodel(TEST_MODEL_FBC, flavor='fbc2')
        save_sbml_model(model, TEST_MODEL_COPY2, flavor='fbc2')
        model_copy = load_cbmodel(TEST_MODEL_COPY2, flavor='fbc2')
        self.assertEqual(model.id, model_copy.id)
        self.assertListEqual(model.compartments.keys(), model_copy.compartments.keys())
        self.assertListEqual(model.metabolites.keys(), model_copy.metabolites.keys())
        self.assertListEqual(model.reactions.keys(), model_copy.reactions.keys())
        for r1, r2 in zip(model.reactions.values(), model_copy.reactions.values()):
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.reversible, r2.reversible)
            self.assertDictEqual(r1.stoichiometry, r2.stoichiometry)
            self.assertEqual(r1.lb, r2.lb)
            self.assertEqual(r1.ub, r2.ub)
            self.assertEqual(str(r1.gpr), str(r2.gpr))
        self.assertListEqual(model.genes.keys(), model_copy.genes.keys())



class COBRA2FBCTest(unittest.TestCase):
    """ Test SBML import and export. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')
        save_sbml_model(model, TEST_MODEL_COPY3, flavor='fbc2')
        model_copy = load_cbmodel(TEST_MODEL_COPY3, flavor='fbc2')
        self.assertEqual(model.id, model_copy.id)
        self.assertListEqual(model.compartments.keys(), model_copy.compartments.keys())
        self.assertListEqual(model.metabolites.keys(), model_copy.metabolites.keys())
        self.assertListEqual(model.reactions.keys(), model_copy.reactions.keys())
        for r1, r2 in zip(model.reactions.values(), model_copy.reactions.values()):
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.reversible, r2.reversible)
            self.assertDictEqual(r1.stoichiometry, r2.stoichiometry)
            self.assertEqual(r1.lb, r2.lb)
            self.assertEqual(r1.ub, r2.ub)
            self.assertEqual(str(r1.gpr), str(r2.gpr))
        self.assertListEqual(model.genes.keys(), model_copy.genes.keys())


class PlainTextIOTest(unittest.TestCase):
    """ Test plain text import and export. """

    def testRun(self):
        model = load_cbmodel(SMALL_TEST_MODEL)
        write_model_to_file(model, PLAIN_TEXT_COPY)
        model_copy = read_model_from_file(PLAIN_TEXT_COPY, kind='cb')
        self.assertListEqual(sorted(model.metabolites.keys()),
                             sorted(model_copy.metabolites.keys()))
        self.assertListEqual(model.reactions.keys(), model_copy.reactions.keys())


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
            self.assertDictEqual(r1.regulators, r2.regulators)
        self.assertDictEqual(model.ratelaws, model_copy.ratelaws)
        self.assertDictEqual(model.constant_params, model_copy.constant_params)
        self.assertDictEqual(model.variable_params, model_copy.variable_params)
        for p1, p2 in zip(model.local_params.values(), model_copy.local_params.values()):
            self.assertDictEqual(p1, p2)



def suite():
    tests = [COBRA2COBRATest, COBRA2FBCTest, FBC2FBCTest, PlainTextIOTest, SBMLTestODE]

    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())