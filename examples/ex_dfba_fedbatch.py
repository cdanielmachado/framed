"""
Example: Dynamic Flux Balance Analysis of Diauxic Growth of E. coli in Batch Bioreactor

@Author: Kai Zhuang
"""
__author__ = 'kaizhuang'


import matplotlib.pyplot as plt

from framed.io_utils.sbml import load_sbml_model, CONSTRAINT_BASED
from framed.core.fixes import fix_bigg_model
from framed.analysis.dfba import *
from framed.bioreactor.bioreactor import *


### Basic Setup
SMALL_TEST_MODEL = 'models/ecoli_core_model.xml'
ec_core_model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
fix_bigg_model(ec_core_model)


### Defining the Ecoli class
# In the Ecoli class, the method update(), which is an abstract method in the superclass Organism, is defined.
# The update() method describes how Ecoli will respond to changes in metabolite concentrations in its environment.
# Usually, the relevant FBA uptake constraints are calculated and updated in the update() method.
class Ecoli(Organism):

    def update(self):
        BR = self.environment

        # calculating and updating the glucose uptake constraint
        rid = BR.metabolites.index('R_EX_glc_e')
        vlb_glc = float(-10 * BR.S[rid] / (BR.S[rid] + 1))
        self.fba_constraints['R_EX_glc_e'] = (vlb_glc, 0)

        self.fba_constraints['R_EX_ac_e'] = (0, None)
        self.fba_constraints['R_EX_o2_e'] = (-15, None)


### Main Program
# creating an instance of Ecoli
ec = Ecoli(ec_core_model)


# creating a batch bioreactor containing Ecoli, glucose, acetate, and oxygen
batch_bioreactor = IdealFedbatch(ec, ['R_EX_glc_e', 'R_EX_ac_e'])
