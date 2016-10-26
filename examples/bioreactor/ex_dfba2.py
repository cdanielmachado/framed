"""
Example: Dynamic Flux Balance Analysis of E. coli producing Acetate in Fedbatch

@Author: Kai Zhuang
"""
__author__ = 'kaizhuang'

import matplotlib.pyplot as plt

from framed.bioreactor import *
from framed.model.fixes import fix_cobra_model
from framed.io.sbml import load_sbml_model, CONSTRAINT_BASED


### Basic Setup
SMALL_TEST_MODEL = 'models/Ec_iAF1260_gene_names.xml'
ec1260 = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
fix_cobra_model(ec1260)


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

        # no acetate uptake.
        self.fba_constraints['R_EX_ac_e'] = (0, None)

        # ideal aerobic condition
        self.fba_constraints['R_EX_o2_e'] = (-15, None)


### Main Program
# creating an instance of Ecoli
ec = Ecoli(ec1260, id='ecoli')

# creating a batch bioreactor containing Ecoli, glucose, acetate, and oxygen
fedbatch_bioreactor = IdealFedbatch(ec, ['R_EX_glc_e', 'R_EX_ac_e'], Sfeed=[1000, 0], volume_max=10)

# set initial conditions
Vinit = [1]             # liquid volume
Xinit = [0.005]          # cell concentration
Sinit = [10, 0]      # concentrations of glucose, acetate
fedbatch_bioreactor.set_initial_conditions(Vinit, Xinit, Sinit)

# set simulation time interval
t0 = 0
tf = 20
dt = 0.1

# run dFBA simulation
result = dFBA(fedbatch_bioreactor, t0, tf, dt, solver='dopri5', verbose=True)


plt.figure(1)
plt.subplot(221)
plt.plot(result['time'], result['ecoli'])
plt.title('E. coli')
plt.subplot(222)
plt.plot(result['time'], result['R_EX_glc_e'])
plt.title('Glucose')
plt.subplot(223)
plt.plot(result['time'], result['R_EX_ac_e'])
plt.title('Acetate')
plt.subplot(224)
plt.plot(result['time'], result['volume'])
plt.title('Volume')
plt.show()



