"""
Example: Dynamic Flux Balance Analysis of Diauxic Growth of E. coli in Batch Bioreactor

"""
__author__ = 'kaizhuang'


import matplotlib.pyplot as plt

from framed.io_utils.sbml import load_sbml_model, CONSTRAINT_BASED
from framed.core.fixes import fix_bigg_model
from framed.analysis.dfba import *
from framed.bioreactor.bioreactor import *

### Basic Setup
SMALL_TEST_MODEL = '../misc/ecoli_core_model.xml'
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

        # calculating and updating the acetate uptake constraint
        rid = BR.metabolites.index('R_EX_ac_e')
        vlb_ac = float(-10 * BR.S[rid] / (BR.S[rid] + 1))
        self.fba_constraints['R_EX_ac_e'] = (vlb_ac, None)
        self.fba_constraints['R_EX_o2_e'] = (-10, None)

        # calculating and updating the oxygen uptake constraint
        rid = BR.metabolites.index('R_EX_o2_e')
        vlb_o2 = float(-15 * BR.S[rid] / (BR.S[rid] + 0.001))
        self.fba_constraints['R_EX_o2_e'] = (vlb_o2, None)


### Defining BatchBioreactor class
# In the Ecoli class, the method update(), which is an abstract method in the superclass Bioreactor, is defined.
# The update() method is used to re-calculate the oxygen transfer rate at each time step
class BatchBioreactor(Bioreactor):

    def update(self, t):
        kLa = 7.5
        rid = self.metabolites.index('R_EX_o2_e')
        deltaO2 = kLa*(0.21 - self.S[rid])
        self.deltaS[rid] = deltaO2

### Main Program
# creating an instance of Ecoli
ec = Ecoli(ec_core_model)

# creating a batch bioreactor containing Ecoli, glucose, acetate, and oxygen
batch_bioreactor = BatchBioreactor(ec, ['R_EX_glc_e', 'R_EX_ac_e', 'R_EX_o2_e'])

# set initial conditions
Vinit = [1]             # liquid volume
Xinit = [0.005]          # cell concentration
Sinit = [5, 0, 0.21]      # concentrations of glucose, acetate, and oxygen
batch_bioreactor.set_initial_conditions(Vinit, Xinit, Sinit)

# set simulation time interval
t0 = 0
tf = 10
dt = 0.1

# run dFBA simulation
t, y = dFBA(batch_bioreactor, t0, tf, dt, solver='dopri5', verbose=True)

plt.figure(1)
plt.subplot(221)
plt.plot(t, y[:, 1])
plt.title('E. coli')
plt.subplot(222)
plt.plot(t, y[:, 2])
plt.title('Glucose')
plt.subplot(223)
plt.plot(t, y[:, 3])
plt.title('Acetate')
plt.subplot(224)
plt.plot(t, y[:, 4])
plt.title('Oxygen')
plt.show()



