"""
Example: Flux Envelope Analysis

@Author: Kai Zhuang
"""
__author__ = 'kaizhuang'

from framed.analysis.variability import production_envelope
from framed.io_utils.sbml import load_sbml_model, CONSTRAINT_BASED
from framed.core.fixes import fix_bigg_model

import matplotlib.pyplot as plt

### Basic Setup
SMALL_TEST_MODEL = 'models/ecoli_core_model.xml'
ec_core_model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
fix_bigg_model(ec_core_model)

### Make the production envelope, and plot it
biomass, target_mins1, target_maxs1 = production_envelope(ec_core_model, 'R_EX_ac_e', steps=10)

biomass, target_mins2, target_maxs2 = production_envelope(ec_core_model, 'R_EX_succ_e', steps=10)

plt.figure(1)
plt.hold(True)
plt.plot(biomass, target_maxs1)
plt.plot(biomass, target_mins1)
plt.ylim(ymin=-0.1)
plt.show()