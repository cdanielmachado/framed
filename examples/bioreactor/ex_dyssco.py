__author__ = 'kaizhuang'

import matplotlib.pyplot as plt

from framed.io.sbml import load_sbml_model, CONSTRAINT_BASED
from framed.io.plaintext import add_reaction_from_str
from framed.model.fixes import fix_cobra_model
from framed.bioreactor import *

import framed.bioreactor.dyssco as dyssco


# fixture class
class Ecoli(Organism):

    def update(self):
        BR = self.environment

        # calculating and updating the glucose uptake constraint
        rid = BR.metabolites.index('R_EX_glc_e')
        vlb_glc = float(-10 * BR.S[rid] / (BR.S[rid] + 1))
        self.fba_constraints['R_EX_glc_e'] = (vlb_glc, 0)
        # ideal aerobic condition
        self.fba_constraints['R_EX_o2_e'] = (-15, None)



### Setup
EC_1260_MODEL = 'models/Ec_iAF1260_gene_names.xml'
ec1260 = load_sbml_model(EC_1260_MODEL, kind=CONSTRAINT_BASED)
fix_cobra_model(ec1260)

# Add 1-3-PDO pathway
add_reaction_from_str(ec1260, 'R_glycerol_dehydratase: M_glyc_c + M_nadh_c + M_h_c --> M_1_3_pdo_c + M_nad_c')
add_reaction_from_str(ec1260, 'R_1_3_PDOtr: M_1_3_pdo_c <-> M_1_3_pdo_e')
add_reaction_from_str(ec1260, 'R_EX_1_3_pdo_e: M_1_3_pdo_e -->')


ec = Ecoli(ec1260, id='ecoli')

# make envelope strains
r_substrate, r_target1, r_target2, r_oxygen = 'R_EX_glc_e', 'R_EX_1_3_pdo_e', 'R_EX_3_hp_e', 'R_EX_o2_e'
aerobic_constraints = {r_substrate: (-10, 0), r_oxygen: (-15, None)}
strains = dyssco.make_envelope_strains(ec, 'R_EX_glc_e', 'R_EX_1_3_pdo_e', N=5, constraints=aerobic_constraints)

# define bioreactor
br = IdealFedbatch(metabolites=['R_EX_glc_e', 'R_EX_1_3_pdo_e'], Sfeed=[1000, 0], volume_max=10)
Vinit = [1]
Xinit = [0.01]
Sinit = [10, 0]
br.set_initial_conditions(Vinit, Xinit, Sinit)

#calculating performances

performances = dyssco.calculate_performances(strains, br, 'R_EX_glc_e', 'R_EX_1_3_pdo_e', 0, 20, 0.1, verbose=True,
                                             additional_yields=['R_EX_co2_e', 'R_EX_nh4_e'])


metrics = dyssco.performances2metrics(performances)
