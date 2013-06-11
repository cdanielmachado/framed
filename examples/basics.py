from framed.io_utils.sbml import *
from framed.core.fixes import fix_bigg_model
from framed.analysis.simulation import FBA

ec_core_model = 'models/ecoli_core_model.xml'

#load SBML model, and fix BIGG model's perculiarities
model = load_sbml_model('models/ecoli_core_model.xml', kind=CONSTRAINT_BASED)
fix_bigg_model(model)

solution = FBA(model)

r_biomass = model.detect_biomass_reaction()

print solution.values[r_biomass]