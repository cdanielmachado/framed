'''

@author: Daniel Machado

   Copyright 2013 Novo Nordisk Foundation Center for Biosustainability,
   Technical University of Denmark.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
   
'''

from framed.design.combinatorial import combinatorial_gene_deletion
from framed.io_utils.sbml import load_sbml_model, GPR_CONSTRAINED
from framed.core.fixes import fix_bigg_model
from time import time

SMALL_TEST_MODEL = '../../../examples/models/ecoli_core_model.xml'
LARGE_TEST_MODEL = '../../../examples/models/Ec_iAF1260_gene_names.xml'


def benchmark(method, model):
    print 'benchmarking', method,
    tstart = time()
    fobj = lambda v: v['R_EX_succ_e']
    max_dels = 1
    combinatorial_gene_deletion(model, fobj, max_dels, method=method)
    tend = time()
    print 'took', tend - tstart


def main():
    model = load_sbml_model(SMALL_TEST_MODEL, GPR_CONSTRAINED)
    #model = load_sbml_model(LARGE_TEST_MODEL, GPR_CONSTRAINED)
    fix_bigg_model(model)
    benchmark('FBA', model)
    benchmark('pFBA', model)
#    benchmark('qpFBA', model)
#    benchmark('MOMA', model)
#    benchmark('lMOMA', model)


if __name__ == '__main__':
    main()
