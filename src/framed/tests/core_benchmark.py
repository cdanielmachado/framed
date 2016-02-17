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

from framed.analysis.deletion import reaction_deletion
from framed.design.combinatorial import combinatorial_gene_deletion
from framed.io_utils.sbml import load_cbmodel
from framed.solvers import solver_instance
from framed.analysis.simulation import FBA
from time import time

SMALL_TEST_MODEL = '../../../examples/models/ecoli_core_model.xml'
LARGE_TEST_MODEL = '../../../examples/models/Ec_iAF1260_flux1.xml'


def benchmark_combinatorial(method, model):
    print 'benchmarking', method,
    tstart = time()
    fobj = lambda v: v['R_EX_succ_e']
    max_dels = 1
    combinatorial_gene_deletion(model, fobj, max_dels, method=method)
    tend = time()
    print 'took', tend - tstart


def benchmark_methods_combinatorial(modelpath):
    model = load_cbmodel(modelpath, flavor='bigg')
    benchmark_combinatorial('FBA', model)
    benchmark_combinatorial('pFBA', model)
    benchmark_combinatorial('qpFBA', model)
    benchmark_combinatorial('MOMA', model)
    benchmark_combinatorial('lMOMA', model)
    benchmark_combinatorial('ROOM', model)


def benchmark_method(method, model):
    print 'benchmarking', method,
    tstart = time()
    reaction_deletion(model, [], method)
    tend = time()
    print 'took', tend - tstart


def benchmark_methods(modelpath):
    model = load_cbmodel(modelpath, flavor='bigg')
    benchmark_method('FBA', model)
    benchmark_method('pFBA', model)
    benchmark_method('qpFBA', model)
    benchmark_method('MOMA', model)
    benchmark_method('lMOMA', model)
    benchmark_method('ROOM', model)


def benchmark_build_problem(modelpath, n=10):
    model = load_cbmodel(modelpath, flavor='bigg')
    print 'benchmarking build problem for', n, 'instances:',
    tstart = time()
            
    for i in range(n):
        solver = solver_instance()
        solver.build_problem(model)
    tend = time()
    print 'took', tend - tstart


def benchmark_solving_stage(modelpath, n=100):
    model = load_cbmodel(modelpath, flavor='bigg')
    print 'benchmarking solving stage for', n, 'repetitions:',
    solver = solver_instance()
    solver.build_problem(model)
        
    tstart = time()
            
    for i in range(n):
        FBA(model, solver=solver)

    tend = time()
    print 'took', tend - tstart


def main():
#    benchmark_methods_combinatorial(SMALL_TEST_MODEL)
#    benchmark_methods(LARGE_TEST_MODEL)
#    benchmark_build_problem(LARGE_TEST_MODEL)
    benchmark_solving_stage(LARGE_TEST_MODEL)


if __name__ == '__main__':
    main()
