'''
Created on May 17, 2013

@author: daniel
'''
from framed.io_utils.sbml import load_sbml_model, GPR_CONSTRAINED
from framed.core.fixes import fix_bigg_model
from framed.design.combinatorial import combinatorial_gene_deletion

from cobra.io import read_sbml_model
from cobra.flux_analysis.double_deletion import double_deletion

from time import time

#MODEL = '/Users/daniel/Dropbox/models/ecoli/constraint-based/Reed2003/Ec_iJR904_flux1.xml'
MODEL = 'misc/Ec_iAF1260_flux1.xml'

N = 10

def test_framed_single():
    model = load_sbml_model(MODEL, GPR_CONSTRAINED)
    fix_bigg_model(model)

    tstart = time()
    objective = {}#{'R_EX_succ_e_': 1}
    combinatorial_gene_deletion(model, objective, 2, method='FBA', targets=model.genes.keys())
    tend = time()   
    print 'framed single:', (tend-tstart)
 
        
def test_cobrapy_single():
    model = read_sbml_model(MODEL)
    
    tstart = time()
    double_deletion(model, element_list_1=model.genes, solver='gurobi', number_of_processes=1)
    tend = time()
    
    print 'cobrapy single:', (tend-tstart)
    
def test_cobrapy_parallel():
    model = read_sbml_model(MODEL)

    tstart = time()
    double_deletion(model, element_list_1=model.genes, solver='gurobi', number_of_processes=2)
    tend = time()
    
    print 'cobrapy parallel:', (tend-tstart)


if __name__ == '__main__':
    test_framed_single()
    test_cobrapy_single()
#    test_cobrapy_parallel()