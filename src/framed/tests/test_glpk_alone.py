'''
Created on Oct 25, 2013

@author: Marta Matos
'''
"""solves the following LP problem:
    
   max  x + y - 50
   s.t.: 50x + 24y <= 2400
         30x + 33y <= 2100
         x >= 45
         y >= 5
         
   the solution is:
        obj value: 1.25
    x = 45
    y = 6.25"""
    
from glpk.glpkpi import *


def solve_lp_prob_glpk():
    
    # create problem
    prob = glp_create_prob()
    
    # so that GLPK does not print an output message on solving the problem
    smcp = glp_smcp()
    glp_init_smcp(smcp)
    smcp.msg_lev = GLP_MSG_OFF
    
    # create bounded variables
    glp_add_cols(prob, 2)
    glp_set_col_bnds(prob, 1, GLP_LO, 45, 10e6)
    glp_set_col_bnds(prob, 2, GLP_LO, 6, 10e6)
    
    # create bounded constraints
    glp_add_rows(prob, 2)
    
    coef_ind = intArray(2 + 1)
    coef_val = doubleArray(2 + 1)
    
    coef_ind[1] = 1
    coef_val[1] = 50
    coef_ind[2] = 2
    coef_val[2] = 24
    glp_set_mat_row(prob, 1, 2, coef_ind, coef_val)
    glp_set_row_bnds(prob, 1, GLP_UP, -10e6, 2400)
    
    coef_ind[1] = 1
    coef_val[1] = 30
    coef_ind[2] = 2
    coef_val[2] = 33
    glp_set_mat_row(prob, 2, 2, coef_ind, coef_val)
    glp_set_row_bnds(prob, 2, GLP_UP, -10e6, 2100)
    
    # create objective function
    glp_set_obj_coef(prob, 0, -50)
    glp_set_obj_coef(prob, 1, 1)
    glp_set_obj_coef(prob, 2, 1)
    glp_set_obj_dir(prob, GLP_MAX)
    
    # solve lp problem
    glp_simplex(prob, smcp)
    
    return [glp_get_obj_val(prob), glp_get_col_prim(prob, 1), glp_get_col_prim(prob, 2)]