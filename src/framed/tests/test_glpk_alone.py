
from glpk.glpkpi import *


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



def solve_milp_int_prob_glpk():

    # create problem
    prob = glp_create_prob()

    # so that GLPK does not print an output message on solving the problem
    iocp = glp_iocp()
    glp_init_iocp(iocp)
    iocp.msg_lev = GLP_MSG_OFF
    iocp.presolve = GLP_ON

    # create bounded variables
    glp_add_cols(prob, 3)
    glp_set_col_bnds(prob, 1, GLP_FR, -10e6, 10e6)
    glp_set_col_kind(prob, 1, GLP_IV)
    glp_set_col_bnds(prob, 2, GLP_FR, -10e6, 10e6)
    glp_set_col_kind(prob, 2, GLP_IV)
    glp_set_col_bnds(prob, 3, GLP_FR, -10e6, 10e6)
    glp_set_col_kind(prob, 3, GLP_IV)

    # create bounded constraints
    glp_add_rows(prob, 3)

    coef_ind = intArray(3 + 1)
    coef_val = doubleArray(3 + 1)

    coef_ind[1] = 1
    coef_val[1] = 0.8
    coef_ind[2] = 2
    coef_val[2] = 0.2
    coef_ind[3] = 3
    coef_val[3] = 0.3
    glp_set_mat_row(prob, 1, 3, coef_ind, coef_val)
    glp_set_row_bnds(prob, 1, GLP_UP, -10e6, 20)

    coef_ind[1] = 1
    coef_val[1] = 0.4
    coef_ind[2] = 2
    coef_val[2] = 0.3
    coef_ind[3] = 3
    coef_val[3] = 0
    glp_set_mat_row(prob, 2, 3, coef_ind, coef_val)
    glp_set_row_bnds(prob, 2, GLP_UP, -10e6, 10)

    coef_ind[1] = 1
    coef_val[1] = 0.2
    coef_ind[2] = 2
    coef_val[2] = 0
    coef_ind[3] = 3
    coef_val[3] = 0.1
    glp_set_mat_row(prob, 3, 3, coef_ind, coef_val)
    glp_set_row_bnds(prob, 3, GLP_UP, -10e6, 5)

    # create objective function
    glp_set_obj_coef(prob, 0, 0)
    glp_set_obj_coef(prob, 1, 20)
    glp_set_obj_coef(prob, 2, 6)
    glp_set_obj_coef(prob, 3, 8)
    glp_set_obj_dir(prob, GLP_MAX)

    glp_write_lp(prob, None, "bla.lp")
    # solve lp problem
    code = glp_intopt(prob, iocp)
    print code

    return [glp_mip_obj_val(prob), glp_mip_col_val(prob, 1), glp_mip_col_val(prob, 2), glp_mip_col_val(prob, 3)]



"""solves the following MILP problem:

   min  -9x_1 - 5x_2 -6x_3 - 4x_4
   s.t.: 6x_1 + 3x_2 + 5x_3 + 2x_4 <= 9
         0x_1 + 0x_2 + 1x_3 + 1x_4 <= 1
         -1x_1 + 0x_2 + 1x_3 + 0x_4 <= 0
         0x_1 + -1x_2 + 0x_3 + 1x_4 <= 0

         x1, x2, x3, x4 are binary variables

   the solution is:
        obj value: -14
    x_1 = x_2 = 1, x_3 = x_4 = 0"""

def solve_milp_bin_prob_glpk():

    # create problem
    prob = glp_create_prob()

    # so that GLPK does not print an output message on solving the problem
    iocp = glp_iocp()
    glp_init_iocp(iocp)
    iocp.msg_lev = GLP_MSG_OFF
    iocp.presolve = GLP_ON

    # create bounded variables
    glp_add_cols(prob, 4)
    glp_set_col_bnds(prob, 1, GLP_DB, 0, 1)
    glp_set_col_kind(prob, 1, GLP_BV)
    glp_set_col_bnds(prob, 2, GLP_DB, 0, 1)
    glp_set_col_kind(prob, 2, GLP_BV)
    glp_set_col_bnds(prob, 3, GLP_DB, 0, 1)
    glp_set_col_kind(prob, 3, GLP_BV)
    glp_set_col_bnds(prob, 4, GLP_DB, 0, 1)
    glp_set_col_kind(prob, 4, GLP_BV)


    # create bounded constraints
    glp_add_rows(prob, 4)

    coef_ind = intArray(4 + 1)
    coef_val = doubleArray(4 + 1)

    coef_ind[1] = 1
    coef_val[1] = 6
    coef_ind[2] = 2
    coef_val[2] = 3
    coef_ind[3] = 3
    coef_val[3] = 5
    coef_ind[4] = 4
    coef_val[4] = 2
    glp_set_mat_row(prob, 1, 4, coef_ind, coef_val)
    glp_set_row_bnds(prob, 1, GLP_UP, -10e6, 9)

    coef_ind[1] = 1
    coef_val[1] = 0
    coef_ind[2] = 2
    coef_val[2] = 0
    coef_ind[3] = 3
    coef_val[3] = 1
    coef_ind[4] = 4
    coef_val[4] = 1
    glp_set_mat_row(prob, 2, 4, coef_ind, coef_val)
    glp_set_row_bnds(prob, 2, GLP_UP, -10e6, 1)

    coef_ind[1] = 1
    coef_val[1] = -1
    coef_ind[2] = 2
    coef_val[2] = 0
    coef_ind[3] = 3
    coef_val[3] = 1
    coef_ind[4] = 4
    coef_val[4] = 0
    glp_set_mat_row(prob, 3, 4, coef_ind, coef_val)
    glp_set_row_bnds(prob, 3, GLP_UP, -10e6, 0)

    coef_ind[1] = 1
    coef_val[1] = 0
    coef_ind[2] = 2
    coef_val[2] = -1
    coef_ind[3] = 3
    coef_val[3] = 0
    coef_ind[4] = 4
    coef_val[4] = 1
    glp_set_mat_row(prob, 4, 4, coef_ind, coef_val)
    glp_set_row_bnds(prob, 4, GLP_UP, -10e6, 0)
    # create objective function
    glp_set_obj_coef(prob, 0, 0)
    glp_set_obj_coef(prob, 1, -9)
    glp_set_obj_coef(prob, 2, -5)
    glp_set_obj_coef(prob, 3, -6)
    glp_set_obj_coef(prob, 3, -4)
    glp_set_obj_dir(prob, GLP_MIN)

    glp_write_lp(prob, None, "bla.lp")
    # solve lp problem
    glp_intopt(prob, iocp)
    
    return [glp_mip_obj_val(prob), glp_mip_col_val(prob, 1), glp_mip_col_val(prob, 2), glp_mip_col_val(prob, 3), glp_mip_col_val(prob, 4)]
