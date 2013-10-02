""" This module implements classes and methods useful for the Dynamic Strain Scanning Optimization strategy

Cite:
Kai Zhuang et al. 2013. Dynamic strain scanning optimization: an efficient strain design strategy for balanced yield,
titer, and productivity.

@author: Kai Zhuang

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
"""

from ..solvers import solver_instance
from ..analysis.dfba import dFBA
from ..analysis.variability import production_envelope
from base import *


def make_envelope_strains(base_organism, r_substrate, r_target, N=10, constraints=None):
    """
    Create N strains along the product envelope.
        (Used for Steps 1 and 2 of DySScO strategy)

    Arguments:
        base_organism: Organism -- the host organism used to product the target product
        r_substrate: str -- the rxn id of the substrate
        r_target: str -- the rxn id of the target product
        N: int -- the number of strains to be generated along the production envelope
        constraints: dict -- custom constraints

    Returns:
        strains: list of Organism -- N strains that are fixed to the production envelope
    """
    strains = []

    # add custom constraints to base_organism
    if constraints:
        base_organism.fba_constraints.update(constraints)

    # create the product envelope
    xvals, ymins, ymaxs = production_envelope(base_organism.model, r_target, steps=N, constraints=constraints)

    # finding the maximum r_substrate uptake rate
    if r_substrate in base_organism.fba_constraints:
        vSmax = base_organism.fba_constraints[r_substrate][0]
    else:
        vSmax = base_organism.model.bounds[r_substrate][0]

    # create new strains along the production envelope
    for i, mu in enumerate(xvals):
        # creating a new strain
        strain = deepcopy(base_organism)                        # create a deepcopy of the base_organism
        strain.fba_constraints[r_target] = (ymaxs[i], ymaxs[i])   # fix target production at ymax[i]
        #strain.Y = float(-ymaxs[i]/vSmax)                       # store the yield of the strain
        #strain.mu = mu                                          # store the growth rate of the strain
        strain.id = base_organism.id + '_mu_' + str(round(mu, 3))
        strains.append(strain)

    return strains


def calculate_performances(strains, bioreactor, r_substrate, r_target, t0, tf, dt, initial_conditions=[],
<<<<<<< HEAD
                          dfba_solver='dopri5', additional_yields=[], verbose=False, get_dfba_solution=False,
                          func_dfba2yield=None):
=======
                          dfba_solver='dopri5', additional_yields=[], verbose=False, save_dfba_solution=False):
>>>>>>> master
    """
    calculates the performances of a list of strains in a given bioreactor

    Arguments:
        strains: list (of Organism)
        bioreactor: Bioreactor
        r_substrate: str -- reaction id of the substrate uptake reaction
        r_target: str -- reaction id of the target metabolite exchange reaction
        t0: float -- initial time
        tf: float -- final time
        dt: float -- time step
        initial_conditions: list (of float) -- the initial conditions in the order of V0, X0, S0 (default: None)
        dfba_solver: str -- ODE solver.  (default: 'dopri5')
        additional_yields: list (of str) -- the reaction ids of the additional yields (yields other than target yield)
                                            to be calculated.
        verbose: bool -- Verbosity control.  (default: False).
        get_dfba_solution: bool -- controls whether dfba solutions are returned (default: False)
        func_dfba2yield: function -- if None, yield is calculated from FBA solutions.
                                     if a function is passed in, yield is calculated from dFBA solutions.
    Returns:
        performances: list (of Dict) -- a list of dictionaries, each entry contains the calculated performance metrics
        of a strain
<<<<<<< HEAD
=======

    Conditional Returns:
        dfba_solutions: list (of Dict) -- a list of dictionaries containing the dFBA solutions
                                         (this is returned only if save_dfba_solution is set to True)
>>>>>>> master
    """

    performances = []

    for strain in strains:
        performance = calculate_performance(strain, bioreactor, r_substrate, r_target, t0, tf, dt,
<<<<<<< HEAD
                                                initial_conditions, dfba_solver, additional_yields, verbose,
                                                get_dfba_solution, func_dfba2yield)
=======
                                                               initial_conditions, dfba_solver, additional_yields,
                                                               verbose, save_dfba_solution)
>>>>>>> master
        performances.append(performance)

    return performances


def calculate_performance(strain, bioreactor, r_substrate, r_target, t0, tf, dt, initial_conditions=[],
<<<<<<< HEAD
                          dfba_solver='vode', additional_yields=[], verbose=False, get_dfba_solution=False,
                          func_dfba2yield=None):
=======
                          dfba_solver='dopri5', additional_yields=[], verbose=False, save_dfba_solution=False):
>>>>>>> master
    """
    calculates the performances of a list of strains in a given bioreactor

    Arguments:
        strain: Organism
        bioreactor: Bioreactor
        r_substrate: str -- reaction id of the substrate uptake reaction
        r_target: str -- reaction id of the target metabolite exchange reaction
        t0: float -- initial time
        tf: float -- final time
        dt: float -- time step
        initial_conditions: list (of float) -- the initial conditions in the order of V0, X0, S0 (default: None)
        dfba_solver: str -- ODE solver.  (default: 'dopri5')
        additional_yields: list (of str) -- the reaction ids of the additional yields (yields other than target yield)
                                            to be calculated.
        verbose: bool -- Verbosity control.  (default: False).
        get_dfba_solution: bool -- controls whether dfba solutions are returned (default: False)
        func_dfba2yield: function -- a function used to calculate yield from dfba solutions

    Returns:
        performance: Dict -- contains the calculated performance metrics of a strain
<<<<<<< HEAD
    """
    performance = {'strain_id': strain.id}
    r_biomass = strain.model.detect_biomass_reaction()

=======

    Conditional Returns:
        dfba_solution: Dict -- contains the dFBA solutions
                            (this is returned only if save_dfba_solution is set to True)
    """
    performance = {'strain': strain.id}
    r_biomass = strain.model.detect_biomass_reaction()

    #print strain.fba_constraints
>>>>>>> master
    # perform FBA simulation
    if verbose:
        print 'Performing FBA simulation.'
    if hasattr(strain, 'solver'):
        fba_solution = strain.solver.solve_lp(strain.fba_objective, constraints=strain.fba_constraints)
    else:
        strain.solver = solver_instance()
        strain.solver.build_problem(strain.model)
        fba_solution = strain.solver.solve_lp(strain.fba_objective, constraints=strain.fba_constraints)

    # growth, substrate uptake, and target production rates from FBA solution
    v_biomass = fba_solution.values[r_biomass]
    v_target = fba_solution.values[r_target]
    v_substrate = fba_solution.values[r_substrate]

<<<<<<< HEAD
    # if the strain does not grow, set growth, titer, productivity to zero, and calculate yield from FBA
    if v_biomass <= 10**-6:
        performance['growth'] = 0
        performance['yield_biomass'] = 0
        performance['titer'] = 0
        performance['productivity'] = 0
        performance['yield'] = - v_target/v_substrate


        for r_id in additional_yields:
            id = 'yield_' + r_id.lstrip('R_EX_').rstrip('_e')
            performance[id] = - fba_solution.values[r_id] / v_substrate
        dfba_solution = None

        if verbose:
            print 'none growing'

    # if the strain grows but does not produce, set yield, titer, productivity to zero, and calculate growth from FBA
    elif v_target <= 10**-6:
        performance['growth'] = v_biomass
        performance['yield_biomass'] = - v_biomass/v_substrate
        performance['titer'] = 0
        performance['productivity'] = 0
        performance['yield'] = 0

=======
    # if the strain does not grow or does not produce target product, set yield, titer, productivity to zero
    if (v_biomass <= 0) or (v_target <= 0):
        performance['growth_rate'] = v_biomass
        performance['product_titer'] = 0
        performance['productivity'] = 0
        performance['product_yield'] = abs(v_target / v_substrate)
        performance['biomass_yield'] = abs(v_biomass / v_substrate)
>>>>>>> master
        for r_id in additional_yields:
            id = 'yield_' + r_id.lstrip('R_EX_').rstrip('_e')
            performance[id] = 0
        dfba_solution = None

<<<<<<< HEAD
        if verbose:
            print 'none producing'

=======
>>>>>>> master
    # if the strain both grows and produces, perform dFBA simulation, and calculate yield, titer, productivity
    else:
        # perform dFBA simulation
        bioreactor.set_organisms([strain])

        if verbose:
            print 'Performing dFBA simulation.'
        dfba_solution = dFBA(bioreactor, t0, tf, dt, initial_conditions, solver=dfba_solver, verbose=verbose)

        # calculating titer and productivity from dFBA solution
        T = dfba_solution[r_target].max()
        index = dfba_solution[r_target].argmax()        # the index at which the production is finished
        P = T/dfba_solution['time'][index]

        # calculate yield from dFBA solution if the method is known, otherwise calculate yield using FBA
<<<<<<< HEAD
        if func_dfba2yield is None:
            Y = - v_target / v_substrate
        else:
            Y = func_dfba2yield(dfba_solution)

        performance['growth'] = v_biomass
        performance['yield_biomass'] = - v_biomass/v_substrate
        performance['titer'] = T
=======
        if hasattr(bioreactor, 'calculate_yield_from_dfba'):
            Y = bioreactor.calculate_yield_from_dfba(dfba_solution, r_substrate, r_target)
        else:
            Y = abs(v_target / v_substrate)


        performance['growth_rate'] = v_biomass
        performance['product_titer'] = T
>>>>>>> master
        performance['productivity'] = P
        performance['product_yield'] = Y
        performance['biomass_yield'] = abs(v_biomass / v_substrate)

        # calculate additional yields
        for r_id in additional_yields:
            id = 'yield_' + r_id.lstrip('R_EX_').rstrip('_e')
            performance[id] = abs(fba_solution.values[r_id] / v_substrate)

<<<<<<< HEAD
    if get_dfba_solution:
        performance['dfba_solution'] = dfba_solution

    return performance


def performances2metrics(performances):
    """
    get a dictionary of metrics from a list of performances
    Arguments:
        performances: Dict (of Dict) -- each list entry contains a dictionary containing the performance metrics
                                        of a strain

    Returns:
        metrics: Dict (of list) -- the performance metrics in form of a dictionary of lists
    """
    performance_metrics = {}
    metrics = performances[0].keys()
=======
    if save_dfba_solution:
        performance['dfba_solution'] = dfba_solution
    else:
        performance['dfba_solution'] = None
>>>>>>> master

    return performance

