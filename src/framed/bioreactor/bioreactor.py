""" This module defines the common classes used for modeling and analyzing bioreactors

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
__author__ = 'kaizhuang'


from copy import deepcopy
from ..solvers import solver_instance
from ..solvers.solver import Status
from scipy.integrate import ode
import numpy
import collections


class Organism(object):
    """
    Organism describes a generic biological organism.
    """

    def __init__(self, model, fba_objective=None, fba_constraints={}, model_deepcopy=True):
        """
        :param model: the mathematical model of the organism
        :param fba_objective: dict -- the FBA objective function.  (only useful if model is a FBA model)
        :param fba_constraints: dict -- none standard FBA constraints.  This can be useful for creating knockout strains
        :param model_deepcopy: bool -- if True, a deepcopy of the model will be created inside the Organism instance,
                                        otherwise, a reference of the model will be created
        :return: none
        """
        if model_deepcopy:
            self.model = deepcopy(model)
        else:
            self.model = model
        self.fba_constraints = fba_constraints
        self.environment = None  # upon initiation, the organism is not placed in any environment

        if fba_objective:
            self.fba_objective = fba_objective
        else:
            self.fba_objective = {model.detect_biomass_reaction(): 1}

    def update(self):
        """
        The update() method is used to change the internal state of the organism.
        - This method is called at each integration step.
        - One usage of this method is to describe how the FBA uptake constraint changes in response to the changes in
            the metabolite concentrations.

        ** this is an abstract method, must be implemented in strain specific subclasses **
         """

        raise NotImplementedError


class Environment(object):
    """
    This class describes a generic environment that contains a number of organisms and metabolites
    """

    def __init__(self):
        self.organisms = []
        self.metabolites = []
        self.initial_conditions = []

    def update(self):
        """
        The update() method is used to change the internal state of the environment.
        This method is called at each integration step.

        ** this is an abstract method, must be implemented for specific environments **
        """
        raise NotImplementedError("update() method must be implemented for the specific environment")

    def set_organisms(self, organisms):
        self.organisms = []
        self.add_organisms(organisms)

    def set_metabolites(self, metabolites):
        self.metabolites = []
        self.add_metabolites(metabolites)

    def add_organism(self, organism):
        organism.environment = self
        self.organisms.append(organism)

    def add_organisms(self, organisms):
        for organism in organisms:
            self.add_organism(organism)

    def add_metabolite(self, metabolite):
        self.metabolites.append(metabolite)

    def add_metabolites(self, metabolites):
        for metabolite in metabolites:
            self.add_metabolite(metabolite)


class DynamicSystem(object):
    """
    This class describes a generic dynamic system
    """
    def ode_RHS(self, y, t):
        """
        this is the Right Hand Side of the system of ODE that describe the dynamic multi-species system
        :param y: state variables such as liquid volume, biomass concentrations, and metabolite concentrations
        :param t: time
        :return:

        ** this is an abstract method, must be implemented for specific environments **
        """
        raise NotImplementedError("the RHS of the ODE must be described for the each specific environment")

    def integrate(self, t0, tf, dt, initial_conditions, solver, verbose=False):
        """
        the integrate() solves the ODE of the dynamic system using the designated solver
        :param t0: initial time
        :param tf: final time
        :param dt: time step
        :param initial_conditions: array-like -- initial conditions of the ODE system
        :param solver: the designated solver
        :return:
        """
        if solver == 'analytical':
            try:
                t, y = self.analytical_integrator(t0, tf, dt, initial_conditions, solver, verbose)
            except NotImplementedError:
                print 'analytical solver have no been implemented yet. will use numerical solver dopri5.'
                t, y = self.numerical_integrator(t0, tf, dt, initial_conditions, solver='dopri5')
        else:
            t, y = self.numerical_integrator(t0, tf, dt, initial_conditions, solver, verbose)
        return t, y

    def numerical_integrator(self, t0, tf, dt, initial_conditions, solver, verbose):
        """
        the numerical_integrator() method integrates the ODE of the dynamic system using a numerical solver
        """
        f = self._ode_RHS
        if initial_conditions:
            y0 = initial_conditions
        else:
            y0 = self.initial_conditions

        MdFBA_ode = ode(f).set_integrator(solver)
        MdFBA_ode.set_initial_value(y0, t0)

        t = [t0]
        y = [y0]

        while MdFBA_ode.successful() and MdFBA_ode.t < tf:
            MdFBA_ode.integrate(MdFBA_ode.t + dt)

            t.append(MdFBA_ode.t)
            y.append(MdFBA_ode.y)

            if verbose:
                print MdFBA_ode.t

        t = numpy.array(t)
        y = numpy.array(y)
        return t, y

    def analytical_integrator(self, t0, tf, dt, initial_conditions, solver, verbose):
        """
        the analytical_integrator() method integrates the ODE of the dynamic system using a user-defined analytical method

        ** this is an abstract method, must be implemented for specific dynamic systems **
        """
        raise NotImplementedError


class Bioreactor(Environment, DynamicSystem):
    """
    This class describes a generic bioreactor with one influent (feed) stream and one effluent stream

    :param organisms: a list of Organism objects
    :param metabolites: a list of string objects containing exchange reactions names.  eg. 'EX_glc(e)'
    Xfeed: concentration of organisms in the feed stream [g/L]
    Sfeed: concentration of metabolites in the feed stream [mmol/L]
    deltaX: custom defined terms to dX/dt [g/L/hr]
    deltaS: custom defined terms to dS/dt [mmol/L/hr]

    """
    def __init__(self, organisms, metabolites, flow_rate_in=0, flow_rate_out=0, volume_max=None, Xfeed=None,
                 Sfeed=None, deltaX=None, deltaS=None, initial_conditions=[]):
        if not isinstance(organisms, collections.Iterable):
            organisms = [organisms]
        if not isinstance(metabolites, collections.Iterable):
            metabolites = [metabolites]
        self.set_organisms(organisms)
        self.set_metabolites(metabolites)

        self.flow_rate_in = flow_rate_in
        self.flow_rate_out = flow_rate_out
        self.volume_max = volume_max

        if Xfeed:
            self.set_Xfeed(Xfeed)
        else:
            self.Xfeed = numpy.zeros(len(organisms))

        if Sfeed:
            self.set_Sfeed(Sfeed)
        else:
            self.Sfeed = numpy.zeros(len(metabolites))

        if deltaX:
            self.deltaX = deltaX
        else:
            self.deltaX = numpy.zeros(len(organisms))

        if deltaS:
            self.deltaS = deltaS
        else:
            self.deltaS = numpy.zeros(len(metabolites))

        self.initial_conditions = initial_conditions

    def set_Xfeed(self, Xfeed):
        assert len(Xfeed) == len(self.organisms), 'The length of Xfeed should equal to the number of organisms'
        self.Xfeed = Xfeed

    def set_Sfeed(self, Sfeed):
        assert len(Sfeed) == len(self.metabolites),  'The length of Sfeed should equal to the number of metabolites'
        self.Sfeed = Sfeed

    def set_initial_conditions(self, Vinit, Xinit, Sinit):
        assert type(Vinit) == type(Xinit) == type(Sinit) == list
        assert len(Vinit) == 1
        assert len(Xinit) == len(self.organisms), 'The length of Xinit should equal to the number of organisms'
        assert len(Sinit) == len(self.metabolites), 'The length of Sinit should equal to the number of metabolites'
        self.initial_conditions = Vinit + Xinit + Sinit

    def update(self):
        pass

    def _ode_RHS(self, t, y):
        """
        the RHS of the ODE that describe the bioreactor system
        :param y:
            y[0]: volume
            y[1] to y[number_of_organisms]: biomass of the organisms
            y[number_of_organisms + 1] to y[-1] concentration of metabolites

        :param t: time
        :return: dy
        """
        number_of_organisms = len(self.organisms)
        number_of_metabolites = len(self.metabolites)
        assert(len(y) == 1 + number_of_organisms + number_of_metabolites)
        dy = numpy.zeros(len(y))

        # creating user-friendly aliases for y
        self.V = y[0]
        self.X = y[1:number_of_organisms + 1]
        self.S = y[number_of_organisms + 1:]

        # assigning growth rates and metabolic production/consumption rates here, the rates are calculated using FBA
        # through CobraPy's optimize method.  These rates can be calculated and assigned using other methods.
        # For example, a Michaelis-Menten kinetic expression can be used

        vs = numpy.zeros([number_of_organisms, number_of_metabolites])     # fluxes through metabolites
        mu = numpy.zeros([number_of_organisms])                          # growth rates of organisms

        for i, organism in enumerate(self.organisms):
            organism.update()       # updates the constraints of the organism model based on environment conditions

            if t == 0:
                organism.solver = solver_instance()
                organism.solver.build_problem(organism.model)

            solution = organism.solver.solve_lp(organism.fba_objective, constraints=organism.fba_constraints)

            #solution = FBA(organism.model, solver=organism.solver, constraints=organism.fba_constraints)

            if solution.status == Status.OPTIMAL:
                mu[i] = solution.fobj

                for j, metabolite in enumerate(self.metabolites):
                    if metabolite in organism.model.reactions.keys():
                        vs[i, j] = solution.values[metabolite]
            else:
                mu[i] = 0
                #print 'no growth'

        self.update()   # updating bioreactor parameters such as flow rates and feed concentration

        # calculating the rates of change of reactor volume[L], biomass [g/L] and metabolite [mmol/L]
        dy[0] = self.flow_rate_in - self.flow_rate_out      # dV/dt [L/hr]
        dy[1:number_of_organisms + 1] = mu * self.X + self.flow_rate_in / self.V * (self.Xfeed - self.X) + self.deltaX  # dX/dt [g/L/hr]
        dy[number_of_organisms + 1:] = numpy.dot(self.X, vs) + self.flow_rate_in / self.V * (self.Sfeed - self.S) +self.deltaS   # dS/dt [mmol/L/hr]

        return dy
