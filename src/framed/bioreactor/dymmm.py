__author__ = 'kaizhuang'

"""
This module implements the dynamic multi-species metabolic modeling (DyMMM) framework

TODO: rewrite the _ode_RHS method in Bioreactor
"""

from copy import deepcopy
import numpy
from ..solvers import solver_instance
from ..solvers.solver import Status
from ..analysis.simulation import FBA


class Organism(object):
    """
    Organism describes a generic biological organism.

    :param model: a mathematical model of the organism
    :param environment: a reference to which environment the organism is placed in
    """

    def __init__(self, model, fba_objective=None):
        self.model = deepcopy(model)
        self.fba_constraints = {}
        self.environment = None  # upon initiation, the organism is not placed in any environment

        if fba_objective:
            self.fba_objective = fba_objective
        else:
            self.fba_objective = {model.detect_biomass_reaction(): 1}

    def update(self):
        """
        This method updates the states of the organism.
        the organism's response to changes in the environmental conditions should be described here

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

    def update(self):
        """
        updates the states of the environment.
        these state changes can be triggered by either natural events or human interventions

        ** this is an abstract method, must be implemented for specific environments **
        """
        raise NotImplementedError("update method for individual organisms must be implemented for DyMMM to work")

    def ode_RHS(self, y, t):
        """
        this is the Right Hand Side of the system of ODE that describe the dynamic multi-species system
        :param y: state variables such as liquid volume, biomass concentrations, and metabolite concentrations
        :param t: time
        :return:

        ** this is an abstract method, must be implemented for specific environments **
        """
        raise NotImplementedError

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

#    def add_metabolite(self, metabolite):
#        assert metabolite[:3] == 'EX_' and metabolite[-3:], 'A reaction name (eg. \'EX_glc(e)\') is expected '
#        self.metabolites.append(metabolite)

#    def add_metabolites(self, metabolites):
#        if isinstance(metabolites, str):
#            self.add_metabolite(metabolites)
#        elif isinstance(metabolites, list):
#            for metabolite in metabolites:
#                self.add_metabolite(metabolite)
#        else:
#            raise TypeError('A reaction name (eg. \'EX_glc(e)\') or a list of reaction names is expected')


class Bioreactor(Environment):
    """
    This class describes a generic bioreactor with one influent (feed) stream and one effluent stream

    :param organisms: a list of Organism objects
    :param metabolites: a list of string objects containing exchange reactions names.  eg. 'EX_glc(e)'
    Xfeed: concentration of organisms in the feed stream
    Sfeed: concentration of metabolites in the feed stream
    """
    def __init__(self, organisms, metabolites):
        self.set_organisms(organisms)
        self.set_metabolites(metabolites)
        self.Xfeed = numpy.zeros(len(organisms))
        self.Sfeed = numpy.zeros(len(metabolites))
        self.flow_rate_in = 0
        self.flow_rate_out = 0
        self.volume_max = None
        self.initial_conditions = []

    def set_Xfeed(self, Xfeed):
        assert len(Xfeed) == len(self.organisms), 'The length of Xfeed should equal to the number of organisms'
        self.Xfeed = Xfeed

    def set_Sfeed(self, Sfeed):
        assert len(Sfeed) == len(self.metabolites),  'The length of Sfeed should equal to the number of metabolites'
        self.Sfeed = Sfeed

    def set_initial_conditions(self, Vinit, Xinit, Sinit):
        assert type(Xinit) == type(Sinit) == list
        if type(Vinit) != list:
            Vinit = [Vinit]
        assert len(Vinit) == 1
        assert len(Xinit) == len(self.organisms), 'The length of Xinit should equal to the number of organisms'
        assert len(Sinit) == len(self.metabolites), 'The length of Sinit should equal to the number of organisms'
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
        dy[1:number_of_organisms + 1] = mu * self.X + self.flow_rate_in / self.V * (self.Xfeed - self.X)   # dX/dt [g/L/hr]
        dy[number_of_organisms + 1:] = numpy.dot(self.X, vs) + self.flow_rate_in / self.V * (self.Sfeed - self.S)    # dS/dt [mmol/L/hr]

        return dy



