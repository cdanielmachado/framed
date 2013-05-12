__author__ = 'kaizhuang'

"""
This module defines the classes and methods used by m-dFBA (multi-species dynamic flux balance analysis)
"""

from copy import deepcopy


class Organism(object):
    """
    Organism describes a generic biological organism.

    :param model: a mathematical model of the organism
    :param environment: a reference to which environment the organism is placed in
    """

    def __init__(self, model):
        self.model = deepcopy(model)
        self.environment = None  # upon initiation, the organism is not placed in any environment

    def update(self, update_function=None):
        """
        This method updates the states of the organism.
        the organism's response to changes in the environmental conditions should be described here

        ** this is an abstract method, must be implemented in strain specific subclasses **
         """

        raise NotImplementedError


