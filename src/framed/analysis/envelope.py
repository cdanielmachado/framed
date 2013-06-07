"""
This module implements classes and methods for analyzing the production envelope

@author: Kai Zhuang, Daniel Machado

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

from variability import FVA
from numpy import linspace


def flux_envelope(model, rxn_x, rxn_y, steps=10):
    """ Calculate the flux envelope for a pair of reactions.

    Arguments:
        model : ConstraintBasedModel -- the model
        rxn_x : str -- reaction on x-axis
        rxn_y : str -- reaction on y-axis
        steps : int -- number of steps to compute (default: 10)

    Returns:
        list (of float), list (of float), list (of float) -- x values, y min values, y max values
    """

    x_range = FVA(model, reactions=[rxn_x])
    xmin, xmax = x_range[rxn_x]
    xvals = linspace(xmin, xmax, steps).tolist()
    ymins, ymaxs = [None]*steps, [None]*steps

    for i, xval in enumerate(xvals):
        constraints = {rxn_x: (xval, xval)}
        y_range = FVA(model, reactions=[rxn_y], constraints=constraints)
        ymins[i], ymaxs[i] = y_range[rxn_y]

    return xvals, ymins, ymaxs


def production_envelope(model, rxn_target, rxn_biomass=None, steps=10):
    """ Calculate the production envelope of the target reaction

    Arguments:
        model : ConstraintBasedModel -- the model
        rxn_target: str -- the target reaction id
        steps: int -- number of steps along the envelope to be calculated (default: 10)
        rxn_biomass: str -- the biomass reaction id (default: None)

    Returns:
        list (of float), list (of float), list (of float) -- biomass values, target minimum values, target maximum values
    """
    if not rxn_biomass:
        rxn_biomass = model.detect_biomass_reaction()

    return flux_envelope(model, rxn_x=rxn_biomass, rxn_y=rxn_target, steps=steps)
