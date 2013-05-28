''' This module implements some basic plotting utilities for common methods.

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

from variability import flux_cone_projection
from matplotlib.pyplot import plot, xlabel, ylabel, show, savefig

def plot_flux_cone_projection(model, r_x, r_y, filename=None, steps=10):
    """ Plots the flux cone projection for a pair of reactions.
    
    Arguments:
        model : ConstraintBasedModel -- the model
        r_x : str -- reaction on x-axis
        r_y : str -- reaction on y-axis
        filename : str -- filename to save image (optional), otherwise display on screen (default)
        steps : int -- number of steps to compute (default: 10)
        
    Returns:
        list (of float), list (of float), list (of float) -- x values, y min values, y max values
    """

    xvals, ymins, ymaxs = flux_cone_projection(model, r_x, r_y, steps)
    plot(xvals, ymins, 'k', xvals, ymaxs, 'k')
    xlabel(model.reactions[r_x].name)
    ylabel(model.reactions[r_y].name)
    
    if filename:
        savefig(filename)
    else:
        show()