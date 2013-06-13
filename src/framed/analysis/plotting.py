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
from matplotlib.pyplot import plot, xlabel, ylabel, show, savefig, xlim, ylim

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
    plot(xvals, ymins, 'b', xvals, ymaxs, 'b',
         [xvals[0], xvals[0]], [ymins[0], ymaxs[0]], 'b',
         [xvals[-1], xvals[-1]], [ymins[-1], ymaxs[-1]], 'b')
    xlabel(model.reactions[r_x].name)
    ylabel(model.reactions[r_y].name)
    xmin, xmax = min(xvals), max(xvals)
    dx = 0.1 * (xmax - xmin)
    xlim((xmin - dx, xmax + dx))
    ymin, ymax = min(ymins), max(ymaxs)
    dy = 0.1 * (ymax - ymin)
    ylim((ymin - dy, ymax + dy))
    
    if filename:
        savefig(filename)
    else:
        show()