""" This module implements some basic plotting utilities for cobra analysis.

Author: Daniel Machado
   
"""

from variability import flux_envelope
from simulation import FBA
import matplotlib.pyplot as plt


def plot_flux_envelope(model, r_x, r_y, substrate=None, constraints=None, reference=None, alternatives=None,
                       label_x=None, label_y=None, filename=None, steps=10):
    """ Plots the flux envelope for a pair of reactions.
    
    Arguments:
        model (CBModel): the model
        r_x (str): reaction on x-axis
        r_y (str): reaction on y-axis
        substrate (str): compute yields instead of rates (optional)
        constraints (dict): additional constraints
        filename (str): filename to save image (optional), otherwise display on screen (default)
        steps (int): number of steps to compute (default: 10)
        
    Returns:
        tuple: x values, y min values, y max values
    """

    xvals, ymins, ymaxs = flux_envelope(model, r_x, r_y, steps, constraints)

    if substrate:
        sol = FBA(model)
        uptk = abs(sol.values[substrate])
        _normalize_list(xvals, uptk)
        _normalize_list(ymins, uptk)
        _normalize_list(ymaxs, uptk)

        if reference:
            _normalize_dict(reference, abs(reference[substrate]))

        if alternatives:
            for fluxes in alternatives:
                _normalize_dict(fluxes, abs(fluxes[substrate]))

    plt.plot(xvals, ymins, 'b', xvals, ymaxs, 'b',
             [xvals[0], xvals[0]], [ymins[0], ymaxs[0]], 'b',
             [xvals[-1], xvals[-1]], [ymins[-1], ymaxs[-1]], 'b')

    if alternatives:
        for fluxes in alternatives:
            plt.plot(fluxes[r_x], fluxes[r_y], 'ro', markersize=5)

    if reference:
        plt.plot(reference[r_x], reference[r_y], 'bo', markersize=5)

    plt.xlabel(label_x) if label_x else plt.xlabel(model.reactions[r_x].name)
    plt.ylabel(label_y) if label_y else plt.ylabel(model.reactions[r_y].name)

    xmin, xmax = min(xvals), max(xvals)
    dx = 0.03 * (xmax - xmin)
    plt.xlim((xmin - dx, xmax + dx))

    ymin, ymax = min(ymins), max(ymaxs)
    dy = 0.03 * (ymax - ymin)
    plt.ylim((ymin - dy, ymax + dy))

    if filename:
        plt.savefig(filename)


def _normalize_list(values, x):
    for i in range(len(values)):
        values[i] = values[i] / x if values[i] is not None else None


def _normalize_dict(fluxes, x):
    for r_id, flux in fluxes.items():
        fluxes[r_id] = flux / x if flux is not None else None 
    