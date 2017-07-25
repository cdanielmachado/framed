""" This module implements some basic plotting utilities for cobra analysis.

Author: Daniel Machado
   
"""

from variability import flux_envelope
from simulation import FBA
import matplotlib.pyplot as plt
import numpy as np


def plot_flux_envelope(model, r_x, r_y, substrate=None, constraints=None, reference=None, alternatives=None,
                       label_x=None, label_y=None, filename=None, steps=10, plot_kwargs=None, fill_kwargs=None):
    """ Plots the flux envelope for a pair of reactions.
    
    Arguments:
        model (CBModel): the model
        r_x (str): reaction on x-axis
        r_y (str): reaction on y-axis
        substrate (str): compute yields instead of rates (optional)
        constraints (dict): additional constraints
        reference (dict): display location of reference flux distribution (eg: wild-type strain) (optional)
        alternatives (list): display location of alternative flux distributions (eg: mutant strains) (optional)
        label_x (str): change x label (optional, uses reaction name by default)
        label_y (str): change y label (optional, uses reaction name by default)
        filename (str): filename to save image (optional), otherwise display on screen (default)
        steps (int): number of steps to compute (default: 10)
        plot_kwargs (dict): keyword parameters to pass along to *pyplot.plot* (optional)
        fill_kwargs (dict): keyword parameters to pass along to *pyplot.fill_between* (optional)

    Returns:
        matplotlib.Axes: axes object
    """

    _, ax = plt.subplots()

    if not plot_kwargs:
        plot_kwargs = {'color': 'b'}

    if not fill_kwargs:
        fill_kwargs = {'color': 'b', 'alpha': 0.2}

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

    ax.plot(xvals, ymins, **plot_kwargs)
    ax.plot(xvals, ymaxs, **plot_kwargs)
    ax.plot([xvals[0], xvals[0]], [ymins[0], ymaxs[0]], **plot_kwargs)
    ax.plot([xvals[-1], xvals[-1]], [ymins[-1], ymaxs[-1]], **plot_kwargs)

    ax.fill_between(xvals, ymins, ymaxs, **fill_kwargs)

    if alternatives:
        for fluxes in alternatives:
            ax.plot(fluxes[r_x], fluxes[r_y], 'ro', markersize=5)

    if reference:
        ax.plot(reference[r_x], reference[r_y], 'bo', markersize=5)

    ax.set_xlabel(label_x) if label_x else ax.set_xlabel(model.reactions[r_x].name)
    ax.set_ylabel(label_y) if label_y else ax.set_ylabel(model.reactions[r_y].name)

    xmin, xmax = min(xvals), max(xvals)
    dx = 0.03 * (xmax - xmin)
    ax.set_xlim((xmin - dx, xmax + dx))

    ymin, ymax = min(ymins), max(ymaxs)
    dy = 0.03 * (ymax - ymin)
    ax.set_ylim((ymin - dy, ymax + dy))

    if filename:
        plt.savefig(filename)
    else:
        return ax


def _normalize_list(values, x):
    for i in range(len(values)):
        values[i] = values[i] / x if values[i] is not None else None


def _normalize_dict(fluxes, x):
    for r_id, flux in fluxes.items():
        fluxes[r_id] = flux / x if flux is not None else None


def plot_flux_bounds(range1, range2=None, keys=None, log=False, unbounded=1000):
    """ Plot and compare flux ranges (although other types of ranges also supported).

    Args:
        range1 (dict): flux ranges
        range2 (dict): alternative flux ranges (optional)
        keys (list): only display reactions in this list (optional)
        log (bool): log scale (default: False)
        unbounded (float): threshold to display unbounded values (default: 1000)

    """

    if not keys:
        keys = range1.keys()

    def bounded_left(x):
        return -unbounded if x is None else x

    def bounded_right(x):
        return unbounded if x is None else x

    lb1 = np.array([bounded_left(range1[key][0]) for key in keys])
    ub1 = np.array([bounded_right(range1[key][1]) for key in keys])

    if log:
        lb1 = np.log10(lb1)
        ub1 = np.log10(ub1)

    if not range2:
        idx = np.arange(len(keys))
        plt.barh(idx, left=lb1, width=(ub1 - lb1))
        plt.yticks(idx + 0.5, keys)
        plt.ylim(0, len(keys))

    else:
        lb2 = np.array([bounded_left(range2[key][0]) for key in keys])
        ub2 = np.array([bounded_right(range2[key][1]) for key in keys])
        idx1 = np.arange(1, 2 * len(keys), 2)
        idx2 = np.arange(0, 2 * len(keys), 2)

        if log:
            lb2 = np.log10(lb2)
            ub2 = np.log10(ub2)

        plt.barh(idx1, left=lb1, width=(ub1 - lb1), color='#7ba6ed')
        plt.barh(idx2, left=lb2, width=(ub2 - lb2), color='#43c6c2')
        plt.yticks(idx1, keys)
        plt.ylim(0, 2 * len(keys))
