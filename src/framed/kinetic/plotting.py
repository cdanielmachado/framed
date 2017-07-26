"""
This module implements plotting utilities for kinetic models.

Author: Daniel Machado
"""

from ..kinetic.simulation import time_course
from matplotlib.pyplot import figure, subplot2grid
from numpy import array


def plot_timecourse(model, time, steps=100, parameters=None, metabolites=None, xlabel=None, ylabel=None,
                    data=None, data_steps=None):
    """ Plot a time-course simulation using a kinetic model.

    Args:
        model (ODEModel): kinetic model
        time (float): final simulation time
        steps (int): number of simulations steps (default: 100)
        parameters (dict): override model parameters (optional)
        metabolites (list): only plot specific metabolites (optional)
        xlabel (str): specify label for x axis (optional)
        ylabel (str): specify label for y axis (optional)
        data (dict): metabolomics data (to overlay in the plot)
        data_steps (list): time_steps of metabolomics data

    Returns:
        AxesSubplot: axes handle from matplotlib

    """
    t, X = time_course(model, time, steps=steps, parameters=parameters)
    fig = figure()
    ax = fig.add_subplot(111)
    legend = model.metabolites.keys()
    if metabolites:
        indices = [model.metabolites.keys().index(m_id) for m_id in metabolites]
        X = X[:, indices]
        legend = metabolites
    ax.plot(t, X)
    ax.legend(legend)

    if data is not None and data_steps is not None:
        if not metabolites:
            metabolites = model.metabolites.keys()
        values = [data[m_id] for m_id in metabolites]
        ax.set_color_cycle(None)
        ax.plot(data_steps, array(values).T, '.')

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    return ax


def plot_flux_sampling(model, sample, reactions=None):
    """ Plot flux sampling results.

    Args:
        model (Model): your model (can be of any kind)
        sample (list): flux vector samples
        reactions (list): only plot specified reactions (default: all)

    """
    # TODO: this is generic enough to be kinetic/cobra compatible (move to core.plotting ?)

    try:
        from seaborn import kdeplot
    except:
        raise RuntimeError('to run this method please install seaborn')

    if not reactions:
        reactions = model.reactions.keys()

    sample = array(sample)
    n = len(reactions)
    margin = 0.3

    for i, rxn_y in enumerate(reactions):
        for j, rxn_x in enumerate(reactions):
            ax = subplot2grid((n, n), (n-1-i, j))
            x_data = sample[:, model.reactions.keys().index(rxn_x)]
            y_data = sample[:, model.reactions.keys().index(rxn_y)]

            x_min, x_max = min(x_data), max(x_data)
            x_delta = (x_max - x_min)*margin if x_max > x_min else x_max
            x_lim = (x_min - x_delta, x_max + x_delta)

            y_min, y_max = min(y_data), max(y_data)
            y_delta = (y_max - y_min)*margin if y_max > y_min else y_max
            y_lim = (y_min - y_delta, y_max + y_delta)

            if j > i:
                ax.scatter(x_data, y_data)
            elif j < i:
                kdeplot(x_data, y_data, cmap="Blues_d", ax=ax)
            else:
                kdeplot(x_data, ax=ax)
            if j == 0:
                ax.set_ylabel(rxn_y)
            else:
                ax.set_yticks([])

            if i == 0:
                ax.set_xlabel(rxn_x)
            else:
                ax.set_xticks([])

            ax.set_xlim(x_lim)
            if i != j:
                ax.set_ylim(y_lim)
