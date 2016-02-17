__author__ = 'daniel'

from ..kinetic.simulation import simulate
from matplotlib.pyplot import figure, subplot2grid
from seaborn import kdeplot
from numpy import array


def plot_simulation(model, time, steps=100, parameters=None, metabolites=None, xlabel=None, ylabel=None):
    t, X = simulate(model, time, steps=steps, parameters=parameters)
    fig = figure()
    ax = fig.add_subplot(111)
    legend = model.metabolites.keys()
    if metabolites:
        indices = [model.metabolites.keys().index(m_id) for m_id in metabolites]
        X = X[:, indices]
        legend = metabolites
    ax.plot(t, X)
    ax.legend(legend)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    return ax


def plot_simulation_vs_data(model, t_steps, data, parameters=None, metabolites=None, xlabel=None, ylabel=None):
    if not metabolites:
        metabolites = data.keys()

    ax = plot_simulation(model, t_steps[-1], parameters=parameters, metabolites=metabolites, xlabel=xlabel, ylabel=ylabel)
    data = array(data.values()).T
    ax.set_color_cycle(None)
    ax.plot(t_steps, data, '.')
    return ax


def plot_sampling_results(model, sample, reactions=None):
    # TODO: this is generic enough to be kinetic/cobra compatible (move to core.plotting ?)

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

