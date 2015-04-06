__author__ = 'daniel'

from ..kinetic.simulation import simulate
from matplotlib.pyplot import figure

def plot_simulation(model, time, n_steps=1000, metabolites=None, xlabel=None, ylabel=None):
    t, X = simulate(model, time, n_steps)
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