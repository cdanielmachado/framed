from framed.io_utils.sbml import load_odemodel
import seaborn
from matplotlib.pyplot import show
from framed.kinetic.plotting import plot_simulation
from framed.kinetic.simulation import find_steady_state, simulate
__author__ = 'daniel'


KINETIC_MODEL = '../../../examples/models/Chassagnole2002_fixed.xml'


def main():
    model = load_odemodel(KINETIC_MODEL)
    model.concentrations['cglcex'] = 2.0
    model.local_parameters['vG6PDH']['rmaxG6PDH'] *= 0.1
    plot_simulation(model, 1e3, metabolites=['cglcex', 'cg6p', 'cpep', 'cpyr'],
                    xlabel='time', ylabel='concentration')
    show()

    fluxes = find_steady_state(model)
    print fluxes.items()



if __name__ == '__main__':
    main()
