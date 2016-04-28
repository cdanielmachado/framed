from framed.io_utils.sbml import load_odemodel
from matplotlib.pyplot import show
from framed.kinetic.fitting import fit_from_metabolomics
from framed.kinetic.plotting import plot_simulation, plot_sampling_results, plot_simulation_vs_data
from framed.kinetic.simulation import simulate

from framed.kinetic.sampling import sample
__author__ = 'daniel'


KINETIC_MODEL = '../../../examples/models/BIOMD0000000051.xml'
PREY_PREDATOR = '../../../examples/models/prey_predator.xml'


def run_simulation_and_plot():
    model = load_odemodel(KINETIC_MODEL)
    plot_simulation(model, 1e3, metabolites=['cglcex', 'cg6p', 'cpep', 'cpyr'],
                    xlabel='time', ylabel='concentration', parameters={'Dil': 0.2/3600})
    show()


def run_sampling():
    model = load_odemodel(KINETIC_MODEL)
    params = model.merge_constants()
    vmaxs = [p_id for p_id in params.keys() if 'max' in p_id]
    p_sample, v_sample = sample(model, 10,  parameters=vmaxs, scale='log2', dist='normal', dist_args=(0, 1))
    reactions = model.reactions.keys()[:3]
    plot_sampling_results(model, v_sample, reactions)
    show()


def run_calibration():
    model = load_odemodel(PREY_PREDATOR)
    perturbed = {'k1': 0.8, 'k2': 0.9, 'k3': 1.5}
    t, X = simulate(model, 10, steps=100, parameters=perturbed)
    all_data = dict(zip(model.metabolites.keys(), X.T))
    bounds = [(0, 10)]*3
    fitted = fit_from_metabolomics(model, t, all_data, bounds=bounds)
    plot_simulation_vs_data(model, t, all_data, parameters=fitted)
    show()


def main():
    run_simulation_and_plot()
    run_sampling()
    run_calibration()


if __name__ == '__main__':
    main()
