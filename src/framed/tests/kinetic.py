from framed.io_utils.sbml import load_odemodel
from scipy.integrate import odeint
from numpy import arange
from matplotlib.pyplot import plot, show
__author__ = 'daniel'


KINETIC_MODEL = '../../../examples/models/BIOMD0000000051.xml'


def main():
    model = load_odemodel(KINETIC_MODEL)
    x0 = model.concentrations.values()
    x0[1] = 2
    t = arange(0, 100, 0.1)
    p = model.get_p()
    p[model.indexed_params.keys().index(('vPGI', 'rmaxPGI'))] *= 0.01
    x = odeint(model.get_ODEs(p), x0, t)
    plot(t, x)
    show()



if __name__ == '__main__':
    main()
