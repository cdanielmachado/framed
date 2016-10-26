=============
Installation
=============

*framed* requires Python 2.7, which is available for all major operating systems. We recommend the `Anaconda python
distribution <https://www.continuum.io/downloads>`_.

*framed* can be easily installed using **pip**:

::

    pip install framed


Additional dependencies:

- libsbml (import/export of models in SBML format)
- sympy   (parsing GPR associations)
- scipy   (kinetic modeling)
- seaborn (plotting methods)

Also, for any kind of constraint-based modeling you will need at least one solver installed. *framed* currently supports:

- Gurobi
- CPLEX
- GLPK (note: this interface is no longer maintained)


