=============
Installation
=============

*framed* requires Python 2.7, which is available for all major operating systems. We recommend the `Anaconda python
distribution <https://www.continuum.io/downloads>`_.

*framed* can be easily installed using **pip**:

::

    pip install framed


Dependencies (automatically installed):

- libsbml (import/export of models in SBML format)
- sympy   (parsing GPR associations)
- scipy   (kinetic modeling)
- matplotlib (plotting methods)

Additional dependences (not automatically installed):

- seaborn (for some types of plots)

For constraint-based modeling you will need at least one solver installed. *framed* currently supports:

- Gurobi
- CPLEX


