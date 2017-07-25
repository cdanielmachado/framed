[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.240430.svg)](https://doi.org/10.5281/zenodo.240430) [![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0) [![PyPI version](https://badge.fury.io/py/framed.svg)](https://badge.fury.io/py/framed) [![Documentation Status](http://readthedocs.org/projects/framed/badge/?version=latest)](http://framed.readthedocs.io/en/latest/?badge=latest)


FRAMED
======

*framed* is a python package for analysis and simulation of metabolic models. The main focus is to provide support for different modeling approaches. 

* Modeling: Constraint-based models, Kinetic models, Bioprocess models
* I/O: Import/Export from SBML and/or plain text formats
* Solver support: Gurobi, CPLEX, GLPK
* COBRA tools:
    * Simulation: FBA, pFBA, loopless-FBA, MOMA, linearMOMA, ROOM
    * Analysis: FVA, gene essentiality, PhPP, flux envelope plots
    * Omics integration: GIMME, E-Flux
    * Strain design: brute force, hill climbing
* Kinetic tools:
    * Time-course and steady-state simulation
    * Parameter and flux sampling
    * Calibration from metabolomics data
* Bioprocess modeling: Dynamic FBA (single and multi-species)

### Documentation

For documentation and API please check: http://framed.readthedocs.io/

### Instalation

```
pip install framed
```

### Credits and License

Developed at:

* The Novo Nordisk Fundation Center for Biosustainability (2013)
* Centre of Biological Engineering, University of Minho (2014-2015)
* European Molecular Biology Laboratory (2016-2017)

Released under an Apache License.

