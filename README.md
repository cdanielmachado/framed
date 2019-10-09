[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1048261.svg)](https://doi.org/10.5281/zenodo.1048261) [![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0) [![PyPI version](https://badge.fury.io/py/framed.svg)](https://badge.fury.io/py/framed) [![Documentation Status](http://readthedocs.org/projects/framed/badge/?version=latest)](http://framed.readthedocs.io/en/latest/?badge=latest)

### Note: this package will soon be deprecated! Please check the new [ReFramed](https://github.com/cdanielmachado/reframed) package.

![FRAMED](logo_300px.png)

## A python FRAmework for Metabolic Engineering and Design

**FRAMED** is a python package for analysis and simulation of metabolic models. The main focus is to provide support for different modeling approaches. 

* Modeling: Constraint-based models, Kinetic models, Bioprocess models
* I/O: Import/Export from multiple SBML flavors and plain text formats (including BioOpt)
* Solver support: Gurobi, CPLEX
* COBRA tools:
    * Simulation: FBA, pFBA, loopless-FBA, MOMA, linearMOMA, ROOM
    * Gene-wise simulation: gene-pFBA, gene-MOMA, gene-lMOMA, gene-ROOM
    * Analysis: FVA, gene essentiality, PhPP, flux envelope plots
    * Ensemble-based simulation (includes import/export of ensemble models in SBML)
    * Omics integration: GIMME, E-Flux
* Kinetic tools:
    * Time-course and steady-state simulation
    * Steady-state flux sampling
    * Calibration from metabolomics data
* Bioprocess modeling:
    * Dynamic FBA (single and multi-species)


### Credits and License

Developed at:

* The Novo Nordisk Fundation Center for Biosustainability (2013)
* Centre of Biological Engineering, University of Minho (2014-2015)
* European Molecular Biology Laboratory (2016-2018)

Released under an Apache License.

