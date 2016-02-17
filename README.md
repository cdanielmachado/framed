FRAMED
======

FRAMED (**FRA**mework for **M**etabolic **E**ngineering and **D**esign)
is a python package for analysis and simulation of metabolic models.

This package is under development (no stable release yet). Some of the current functionality includes:

* Import/Export of SBML files (constraint-based and kinetic models)
* Import/Export from plain text (constraint-based only)
* Simulation methods: FBA, pFBA, MOMA, linearMOMA, ROOM, loopless-FBA, dynamicFBA (single and multi-species)
* Analysis tools: FVA, PhPP, production envelope plotting, gene essentiality analysis
* Reconstruction: gapFind, gapFill
* Rational strain design: combinatorial, greedy (gene and reaction-based)
* Omics integration: GIMME, E-Flux
* Kinetic models (experimental support): time-course and steady-state simulation, sampling, plotting utils

This is an overview of the package structure:

![package overview](https://raw.githubusercontent.com/cdanielmachado/framed/master/docs/package_overview.png)

### Credits and License

Initial development at the Novo Nordisk Fundation Center for Biosustainability (Technical University of Denmark)
and currently at the Centre of Biological Engineering (University of Minho). Released under an Apache License.

