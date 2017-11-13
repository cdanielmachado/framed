|DOI| |License| |PyPI version| |Documentation Status|

FRAMED
======

*framed* is a python package for analysis and simulation of metabolic
models. The main focus is to provide support for different modeling
approaches.

-  Modeling: Constraint-based models, Kinetic models, Bioprocess models
-  I/O: Import/Export from SBML and other plain text formats (including
   BioOpt)
-  Solver support: Gurobi, CPLEX
-  COBRA tools:

   -  Simulation: FBA, pFBA, loopless-FBA, MOMA, linearMOMA, ROOM
   -  Gene-wise simulation: gene-pFBA, gene-MOMA, gene-lMOMA, gene-ROOM
   -  Analysis: FVA, gene essentiality, PhPP, flux envelope plots
   -  Ensemble-based simulation (and SBML import/export of ensemble
      models)
   -  Omics integration: GIMME, E-Flux
   -  Strain design: brute force, hill climbing

-  Kinetic tools:

   -  Time-course and steady-state simulation
   -  Steady-state flux sampling
   -  Calibration from metabolomics data

-  Bioprocess modeling:

   -  Dynamic FBA (single and multi-species)

-  Microbial community modeling:

   -  SMETANA

Documentation
~~~~~~~~~~~~~

For documentation and API please check: http://framed.readthedocs.io/

Instalation from PyPI (stable releases)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    pip install framed

Instalation from github (latest development release)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    pip install https://github.com/cdanielmachado/framed/archive/master.zip

Credits and License
~~~~~~~~~~~~~~~~~~~

Developed at:

-  The Novo Nordisk Fundation Center for Biosustainability (2013)
-  Centre of Biological Engineering, University of Minho (2014-2015)
-  European Molecular Biology Laboratory (2016-2017)

Released under an Apache License.

.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.240430.svg
   :target: https://doi.org/10.5281/zenodo.240430
.. |License| image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
   :target: https://opensource.org/licenses/Apache-2.0
.. |PyPI version| image:: https://badge.fury.io/py/framed.svg
   :target: https://badge.fury.io/py/framed
.. |Documentation Status| image:: http://readthedocs.org/projects/framed/badge/?version=latest
   :target: http://framed.readthedocs.io/en/latest/?badge=latest
