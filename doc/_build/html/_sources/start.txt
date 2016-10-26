===============
Getting started
===============


Loading models
--------------

*framed* supports different kinds of metabolic models. In any case, loading a model is quite simple.

For constraint-based (cobra) models:

::

    from framed import load_cbmodel
    model = load_cbmodel('my_model.xml')

Note that different people have been using different modeling conventions to store constraint-based models in SBML format.
*framed* handles this problem using *flavors*. We currently support two flavors:

- **cobra**: This is the format adopted by the cobra toolbox.
- **fbc2**: This is the new fbc2 extension for SBML.

If can optionally specify the *flavor* of the model you are loading and *framed* will automatically apply a set of rules
to improve compatibility (for example, removing boundary metabolites):

::

    model = load_cbmodel('my_model.xml', flavor='fbc2')


For kinetic models:

::

    from framed import load_odemodel
    model = load_odemodel('my_model.xml')


Saving models
-------------

If you apply any kind of changes to a model, you can easily save your model as follows:

::

    from framed import save_sbml_model
    save_sbml_model(model, 'my_output_file.xml')

For constraint-based models you can also specify a specific *flavor* (this is useful for converting between different *flavors*):

::

    save_sbml_model(model, 'my_output_file.xml', flavor='cobra')


Human readable format
~~~~~~~~~~~~~~~~~~~~~

Additionally, *framed* introduces a new shorthand notation for import/export of models in an human-readable format
(You can see this format by typing ``print model``).

::

    from framed import read_cbmodel_from_file, write_model_to_file
    model = read_cbmodel_from_file('my_file.txt')
    write_model_to_file(model, 'my_new_file.txt')

Note: This format is not yet available for kinetic models.


Model manipulation
------------------

Once you load a model you can easily access the model's attributes perform different kinds of operations using the *model* object.
Note that some operations might be specific to constraint-based or kinetic models.

You can quickly look at your model by printing it:

::

    print model

::

    R_ACALD: M_acald_c + M_coa_c + M_nad_c <-> M_accoa_c + M_h_c + M_nadh_c
    R_ACALDt: M_acald_e <-> M_acald_c
    R_ACKr: M_ac_c + M_atp_c <-> M_actp_c + M_adp_c
    R_ACONTa: M_cit_c <-> M_acon_C_c + M_h2o_c
    R_ACONTb: M_acon_C_c + M_h2o_c <-> M_icit_c
    (...)

You can also inspect particular metabolites and reactions.

::

    print model.reactions.R_PGI

::

    R_PGI: M_g6p_c <-> M_f6p_c

::

    print model.metabolites.M_g6p_c

::

    D-Glucose-6-phosphate

::

    print model.reactions.R_TKT1.gpr

::

    (G_b2465 or G_b2935)

You can easily make changes to a model. For instance, let's change the maximum glucose uptake rate:

::

    model.reactions.R_EX_glc_e.lb = -10

But for programmatic access, you may want to use reactions as a dictionary instead:

::

    for reaction in model.reactions.values():
        reaction.lb = 0

For convenience, you can access some methods from the model class directly:

::

    for r_id in model.reactions:
        model.set_lower_bound(r_id, 0)

You can easily add a new reaction to a model (or replace an existing reaction with the same *id*):

::

    model.add_reaction_from_str('R_new: A + 0.5 B --> 2 C')

For constraint-based models you can optionally define the flux bounds as well:

::

    model.add_reaction_from_str('R_new: S + M_atp_c <-> P + M_adp_c [-10, 10]')

You can also remove reactions (as well as metabolites, genes or compartments):

::

    model.remove_reaction('R_PGI')

::

    model.remove_metabolites(['M_h2o_c', 'M_h_c'])

There is a lot more you can try! Just take a look into our API.


Model transformations
~~~~~~~~~~~~~~~~~~~~~

There are some utilities to make some overall model transformations. One obvious one is to remove blocked reactions and
dead-end metabolites (only for constraint-based models):

::

    from framed import simplify
    simplify(model)

You can also create a fully irreversible model by splitting reversible reactions in two directions (useful for some applications):

::

    from framed import make_irreversible
    make_irreversible(model)

If you don't want to change your original model, in both cases you can generate a new model:

::

    new_model = simplify(model, inplace=False)
