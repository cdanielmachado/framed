import pandas as pd
from warnings import warn


class MetaNetX:
    """ This class implements a resource for mapping idenfitiers based on the MetaNetX database.

        It must be instantiated with the path to the local directory where MNX is stored.

    """

    def __init__(self, path_to_files, version=2.0):
        """ Build an instance of the MetaNetX translation service.

        Args:
            path_to_files (str): path to the local directory where MNX is stored.
        """
        self.path_to_files = path_to_files
        self.reac_xref = None
        self.chem_xref = None
        self.version = version

    def _get_reac_xref(self):
        if self.reac_xref is None:
            filename = self.path_to_files + '/reac_xref.tsv'
            self.reac_xref = pd.read_csv(filename, sep='\t', header=None, comment='#',
                                         names=['external_id', 'mnx_id'])
            if self.version >= 3.0:
                self.reac_xref = self.reac_xref[self.reac_xref['external_id'].str.contains(':')]

            self.reac_xref['db'] = self.reac_xref['external_id'].apply(lambda x: x.split(':')[0])
            self.reac_xref['local_id'] = self.reac_xref['external_id'].apply(lambda x: x.split(':')[1])

        return self.reac_xref

    def _get_chem_xref(self):
        if self.chem_xref is None:
            filename = self.path_to_files + '/chem_xref.tsv'
            self.chem_xref = pd.read_csv(filename, sep='\t', header=None, comment='#', usecols=[0, 1],
                                         names=['external_id', 'mnx_id'])

            if self.version >= 3.0:
                self.chem_xref = self.chem_xref[self.chem_xref['external_id'].str.contains(':')]

            self.chem_xref['db'] = self.chem_xref['external_id'].apply(lambda x: x.split(':')[0])
            self.chem_xref['local_id'] = self.chem_xref['external_id'].apply(lambda x: x.split(':')[1])
        return self.chem_xref

    def _slice_reac_xref(self, db):
        slice_name = 'reac_xref_' + db
        if not hasattr(self, slice_name):
            reac_xref = self._get_reac_xref()
            reac_xref_db = reac_xref.query('db == "{}"'.format(db))
            setattr(self, slice_name, reac_xref_db)
        return getattr(self, slice_name)

    def _slice_chem_xref(self, db):
        slice_name = 'chem_xref_' + db
        if not hasattr(self, slice_name):
            chem_xref = self._get_chem_xref()
            chem_xref_db = chem_xref.query('db == "{}"'.format(db))
            setattr(self, slice_name, chem_xref_db)
        return getattr(self, slice_name)

    def translate_reaction_id(self, rxn_id, from_db, to_db):
        """ Translate a reaction identifier.

        Args:
            rxn_id (str): original identifier
            from_db (str): source database
            to_db (str): target database

        Returns:
            list: zero or more possible translations for the identifier in the target database
        """

        dbs = ['rhea', 'kegg', 'metacyc', 'upa', 'seed', 'bigg', 'biopath', 'reactome']

        from_db = from_db.lower()
        to_db = to_db.lower()

        if from_db not in dbs:
            warn('Invalid source database: ' + from_db)
            return

        if to_db not in dbs:
            warn('Invalid source database: ' + to_db)
            return

        reac_xref_from_db = self._slice_reac_xref(from_db)
        reac_xref_to_db = self._slice_reac_xref(to_db)

        query1 = reac_xref_from_db.query('local_id == "{}"'.format(rxn_id))

        if len(query1) > 0:
            mnx_id = query1['mnx_id'].iloc[0]
            query2 = reac_xref_to_db.query('mnx_id == "{}"'.format(mnx_id))

            results = query2['local_id'].tolist()
        else:
            results = []

        return results

    def translate_metabolite_id(self, met_id, from_db, to_db):
        """ Translate a metabolite identifier.

        Args:
            met_id (str): original identifier
            from_db (str): source database
            to_db (str): target database

        Returns:
            list: zero or more possible translations for the identifier in the target database
        """
        dbs = ['chebi', 'kegg', 'metacyc', 'upa', 'seed', 'bigg', 'biopath', 'lipidmaps', 'hmdb', 'reactome', 'umbbd']

        from_db = from_db.lower()
        to_db = to_db.lower()

        if from_db not in dbs:
            warn('Invalid source database: ' + from_db)
            return

        if to_db not in dbs:
            warn('Invalid source database: ' + to_db)
            return

        chem_xref_from_db = self._slice_chem_xref(from_db)
        chem_xref_to_db = self._slice_chem_xref(to_db)

        query1 = chem_xref_from_db.query('local_id == "{}"'.format(met_id))

        if len(query1) > 0:
            mnx_id = query1['mnx_id'].iloc[0]
            query2 = chem_xref_to_db.query('mnx_id == "{}"'.format(mnx_id))

            results = query2['local_id'].tolist()
        else:
            results = []

        return results

    def model_translation(self, model, from_db, to_db, translate_reactions=True, translate_metabolites=True,
                          from_reaction_parser=None, from_metabolite_parser=None,
                          to_reaction_parser=None, to_metabolite_parser=None, inplace=False):
        """ Translate all possible identifiers in a model from one database to another

        Args:
            model (Model): original model
            from_db (str): source database
            to_db (str): target database
            translate_reactions (bool): translate reaction identifiers (default: True)
            translate_metabolites (bool): translate metabolite identifiers (default: True)
            from_reaction_parser (function): apply transformation to original reaction id (optional)
            from_metabolite_parser (function): apply transformation to original metabolite id (optional)
            to_reaction_parser (function): apply transformation to target reaction id (optional)
            to_metabolite_parser (function): apply transformation to target metabolite id (optional)
            inplace (bool): apply changes directly in the original model (default: False)

        Returns:
            Model: translated model

        Notes:
            If there is not a one-to-one mapping between the original and target identifiers, the substitutions
            will be unique and applied in alphabetic order.

            For example, if MNX reaction A matches [R1, R2] in *db1*, and [Rx, Ry, Rz] in *db2*,
            the applied translation from *db1* to *db2* will be {R1: Rx, R2: Ry}.
        """

        if not inplace:
            model = model.copy()

        if from_reaction_parser is None:
            from_reaction_parser = lambda x: x

        if from_metabolite_parser is None:
            from_metabolite_parser = lambda x: x

        if to_reaction_parser is None:
            to_reaction_parser = lambda _, x: x

        if to_metabolite_parser is None:
            to_metabolite_parser = lambda _, x: x

        met_dict = {}
        rxn_dict = {}

        def key_replace(ord_dict, key, new_key):
            if new_key in ord_dict:
                warn('Key already in dict, ignored:' + new_key)
            else:
                item = ord_dict[key]
                del ord_dict[key]
                ord_dict[new_key] = item

        if translate_metabolites:

            for m_id in sorted(model.metabolites):
                met_id = from_metabolite_parser(m_id)
                new_ids = self.translate_metabolite_id(met_id, from_db, to_db)
                new_ids = [to_metabolite_parser(m_id, new_id) for new_id in new_ids]
                new_ids = sorted(set(new_ids) - set(met_dict.values()))

                if len(new_ids) > 0:
                    model.metabolites[m_id].id = new_ids[0]
                    met_dict[m_id] = new_ids[0]

            for m_id, new_id in met_dict.items():
                key_replace(model.metabolites, m_id, new_id)

        for r_id in sorted(model.reactions):

            if translate_reactions:
                rxn_id = from_reaction_parser(r_id)
                new_ids = self.translate_reaction_id(rxn_id, from_db, to_db)
                new_ids = [to_reaction_parser(r_id, new_id) for new_id in new_ids]
                new_ids = sorted(set(new_ids) - set(rxn_dict.values()))

                model.reactions[r_id].synonyms = len(new_ids)

                if len(new_ids) > 0:
                    model.reactions[r_id].id = new_ids[0]
                    rxn_dict[r_id] = new_ids[0]

            if translate_metabolites:
                reaction = model.reactions[r_id]
                for m_id in reaction.stoichiometry.keys():
                    if m_id in met_dict:
                        key_replace(reaction.stoichiometry, m_id, met_dict[m_id])

        for r_id, new_id in rxn_dict.items():
            key_replace(model.reactions, r_id, new_id)

        if not inplace:
            return model

    def bigg2seed_model(self, model, inplace=False):
        """ Translate a model with BiGG identifiers to SEED identifiers

        Args:
            model (Model): original model
            inplace (bool): apply changes directly in the original model (default: False)

        Returns:
            Model: translated model
        """
        def from_reaction_parser(bigg_id):
            return bigg_id[2:]

        def from_metabolite_parser(bigg_id):
            return bigg_id[2:-2]

        def to_reaction_parser(_, seed_id):
            return seed_id + '_c0'

        def to_metabolite_parser(bigg_id, seed_id):
            return seed_id + bigg_id[-2:] + '0'

        return self.model_translation(model, 'bigg', 'seed',
                                      from_reaction_parser=from_reaction_parser,
                                      from_metabolite_parser=from_metabolite_parser,
                                      to_reaction_parser=to_reaction_parser,
                                      to_metabolite_parser=to_metabolite_parser,
                                      inplace=inplace)

    def seed2bigg_model(self, model, inplace=False):
        """ Translate a model with SEED identifiers to BiGG identifiers

        Args:
            model (Model): original model
            inplace (bool): apply changes directly in the original model (default: False)

        Returns:
            Model: translated model
        """

        def from_reaction_parser(seed_id):
            return seed_id[:-3]

        def from_metabolite_parser(seed_id):
            return seed_id[:-3]

        def to_reaction_parser(_, bigg_id):
            return 'R_' + bigg_id

        def to_metabolite_parser(seed_id, bigg_id):
            return 'M_' + bigg_id + '_' + seed_id[-2]

        return self.model_translation(model, 'seed', 'bigg',
                                      from_reaction_parser=from_reaction_parser,
                                      from_metabolite_parser=from_metabolite_parser,
                                      to_reaction_parser=to_reaction_parser,
                                      to_metabolite_parser=to_metabolite_parser,
                                      inplace=inplace)
