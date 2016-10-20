""" This module implements methods for handling multi-omics datasets. 

Author: Daniel Machado

"""

import json
import pandas as pd
import os

__author__ = 'daniel'

MANDATORY_FIELDS = ['conditions',
                    'carbon_source',
                    'aerobic',
                    'gene_deletions',
                    'growth_rate']


class DataSet:
    """ Represents a multi-omics dataset for multiple experimental conditions.

    The dataset must be loaded from a json file with an acceptable format (see examples).

    Mandatory fields include:
        - conditions
        - carbon source
        - mode (aerobic/anaerobic)
        - gene deletions (if any)
        - measured growth rate

    Optionally these omics data may be available:
        - transcriptomics
        - proteomics
        - fluxomics
        - metabolomics
    """


    def __init__(self, jsonfile):
        try:
            with open(jsonfile, 'r') as f:
                self._data = json.load(f)
                self.path = os.path.dirname(jsonfile) + '/'

            for field in MANDATORY_FIELDS:
                setattr(self, field, pd.Series(self._data[field],
                                               index=self._data['conditions']))
        except Exception as e:
            print 'Failed to load dataset:', e

    def get_omics_data(self, kind, elements=None, conditions=None):
        if not hasattr(self, kind):
            try:
                datafile = self._data[kind]['datafile']
                data = pd.read_csv(self.path + datafile, index_col=0, header=0)
                setattr(self, kind, data)
            except:
                print 'Failed to load', kind, 'data'
                return

        data = getattr(self, kind)

        if elements and conditions:
            return data.loc[elements, conditions]
        elif conditions:
            return data.loc[:, conditions]
        elif elements:
            return data.loc[elements, :]
        else:
            return data

    def get_gene_deletions(self, condition):
        deletions = self.gene_deletions[condition]

        if deletions:
            return deletions.split()
        else:
            return []

    def get_fluxomics(self, reactions=None, conditions=None):
        return self.get_omics_data('fluxomics', reactions, conditions)

    def get_proteomics(self, proteins=None, conditions=None):
        return self.get_omics_data('proteomics', proteins, conditions)

    def get_metabolomics(self, metabolites=None, conditions=None):
        return self.get_omics_data('metabolomics', metabolites, conditions)

    def get_transcriptomics(self, genes=None, conditions=None):
        return self.get_omics_data('transcriptomics', genes, conditions)
