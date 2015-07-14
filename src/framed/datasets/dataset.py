from collections import OrderedDict
from json import load, dump

__author__ = 'daniel'

class DataSet():

    def __init__(self):
        self.description = ''
        self.organism = ''
        self.timecourse = False
        self.conditions = []
        self.experiment_details = []
        self.files = {}


    def load(self, jsonfile):
        with open(jsonfile, 'r') as f:
            data = load(f)
            for key, value in data:
                setattr(self, key, value)


