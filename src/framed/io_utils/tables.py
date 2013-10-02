''' This module implements methods for reading and writing tables as CSV files.

@author: Daniel Machado

   Copyright 2013 Novo Nordisk Foundation Center for Biosustainability,
   Technical University of Denmark.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
   
'''

from csv import reader, writer

def read_csv_table(filename, row_names=False, column_names=False, datatype=str, delimiter=','):
    """ Reads a table from a csv file.
    
    Arguments:
        filename : str -- file path
        row_names : bool -- first column contains row names (default: False)
        column_names : bool -- first row contains column names (default: False)
        datatype : type -- type to cast table items (default: str)
        delimiter : str -- delimiter character (default: ',')

    Returns:
        [row names,] [column names,] data -- the data and optionally row and column names
    """  
    
    def cast(x):
        try:
            return datatype(x)
        except:
            return None
    
    with open(filename) as f:
        r = reader(f, delimiter=delimiter)
        table = [line for line in r]
    if row_names and column_names:
        rows = [line[0] for line in table[1:]]
        columns = table[0][1:]
        data = [map(cast, line[1:]) for line in table[1:]]
        return rows, columns, data
    elif row_names:
        rows = [line[0] for line in table]
        data = [map(cast, line[1:]) for line in table]
        return rows, data
    elif column_names:
        columns = table[0]
        data = [map(cast, line) for line in table[1:]]
        return columns, data
    else:
        data = [map(cast, line) for line in table]
        return data

def write_csv_table(filename, data, row_names=None, column_names=None, delimiter=','):
    """ Reads a table from a csv file.
    
    Arguments:
        filename : str -- file path
        data : list of list -- data to write
        row_names : list of str -- row names (optional)
        column_names : list of str -- column names (optional)
        delimiter : str -- delimiter character (default: ',')
    """  
    with open(filename, 'w') as f:
        w = writer(f, delimiter=delimiter)
        if column_names:
            if row_names:
                column_names = [''] + column_names
            w.writerow(column_names)
        if row_names:
            data = [[x] + y for x, y in zip(row_names, data)]
        for row in data:
            w.writerow(row)

