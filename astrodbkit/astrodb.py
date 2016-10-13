#!/usr/bin/python
# encoding: utf-8
# Author: Joe Filippazzo, jcfilippazzo@gmail.com

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import io
import os
import sys
import sqlite3
import warnings
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import astropy.io.ascii as ii
import astropy.table as at
from astropy.utils.data import download_file
from . import votools

warnings.simplefilter('ignore')

# Set get_input to be either raw_input or input depending on the python version
if sys.version_info[0] >= 3:
    get_input = input  # Python 3
else:
    get_input = raw_input


def create_database(dbpath):
    """
    Create a new database at the given dbpath

    Parameters
    ----------
    dbpath: str
        The full path for the new database, including the filename and .db file extension.

  """
    if dbpath.endswith('.db'):
        sources_table = "CREATE TABLE sources (id INTEGER PRIMARY KEY, ra REAL, dec REAL, designation TEXT, " \
                        "publication_id INTEGER, shortname TEXT, names TEXT, comments TEXT)"
        os.system("sqlite3 {} '{}'".format(dbpath, sources_table))
        if os.path.isfile(dbpath):
            print(
                "\nDatabase created! To load, run\n\ndb = astrodb.Database('{}')"
                "\n\nThen run db.modify_table() method to create tables.".format(dbpath))
    else:
        print("Please provide a path and file name with a .db file extension, e.g. /Users/<username>/Desktop/test.db")


class Database:
    def __init__(self, dbpath, directory='tabledata'):
        """
        Initialize the database.

        Parameters
        ----------
        dbpath: str
            The path to the .db or .sql database file.
        directory: str
            Folder in which individual tables are stored (Default: tabledata)

        Returns
        -------
        object
            The database object

        """
        if os.path.isfile(dbpath):

            # Alternatively, just list the directory with the schema and .sql files and require that dbpath
            # is the schema file, then load the tables individually

            # If it is a .sql file, create an empty database in the
            # working directory and generate the database from file
            if dbpath.endswith('.sql'):
                self.sqlpath = dbpath
                self.dbpath = dbpath.replace('.sql', '.db')

                # If the .db file already exists, rename it with the date
                if os.path.isfile(self.dbpath):
                    import datetime
                    date = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M")
                    print("Renaming existing file {} to {}".format(self.dbpath, self.dbpath.replace('.db', date+'.db')))
                    os.system("mv {} {}".format(self.dbpath, self.dbpath.replace('.db', date+'.db')))

                # Make the new database from the .sql files
                # First the schema...
                os.system("sqlite3 {} < {}".format(self.dbpath, self.sqlpath))

                # Then load the table data...
                print('Populating database...')
                tables = os.popen('sqlite3 {} ".tables"'.format(self.dbpath)).read().replace('\n',' ').split()
                for table in tables:
                    os.system('sqlite3 {0} ".read {1}/{2}.sql"'.format(self.dbpath, directory, table))
            elif dbpath.endswith('.db'):
                self.sqlpath = dbpath.replace('.db', '.sql')
                self.dbpath = dbpath
            else:
                self.sqlpath = dbpath + '.sql'
                self.dbpath = dbpath

            # Create .sql schema file if it doesn't exist
            os.system('touch {}'.format(self.sqlpath.replace(' ', '\ ')))

            # Create connection
            con = sqlite3.connect(self.dbpath, isolation_level=None, detect_types=sqlite3.PARSE_DECLTYPES)
            con.text_factory = sqlite3.OptimizedUnicode
            self.conn = con
            self.curs = con.cursor()
            self.list = con.cursor().execute

            # Make dictionary
            def dict_factory(cursor, row):
                d = {}
                for idx, col in enumerate(cursor.description):
                    d[col[0]] = row[idx]
                return d

            self.dict = con.cursor()
            self.dict.row_factory = dict_factory
            self.dict = self.dict.execute

            # Make sure the ignore table exists
            self.list(
                "CREATE TABLE IF NOT EXISTS ignore (id INTEGER PRIMARY KEY, id1 INTEGER, id2 INTEGER, tablename TEXT)")

            # Activate foreign key support
            self.list('PRAGMA foreign_keys=ON')

            print("Database ready for use")

        else:
            raise AttributeError("Sorry, no such file '{}'".format(dbpath))

    def __repr__(self):
        self.info()
        print("\nFor a quick summary of how to use astrodb.Database, type db.help(), \n"
              "where 'db' corresponds to the name of the astrodb.Database instance.")
        return ''

    def add_data(self, data, table, delimiter='|', bands='', verbose=False):
        """
    Adds data to the specified database table. Column names must match table fields to insert,
    however order and completeness don't matter.

    Parameters
    ----------
    data: str, sequence
      The path to an ascii file or a list of lists. The first row or element must
      be the list of column names
    table: str
      The name of the table into which the data should be inserted
    delimiter: str
      The string to use as the delimiter when parsing the ascii file
    bands: sequence
      Sequence of band to look for in the data header when digesting columns of
      multiple photometric measurements (e.g. ['MKO_J','MKO_H','MKO_K']) into individual
      rows of data for database insertion
    verbose: bool
      Print diagnostic messages

    """
        # Store raw entry
        entry, del_records = data, []

        # Digest the ascii file into table
        if isinstance(data, str) and os.path.isfile(data):
            data = ii.read(data, delimiter=delimiter)

        # Or read the sequence of data elements into a table
        elif isinstance(data, (list, tuple, np.ndarray)):
            data = ii.read(['|'.join(map(str, row)) for row in data], data_start=1, delimiter='|')

        else:
            data = None

        if data:

            # Get list of all columns and make an empty table for new records
            metadata = self.query("PRAGMA table_info({})".format(table), fmt='table')
            columns, types, required = [np.array(metadata[n]) for n in ['name', 'type', 'notnull']]
            new_records = at.Table(names=columns, dtype=[type_dict[t] for t in types])

            # Convert data dtypes to those of the existing table
            for col in data.colnames:
                try:
                    temp = data[col].astype(new_records[col].dtype)
                    data.replace_column(col, temp)
                except (KeyError, AttributeError):
                    continue

            # If a row contains photometry for multiple bands, use the *multiband argument and execute this
            if bands and table.lower() == 'photometry':

                # Pull out columns that are band names
                for b in list(set(bands) & set(data.colnames)):
                    try:
                        # Get the repeated data plus the band data and rename the columns
                        band = data[list(set(columns) & set(data.colnames)) + [b, b + '_unc']]
                        for suf in ['', '_unc']:
                            band.rename_column(b + suf, 'magnitude' + suf)
                        band.add_column(at.Column([b] * len(band), name='band'))

                        # Add the band data to the list of new_records
                        new_records = at.vstack([new_records, band])
                    except IOError:
                        pass

            else:
                # Inject data into full database table format
                new_records = at.vstack([new_records, data])[new_records.colnames]

            # Reject rows that fail column requirements, e.g. NOT NULL fields like 'source_id'
            for r in columns[np.where(np.logical_and(required, columns != 'id'))]:
                # Null values...
                new_records = new_records[np.where(new_records[r])]

                # Masked values...
                new_records = new_records[~new_records[r].mask]

                # NaN values...
                if new_records.dtype[r] in (int, float):
                    new_records = new_records[~np.isnan(new_records[r])]

            # For spectra, try to populate the table by reading the FITS header
            if table.lower() == 'spectra':
                for n, new_rec in enumerate(new_records):

                    # Convert relative path to absolute path
                    relpath = new_rec['spectrum']
                    if relpath.startswith('$'):
                        abspath = os.popen('echo {}'.format(relpath.split('/')[0])).read()[:-1]
                        if abspath:
                            new_rec['spectrum'] = relpath.replace(relpath.split('/')[0], abspath)

                    # Test if the file exists and try to pull metadata from the FITS header
                    if os.path.isfile(new_rec['spectrum']):
                        new_records[n]['spectrum'] = relpath
                        new_records[n] = _autofill_spec_record(new_rec)
                    else:
                        print('Error adding the spectrum at {}'.format(new_rec['spectrum']))
                        del_records.append(n)

                # Remove bad records from the table
                new_records.remove_rows(del_records)

            # Get some new row ids for the good records
            rowids = self._lowest_rowids(table, len(new_records))

            # Add the new records
            for N, new_rec in enumerate(new_records):
                new_rec = list(new_rec)
                new_rec[0] = rowids[N]
                for n, col in enumerate(new_rec):
                    if type(col) == np.int64 and sys.version_info[0] >= 3:
                        # Fix for Py3 and sqlite3 issue with numpy types
                        new_rec[n] = col.item()
                    if type(col) == np.ma.core.MaskedConstant:
                        new_rec[n] = None
                self.modify("INSERT INTO {} VALUES({})".format(table, ','.join('?' * len(columns))), new_rec)
                new_records[N]['id'] = rowids[N]

            # print a table of the new records or bad news
            if new_records:
                pprint(new_records, names=columns,
                       title="{} new records added to the {} table.".format(len(new_records), table.upper()))
            else:
                print('No new records added to the {} table. Please check your input: {}'.format(table, entry))

            # Run table clean up
            self.clean_up(table, verbose)

        else:
            print('Please check your input: {}'.format(entry))

    def add_foreign_key(self, table, parent, key_child, key_parent, verbose=True):
        """
        Add foreign key (**key_parent** from **parent**) to **table** column **key_child**

        Parameters
        ----------
        table: string
            The name of the table to modify. This is the child table.
        parent: string or list of strings
            The name of the reference table. This is the parent table.
        key_child: string or list of strings
            Column in **table** to set as foreign key. This is the child key.
        key_parent: string or list of strings
            Column in **parent** that the foreign key refers to. This is the parent key.
        verbose: bool, optional
            Verbose output
        """

        # Temporarily turn off foreign keys
        self.list('PRAGMA foreign_keys=OFF')

        metadata = self.query("PRAGMA table_info({})".format(table), fmt='table')
        columns, types, required, pk = [np.array(metadata[n]) for n in ['name', 'type', 'notnull', 'pk']]

        # Set constraints
        constraints = []
        for elem in required:
            if elem > 0:
                constraints.append('NOT NULL')
            else:
                constraints.append('')

        # Set PRIMARY KEY columns
        ind, = np.where(pk >= 1)
        for i in ind:
            constraints[i] += ' UNIQUE'  # Add UNIQUE constraint to primary keys
        pk_names = columns[ind]

        try:
            # Rename the old table and create a new one
            self.list("DROP TABLE IF EXISTS TempOldTable_foreign")
            self.list("ALTER TABLE {0} RENAME TO TempOldTable_foreign".format(table))

            # Re-create the table specifying the FOREIGN KEY
            sqltxt = "CREATE TABLE {0} ({1}".format(table, ', '.join(['{} {} {}'.format(c, t, r)
                                                                      for c, t, r in zip(columns, types, constraints)]))
            sqltxt += ', PRIMARY KEY({})'.format(', '.join([elem for elem in pk_names]))
            if isinstance(key_child, type(list())):
                for kc, p, kp in zip(key_child, parent, key_parent):
                    sqltxt += ', FOREIGN KEY ({0}) REFERENCES {1} ({2}) ON UPDATE CASCADE'.format(kc, p, kp)
            else:
                sqltxt += ', FOREIGN KEY ({0}) REFERENCES {1} ({2}) ON UPDATE CASCADE'.format(key_child, parent, key_parent)
            sqltxt += ' )'

            self.list(sqltxt)

            # Populate the new table and drop the old one
            tempdata = self.query("PRAGMA table_info(TempOldTable_foreign)", fmt='table')
            old_columns = [c for c in tempdata['name'] if c in columns]
            self.list("INSERT INTO {0} ({1}) SELECT {1} FROM TempOldTable_foreign".format(table, ','.join(old_columns)))
            self.list("DROP TABLE TempOldTable_foreign")

            if verbose:
                # print('Successfully added foreign key.')
                t = self.query('SELECT name, sql FROM sqlite_master', fmt='table')
                print(t[t['name'] == table]['sql'][0].replace(',', ',\n'))

        except:
            print('Error attempting to add foreign key.')
            self.list("DROP TABLE IF EXISTS {0}".format(table))
            self.list("ALTER TABLE TempOldTable_foreign RENAME TO {0}".format(table))
            raise sqlite3.IntegrityError('Failed to add foreign key')

        # Reactivate foreign keys
        self.list('PRAGMA foreign_keys=ON')

    def clean_up(self, table, verbose=False):
        """
        Removes exact duplicates, blank records or data without a *source_id* from the specified **table**.
        Then finds possible duplicates and prompts for conflict resolution.

        Parameters
        ----------
        table: str
            The name of the table to remove duplicates, blanks, and data without source attributions.
        verbose: bool
            Print out some diagnostic messages

        """
        # Get the table info and all the records
        metadata = self.query("PRAGMA table_info({})".format(table), fmt='table')
        columns, types, required = [np.array(metadata[n]) for n in ['name', 'type', 'notnull']]
        # records = self.query("SELECT * FROM {}".format(table), fmt='table', use_converters=False)
        ignore = self.query("SELECT * FROM ignore WHERE tablename LIKE ?", (table,))
        duplicate, command = [1], ''

        # Remove records with missing required values
        req_keys = columns[np.where(required)]
        try:
            self.modify("DELETE FROM {} WHERE {}".format(table, ' OR '.join([i + ' IS NULL' for i in req_keys])),
                        verbose=False)
            self.modify(
                "DELETE FROM {} WHERE {}".format(table, ' OR '.join([i + " IN ('null','None','')" for i in req_keys])),
                verbose=False)
        except:
            pass

        # Remove exact duplicates
        self.modify("DELETE FROM {0} WHERE id NOT IN (SELECT min(id) FROM {0} GROUP BY {1})".format(table, ', '.join(
                columns[1:])), verbose=False)

        # Check for records with identical required values but different ids.
        if table.lower() != 'sources': req_keys = columns[np.where(np.logical_and(required, columns != 'id'))]

        # List of old and new pairs to ignore
        if not type(ignore) == np.ndarray: ignore = np.array([])
        new_ignore = []

        while any(duplicate):
            # Pull out duplicates one by one
            if 'source_id' not in columns:  # Check if there is a source_id in the columns
                SQL = "SELECT t1.id, t2.id FROM {0} t1 JOIN {0} t2 ON t1.id=t2.id WHERE ".format(table)
            else:
                SQL = "SELECT t1.id, t2.id FROM {0} t1 JOIN {0} t2 ON t1.source_id=t2.source_id " \
                      "WHERE t1.id!=t2.id AND ".format(table)

            if any(req_keys):
                SQL += ' AND '.join(['t1.{0}=t2.{0}'.format(i) for i in req_keys]) + ' AND '

            if any(ignore):
                SQL += ' AND '.join(
                    ["(t1.id NOT IN ({0}) AND t2.id NOT IN ({0}))".format(','.join(map(str, [id1, id2])))
                        for id1, id2 in zip(ignore['id1'], ignore['id2'])]
                    if any(ignore) else '') + ' AND '

            if any(new_ignore):
                SQL += ' AND '.join(
                    ["(t1.id NOT IN ({0}) AND t2.id NOT IN ({0}))".format(','.join(map(str, ni)))
                     for ni in new_ignore] if new_ignore else '') + ' AND '

            # Clean up empty WHERE at end if it's present (eg, for empty req_keys, ignore, and new_ignore)
            if SQL[-6:] == 'WHERE ':
                SQL = SQL[:-6]

            # Clean up hanging AND if present
            if SQL[-5:] == ' AND ':
                SQL = SQL[:-5]

            if verbose:
                print('\nSearching for duplicates with: {}\n'.format(SQL))

            duplicate = self.query(SQL, fetch='one')

            # Compare potential duplicates and prompt user for action on each
            try:
                # Run record matches through comparison and return the command
                command = self._compare_records(table, duplicate)

                # Add acceptable duplicates to ignore list or abort
                if command == 'keep':
                    new_ignore.append([duplicate[0], duplicate[1]])
                    self.list("INSERT INTO ignore VALUES(?,?,?,?)", (None, duplicate[0], duplicate[1], table.lower()))
                elif command == 'undo':
                    pass  # TODO: Add this functionality!
                elif command == 'abort':
                    break
                else:
                    pass
            except:
                break

        # Finish or abort table clean up
        if command == 'abort':
            print('\nAborted clean up of {} table.'.format(table.upper()))
            return 'abort'
        else:
            print('\nFinished clean up on {} table.'.format(table.upper()))

    # @property
    def close(self, silent=False, directory='tabledata'):
        """
        Close the database and ask to save and delete the file

        Parameters
        ----------
        silent: bool
            Close quietly without saving or deleting (Default: False).
        """
        if not silent:
            saveme = get_input("Save database contents to '{}/'? (y, [n]) \n"
                               "To save under a folder name, run db.save() before closing. ".format(directory))
            if saveme.lower() == 'y':
                self.save()

            delete = get_input("Do you want to delete {0}? (y,[n]) \n"
                               "Don't worry, a new one will be generated if you run astrodb.Database('{1}') "
                               .format(self.dbpath, self.sqlpath))
            if delete.lower() == 'y':
                print("Deleting {}".format(self.dbpath))
                os.system("rm {}".format(self.dbpath))

        print('Closing connection')
        self.conn.close()

    def _compare_records(self, table, duplicate, options=['r', 'c', 'k', 'sql']):
        """
        Compares similar records and prompts the user to make decisions about keeping, updating, or modifying records in question.

        Parameters
        ----------
        table: str
            The name of the table whose records are being compared.
        duplicate: sequence
            The ids of the potentially duplicate records
        options: list
            The allowed options: 'r' for replace, 'c' for complete, 'k' for keep, 'sql' for raw SQL input.

        """
        # print the old and new records suspected of being duplicates
        verbose = True
        if duplicate[0] == duplicate[1]:  # No need to display if no duplicates were found
            verbose = False

        data = self.query("SELECT * FROM {} WHERE id IN ({})".format(table, ','.join(map(str, duplicate))), \
                          fmt='table', verbose=verbose, use_converters=False)
        columns = data.colnames[1:]
        old, new = [[data[n][k] for k in columns[1:]] for n in [0, 1]]

        # Prompt the user for action
        replace = get_input(
            "\nKeep both records [k]? Or replace [r], complete [c], or keep only [Press *Enter*] record {}? (Type column name to inspect or 'help' for options): ".format(
                    duplicate[0])).lower()
        replace = replace.strip()

        while replace in columns or replace == 'help':
            if replace in columns:
                pprint(np.asarray([[i for idx, i in enumerate(old) if idx in [0, columns.index(replace)]], \
                                   [i for idx, i in enumerate(new) if idx in [0, columns.index(replace)]]]), \
                       names=['id', replace])

            elif replace == 'help':
                _help()

            replace = get_input(
                "\nKeep both records [k]? Or replace [r], complete [c], or keep only [Press *Enter*] record {}? (Type column name to inspect or 'help' for options): ".format(
                        duplicate[0])).lower()

        if replace and replace.split()[0] in options:

            # Replace the entire old record with the new record
            if replace == 'r':
                sure = get_input(
                    'Are you sure you want to replace record {} with record {}? [y/n] : '.format(*duplicate))
                if sure.lower() == 'y':
                    self.modify("DELETE FROM {} WHERE id={}".format(table, duplicate[0]), verbose=False)
                    self.modify("UPDATE {} SET id={} WHERE id={}".format(table, duplicate[0], duplicate[1]),
                                verbose=False)

            # Replace specific columns
            elif replace.startswith('r'):
                replace_cols = replace.split()[1:]
                if all([i in columns for i in replace_cols]):
                    empty_cols, new_vals = zip(
                        *[['{}=?'.format(e), n] for e, n in zip(columns, new) if e in replace_cols])
                    if empty_cols:
                        self.modify("DELETE FROM {} WHERE id={}".format(table, duplicate[1]), verbose=False)
                        self.modify("UPDATE {} SET {} WHERE id={}".format(table, ','.join(empty_cols), duplicate[0]),
                                    tuple(new_vals), verbose=False)
                else:
                    badcols = ','.join([i for i in replace_cols if i not in columns])
                    print("\nInvalid column names for {} table: {}".format(table, badcols))

            # Complete the old record with any missing data provided in the new record, then delete the new record
            elif replace == 'c':
                try:
                    empty_cols, new_vals = zip(
                        *[['{}=?'.format(e), n] for e, o, n in zip(columns[1:], old, new) if n and not o])
                    self.modify("DELETE FROM {} WHERE id={}".format(table, duplicate[1]), verbose=False)
                    self.modify("UPDATE {} SET {} WHERE id={}".format(table, ','.join(empty_cols), duplicate[0]),
                                tuple(new_vals), verbose=False)
                except:
                    pass

            # Keep both records
            elif replace == 'k':
                return 'keep'

            # Execute raw SQL
            elif replace.startswith('sql') and 'sql' in options:
                self.modify(replace[4:], verbose=False)

        # Abort the current database clean up
        elif replace == 'abort':
            return 'abort'

        # Delete the higher id record
        elif not replace:
            self.modify("DELETE FROM {} WHERE id={}".format(table, max(duplicate)), verbose=False)

        # Prompt again
        else:
            print("\nInvalid command: {}\nTry again or type 'help' or 'abort'.\n".format(replace))

    def _explicit_query(self, SQL, use_converters=True):
        """
        Sorts the column names so they are returned in the same order they are queried. Also turns
        ambiguous SELECT statements into explicit SQLite language in case column names are not unique.

        HERE BE DRAGONS!!! Bad bad bad. This method needs to be reworked.

        Parameters
        ----------
        SQL: str
            The SQLite query to parse
        use_converters: bool
            Apply converters to columns with custom data types

        Returns
        -------
        (SQL, columns): (str, sequence)
            The new SQLite string to use in the query and the ordered column names

        """
        try:
            # If field names are given, sort so that they come out in the same order they are fetched
            if 'select' in SQL.lower() and 'from' in SQL.lower():

                # Make a dictionary of the table aliases
                tdict = {}
                from_clause = SQL.lower().split('from ')[-1].split(' where')[0]
                tables = [j for k in [i.split(' on ') for i in from_clause.split(' join ')] for j in k if '=' not in j]

                for t in tables:
                    t = t.replace('as', '')
                    try:
                        name, alias = t.split()
                        tdict[alias] = name
                    except:
                        tdict[t] = t

                # Get all the column names and dtype placeholders
                columns = \
                SQL.replace(' ', '').lower().split('distinct' if 'distinct' in SQL.lower() else 'select')[1].split(
                    'from')[0].split(',')

                # Replace * with the field names
                for n, col in enumerate(columns):
                    if '.' in col:
                        t, col = col.split('.')
                    else:
                        t = tables[0]

                    if '*' in col:
                        col = np.array(self.list("PRAGMA table_info({})".format(tdict.get(t))).fetchall()).T[1]
                    else:
                        col = [col]

                    columns[n] = ["{}.{}".format(t, c) if len(tables) > 1 else c for c in col]

                # Flatten the list of columns and dtypes
                columns = [j for k in columns for j in k]

                # Get the dtypes
                dSQL = "SELECT " \
                       + ','.join(["typeof({})".format(col) for col in columns]) \
                       + ' FROM ' + SQL.replace('from', 'FROM').split('FROM')[-1]
                if use_converters:
                    dtypes = [None] * len(columns)
                else:
                    dtypes = self.list(dSQL).fetchone()

                # Reconstruct SQL query
                SQL = "SELECT {}".format('DISTINCT ' if 'distinct' in SQL.lower() else '') \
                      + (','.join(["{0} AS '{0}'".format(col) for col in columns]) \
                             if use_converters else ','.join(["{1}{0}{2} AS '{0}'".format(col,
                                                                                          'CAST(' if dt != 'null' else '',
                                                                                          ' AS {})'.format(
                                                                                              dt) if dt != 'null' else '') \
                                                              for dt, col in zip(dtypes, columns)])) \
                      + ' FROM ' \
                      + SQL.replace('from', 'FROM').split('FROM')[-1]

            elif 'pragma' in SQL.lower():
                columns = ['cid', 'name', 'type', 'notnull', 'dflt_value', 'pk']

            return SQL, columns

        except:
            return SQL, ''

    def help(self):
        """
        See a quick summary of the most useful methods in astrodb.Database
        """

        helptext = """
The astrodb.Database class, hereafter db, provides a variety of methods to interact with a SQLite database file.
Docstrings are available for all methods and can be accessed in the usual manner; eg, help(db.query).
We list a few key methods below.

    * db.query() - send SELECT commands to the database. Returns results in a variety of formats
    * db.add_data() - add data to an existing table, either by providing a file or by providing the data itself
    * db.table() - create or modify tables in the database
    * db.modify() - send more general SQL commands to the database
    * db.info() - get a quick summary of the contents of the database
    * db.schema() - quickly examine the columns, types, etc of a specified table
    * db.search() - search through a table to find entries matching the criteria
    * db.references() - search for all entries in all tables matching the criteria. Useful for publication
    * db.save() - export a copy of the database in ascii format, which can then be re-populated by astrodb.Database
    * db.close() - close the database connection, will prompt to save and to delete the binary database file

The full documentation can be found online at: http://astrodbkit.readthedocs.io/en/latest/index.html
        """
        print(helptext)

    def info(self):
        """
        Prints out information for the loaded database, namely the available tables and the number of entries for each.
        """
        t = self.query("SELECT * FROM sqlite_master WHERE type='table'", fmt='table')
        all_tables = t['name'].tolist()
        print('\nDatabase path: {} \nSQL path: {}\n'.format(self.dbpath, self.sqlpath))
        print('Database Inventory')
        print('==================')
        for table in ['sources'] + [t for t in all_tables if
                                    t not in ['sources', 'sqlite_sequence']]:
            x = self.query('select count() from {}'.format(table), fmt='array', fetch='one')
            if x is None: continue
            print('{}: {}'.format(table.upper(), x[0]))

    def inventory(self, source_id, fetch=False, fmt='table'):
        """
        Prints a summary of all objects in the database. Input string or list of strings in **ID** or **unum**
        for specific objects.

        Parameters
        ----------
        source_id: int
            The id from the SOURCES table whose data across all tables is to be printed.
        fetch: bool
            Return the results.
        fmt: str
            Returns the data as a dictionary, array, or astropy.table given 'dict', 'array', or 'table'

        Returns
        -------
        data_tables: dict
            Returns a dictionary of astropy tables with the table name as the keys.

        """
        data_tables = {}

        t = self.query("SELECT * FROM sqlite_master WHERE type='table'", fmt='table')
        all_tables = t['name'].tolist()
        for table in ['sources'] + [t for t in all_tables if
                                    t not in ['sources', 'sqlite_sequence']]:

            try:

                # Get the columns, pull out redundant ones, and query the table for this source's data
                t = self.query("PRAGMA table_info({})".format(table), fmt='table')
                columns = np.array(t['name'])
                types = np.array(t['type'])

                if table == 'sources' or 'source_id' in columns:

                    # Only get simple data types and exclude redundant 'source_id' for nicer printing
                    columns = columns[
                        ((types == 'REAL') | (types == 'INTEGER') | (types == 'TEXT')) & (columns != 'source_id')]

                    # Query the table
                    try:
                        id = 'id' if table.lower() == 'sources' else 'source_id'
                        data = self.query(
                            "SELECT {} FROM {} WHERE {}={}".format(','.join(columns), table, id, source_id),
                            fmt='table')

                        if not data and table.lower() == 'sources':
                            print(
                            'No source with id {}. Try db.search() to search the database for a source_id.'.format(
                                source_id))

                    except:
                        data = None

                    # If there's data for this table, save it
                    if data:
                        if fetch:
                            data_tables[table] = self.query(
                                "SELECT {} FROM {} WHERE {}={}".format(','.join(columns), table, id, source_id), \
                                fetch=True, fmt=fmt)
                        else:
                            data = data[list(columns)]
                            pprint(data, title=table.upper())

                else:
                    pass

            except:
                print('Could not retrieve data from {} table.'.format(table.upper()))

        if fetch: return data_tables

    def lookup(self, criteria, table, columns=''):
        """
        Returns a table of records from *table* the same length as *criteria*
        with the best match for each element.
        
        Parameters
        ----------
        criteria: sequence 
            The search criteria
        table: str
            The table to search
        columns: sequence
            The column name in the sources table to search
        
        Returns
        -------
        results: sequence
            A sequence the same length as objlist with source_ids that correspond 
            to successful matches and blanks where no matches could be made
        """
        results, colmasks = [], []
        
        # Iterate through the list, trying to match objects
        for n,criterion in enumerate(criteria):
            records = self.search(criterion, table, columns=columns, fetch=True)
                
            # If multiple matches, take the first but notify the user of the other matches
            if len(records)>1:
                print("'{}' matched to {} other record{}.".format(criterion, len(records)-1, \
                      's' if len(records)-1>1 else ''))
            
            # If no matches, make an empty row
            if len(records)==0:
                records.add_row(np.asarray(np.zeros(len(records.colnames))).T)
                colmasks.append([True]*len(records.colnames))
            else:
                colmasks.append([False]*len(records.colnames))
            
            # Grab the first row
            results.append(records[0])
        
        # Add all the rows to the results table
        table = at.Table(rows=results, names=results[0].colnames, masked=True)
        
        # Mask the rows with no matches
        for col,msk in zip(records.colnames,np.asarray(colmasks).T): 
            table[col].mask = msk
        
        return table

    def _lowest_rowids(self, table, limit):
        """
        Gets the lowest available row ids for table insertion. Keeps things tidy!

        Parameters
        ----------
        table: str
            The name of the table being modified
        limit: int
            The number of row ids needed

        Returns
            -------
        available: sequence
            An array of all available row ids

        """
        try:
            t = self.query("SELECT id FROM {}".format(table), unpack=True, fmt='table')
            ids = t['id']
            all_ids = np.array(range(1, max(ids)))
        except TypeError:
            ids = None
            all_ids = np.array(range(1, limit+1))

        available = all_ids[np.in1d(all_ids, ids, assume_unique=True, invert=True)][:limit]

        # If there aren't enough empty row ids, start using the new ones
        if len(available) < limit:
            diff = limit - len(available)
            available = np.concatenate((available, np.array(range(max(ids) + 1, max(ids) + 1 + diff))))

        return available

    def merge(self, conflicted, tables=[], diff_only=True):
        """
        Merges specific **tables** or all tables of **conflicted** databse into the master database.

        Parameters
        ----------
        conflicted: str
            The path of the SQL database to be merged into the master.
        tables: list (optional)
            The list of tables to merge. If None, all tables are merged.
        diff_only: bool
                If True, only prints the differences of each table and doesn't actually merge anything.

        """
        if os.path.isfile(conflicted):
            # Load and attach master and conflicted databases
            con, master, reassign = Database(conflicted), self.list("PRAGMA database_list").fetchall()[0][2], {}
            con.modify("ATTACH DATABASE '{}' AS m".format(master), verbose=False)
            self.modify("ATTACH DATABASE '{}' AS c".format(conflicted), verbose=False)
            con.modify("ATTACH DATABASE '{}' AS c".format(conflicted), verbose=False)
            self.modify("ATTACH DATABASE '{}' AS m".format(master), verbose=False)

            # Drop any backup tables from failed merges
            for table in tables: self.modify("DROP TABLE IF EXISTS Backup_{0}".format(table), verbose=False)

            # Gather user data to add to CHANGELOG table
            import socket, datetime
            if not diff_only: user = get_input('Please enter your name : ')
            machine_name = socket.gethostname()
            date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
            modified_tables = []

            # Merge table by table, starting with SOURCES
            if not isinstance(tables, type(list())):
                tables = [tables]

            tables = tables or ['sources'] + [t for t in zip(*self.list(
                "SELECT * FROM sqlite_master WHERE name NOT LIKE '%Backup%' AND name!='sqlite_sequence' AND type='table'{}".format(
                    " AND name IN ({})".format("'" + "','".join(tables) + "'") if tables else '')))[1] if
                                              t != 'sources']
            for table in tables:
                # Get column names and data types from master table and column names from conflicted table
                metadata = self.query("PRAGMA table_info({})".format(table), fmt='table')
                columns, types, constraints = [np.array(metadata[n]) for n in ['name', 'type', 'notnull']]
                # columns, types, constraints = self.query("PRAGMA table_info({})".format(table), unpack=True)[1:4]
                conflicted_cols = con.query("PRAGMA table_info({})".format(table), unpack=True)[1]

                if any([i not in columns for i in conflicted_cols]):
                    # Abort table merge if conflicted has new columns not present in master. New columns must be added to the master database first via db.edit_columns().
                    print(
                    "\nMerge of {0} table aborted since conflicted copy has columns {1} not present in master.\nAdd new columns to master with astrodb.table() method and try again.\n".format(
                        table.upper(), [i for i in conflicted_cols if i not in columns]))

                else:
                    # Add new columns from master table to conflicted table if necessary
                    if any([i not in conflicted_cols for i in columns]):
                        con.modify("DROP TABLE IF EXISTS Conflicted_{0}".format(table))
                        con.modify("ALTER TABLE {0} RENAME TO Conflicted_{0}".format(table))
                        # TODO: Update to allow multiple primary and foreign keys
                        con.modify("CREATE TABLE {0} ({1})".format(table, ', '.join( \
                                ['{} {} {}{}'.format(c, t, r, ' UNIQUE PRIMARY KEY' if c == 'id' else '') \
                                 for c, t, r in zip(columns, types, constraints * ['NOT NULL'])])))
                        con.modify("INSERT INTO {0} ({1}) SELECT {1} FROM Conflicted_{0}".format(table, ','.join(
                            conflicted_cols)))
                        con.modify("DROP TABLE Conflicted_{0}".format(table))

                    # Pull unique records from conflicted table
                    data = map(list, con.list(
                        "SELECT * FROM (SELECT 1 AS db, {0} FROM m.{2} UNION ALL SELECT 2 AS db, {0} FROM c.{2}) GROUP BY {1} HAVING COUNT(*)=1 AND db=2".format(
                            ','.join(columns), ','.join(columns[1:]), table)).fetchall())

                    if data:

                        # Just print(the table differences
                        if diff_only:
                            pprint(zip(*data)[1:], names=columns, title='New {} records'.format(table.upper()))

                        # Add new records to the master and then clean up tables
                        else:
                            # Make temporary table copy so changes can be undone at any time
                            self.modify("DROP TABLE IF EXISTS Backup_{0}".format(table), verbose=False)
                            self.modify("ALTER TABLE {0} RENAME TO Backup_{0}".format(table), verbose=False)
                            self.modify("CREATE TABLE {0} ({1})".format(table, ', '.join( \
                                    ['{} {} {}{}'.format(c, t, r, ' UNIQUE PRIMARY KEY' if c == 'id' else '') \
                                     for c, t, r in zip(columns, types, constraints * ['NOT NULL'])])), verbose=False)
                            self.modify(
                                "INSERT INTO {0} ({1}) SELECT {1} FROM Backup_{0}".format(table, ','.join(columns)),
                                verbose=False)

                            # Create a dictionary of any reassigned ids from merged SOURCES tables and replace applicable source_ids in other tables.
                            print("\nMerging {} tables.\n".format(table.upper()))
                            try:
                                count = self.query("SELECT MAX(id) FROM {}".format(table), fetch='one')[0] + 1
                            except TypeError:
                                count = 1
                            for n, i in enumerate([d[1:] for d in data]):
                                if table == 'sources':
                                    reassign[i[0]] = count
                                elif 'source_id' in columns and i[1] in reassign.keys():
                                    i[1] = reassign[i[1]]
                                else:
                                    pass
                                i[0] = count
                                data[n] = i
                                count += 1

                            # Insert unique records into master
                            for d in data: self.modify(
                                "INSERT INTO {} VALUES({})".format(table, ','.join(['?' for c in columns])), d,
                                verbose=False)
                            pprint(zip(*data), names=columns,
                                   title="{} records added to {} table at '{}':".format(len(data), table, master))

                            # Run clean_up on the table to check for conflicts
                            abort = self.clean_up(table)

                            # Undo all changes to table if merge is aborted. Otherwise, push table changes to master.
                            if abort:
                                self.modify("DROP TABLE {0}".format(table), verbose=False)
                                self.modify("ALTER TABLE Backup_{0} RENAME TO {0}".format(table), verbose=False)
                            else:
                                self.modify("DROP TABLE Backup_{0}".format(table), verbose=False)
                                modified_tables.append(table.upper())

                    else:
                        print("\n{} tables identical.".format(table.upper()))

            # Add data to CHANGELOG table
            if not diff_only:
                user_description = get_input('\nPlease describe the changes made in this merge: ')
                self.list("INSERT INTO changelog VALUES(?, ?, ?, ?, ?, ?, ?)", \
                          (None, date, str(user), machine_name, ', '.join(modified_tables), user_description,
                           os.path.basename(conflicted)))

            # Finish up and detach
            if diff_only:
                print("\nDiff complete. No changes made to either database. Set `diff_only=False' to apply merge.")
            else:
                print("\nMerge complete!")

            con.modify("DETACH DATABASE c", verbose=False)
            self.modify("DETACH DATABASE c", verbose=False)
            con.modify("DETACH DATABASE m", verbose=False)
            self.modify("DETACH DATABASE m", verbose=False)
        else:
            print("File '{}' not found!".format(conflicted))

    def modify(self, SQL, params='', verbose=True):
        """
        Wrapper for CRUD operations to make them distinct from queries and automatically pass commit() method to cursor.

        Parameters
        ----------
        SQL: str
            The SQL query to execute
        params: sequence
            Mimics the native parameter substitution of sqlite3
        verbose: bool
                Prints the number of modified records
        """
        # Make sure the database isn't locked
        self.conn.commit()

        if SQL.lower().startswith('select'):
            print('Use self.query method for queries.')
        else:
            self.list(SQL, params)
            self.conn.commit()
            if verbose:
                print('Number of records modified: {}'.format(self.list("SELECT changes()").fetchone()[0] or '0'))

    def output_spectrum(self, spectrum, filepath, header={}):
        """
        Prints a file of the given spectrum to an ascii file with specified filepath.

        Parameters
        ----------
        spectrum: int, sequence
            The id from the SPECTRA table or a [w,f,e] sequence
        filepath: str
            The path of the file to print the data to.
        header: dict
                A dictionary of metadata to add of update in the header

        """
        # If an integer is supplied, get the spectrum from the SPECTRA table
        if isinstance(spectrum, int):
            data = self.query("SELECT * FROM spectra WHERE id={}".format(spectrum), fetch='one', fmt='dict')
            try:
                data['header'] = list(map(list, data['spectrum'].header.cards)) + [[k, v, ''] for k, v in
                                                                                   header.items()]
            except:
                data['header'] = ''

        # If a [w,f,e] sequence is supplied, make it into a Spectrum object
        elif isinstance(spectrum, (list, tuple, np.ndarray)):
            data = {'spectrum': Spectrum(spectrum, header=header), 'wavelength_units': '', 'flux_units': ''}
            try:
                data['header'] = list(map(list, data['spectrum'].header.cards))
            except:
                data['header'] = ''

        if data:
            fn = filepath if filepath.endswith('.txt') else filepath + 'spectrum.txt'

            # Write the header
            if data['header']:
                for n, line in enumerate(data['header']):
                    data['header'][n] = ['# {}'.format(str(line[0])).ljust(10)[:10],
                                         '{:50s} / {}'.format(*map(str, line[1:]))]
                try:
                    ii.write([np.asarray(i) for i in np.asarray(data['header']).T], fn, delimiter='\t',
                             format='no_header')
                except IOError:
                    pass

            # Write the data
            names = ['# wavelength [{}]'.format(data['wavelength_units']), 'flux [{}]'.format(data['flux_units'])]
            if len(data['spectrum'].data) == 3:
                if type(data['spectrum'].data[2]) in [np.ndarray, list]:
                    names += ['unc [{}]'.format(data['flux_units'])]
                else:
                    data['spectrum'].data = data['spectrum'].data[:2]

            with open(fn, mode='a') as f:
                ii.write([np.asarray(i, dtype=np.float64) for i in data['spectrum'].data], f, names=names,
                         delimiter='\t')

        else:
            print("Could not output spectrum: {}".format(spectrum))

    def plot_spectrum(self, spectrum_id, table='spectra', column='spectrum', overplot=False, color='b', norm=False):
        """
        Plots a spectrum from the given column and table

        Parameters
        ----------
        spectrum_id: int
            The id from the table of the spectrum to plot.
        overplot: bool
            Overplot the spectrum
        table: str
            The table from which the plot is being made
        column: str
            The column with SPECTRUM data type to plot
        color: str
            The color used for the data
        norm: bool, sequence
                True or (min,max) wavelength range in which to normalize the spectrum

        """
        i = self.query("SELECT * FROM {} WHERE id={}".format(table, spectrum_id), fetch='one', fmt='dict')
        if i:
            try:
                spec = scrub(i[column].data, units=False)
                w, f = spec[:2]
                try:
                    e = spec[2]
                except:
                    e = ''

                # Draw the axes and add the metadata
                if not overplot:
                    fig, ax = plt.subplots()
                    plt.rc('text', usetex=False)
                    ax.set_yscale('log', nonposy='clip')
                    plt.figtext(0.15, 0.88, '\n'.join(['{}: {}'.format(k, v) for k, v in i.items() if k != column]), \
                                verticalalignment='top')
                    try:
                        ax.set_xlabel(r'$\lambda$ [{}]'.format(i.get('wavelength_units')))
                        ax.set_ylabel(r'$F_\lambda$ [{}]'.format(i.get('flux_units')))
                    except:
                        pass
                    ax.legend(loc=8, frameon=False)
                else:
                    ax = plt.gca()

                # Normalize the data
                if norm:
                    try:
                        if isinstance(norm, bool): norm = (min(w), max(w))

                        # Normalize to the specified window
                        norm_mask = np.logical_and(w >= norm[0], w <= norm[1])
                        C = 1. / np.trapz(f[norm_mask], x=w[norm_mask])
                        f *= C
                        try:
                            e *= C
                        except:
                            pass

                    except:
                        print('Could not normalize.')

                # Plot the data
                ax.loglog(w, f, c=color, label='spec_id: {}'.format(i['id']))
                X, Y = plt.xlim(), plt.ylim()
                try:
                    ax.fill_between(w, f - e, f + e, color=color, alpha=0.3), ax.set_xlim(X), ax.set_ylim(Y)
                except:
                    print('No uncertainty array for spectrum {}'.format(spectrum_id))
                plt.ion()

            except IOError:
                print("Could not plot spectrum {}".format(spectrum_id))
                plt.close()

        else:
            print("No spectrum {} in the SPECTRA table.".format(spectrum_id))

    def query(self, SQL, params='', fmt='array', fetch='all', unpack=False, export='', \
              verbose=False, use_converters=True):
        """
        Wrapper for cursors so data can be retrieved as a list or dictionary from same method.

        Parameters
        ----------
        SQL: str
            The SQL query to execute
        params: sequence
            Mimics the native parameter substitution of sqlite3
        fmt: str
            Returns the data as a dictionary, array, astropy.table, or pandas.Dataframe
            given 'dict', 'array', 'table', or 'pandas'
        fetch: str
            String indicating whether to return **all** results or just **one**
        unpack: bool
            Returns the transpose of the data
        export: str
            The file path of the ascii file to which the data should be exported
        verbose: bool
            print the data as well
        use_converters: bool
            Apply converters to columns with custom data types

        Returns
        -------
        result: (array,dict,table)
            The result of the database query
        """
        try:
            # Restrict queries to SELECT and PRAGMA statements
            if SQL.lower().startswith('select') or SQL.lower().startswith('pragma'):

                # Make the query explicit so that column and table names are preserved
                SQL, columns = self._explicit_query(SQL, use_converters=use_converters)

                # Get the data as a dictionary
                dictionary = self.dict(SQL, params).fetchall()

                if any(dictionary):

                    # Fetch one
                    if fetch == 'one': dictionary = [dictionary.pop(0)]

                    # Make an Astropy table
                    table = at.Table(dictionary)

                    # Reorder the columns
                    try:
                        table = table[columns]
                    except:
                        pass

                    # Make an array
                    array = np.asarray(table)

                    # Unpack the results if necessary (data types are not preserved)
                    if unpack: array = np.array(zip(*array))

                    # print on screen
                    if verbose: pprint(table)

                    # print the results to file
                    if export:
                        # If .vot or .xml, assume VOTable export with votools
                        if export.lower().endswith('.xml') or export.lower().endswith('.vot'):
                            votools.dict_tovot(dictionary, export)
                        # Otherwise print as ascii
                        else:
                            ii.write(table, export, Writer=ii.FixedWidthTwoLine, fill_values=[('None', '-')])

                    # Or return the results
                    else:
                        if fetch == 'one':
                            dictionary, array = dictionary[0], array if unpack else np.array(list(array[0]))

                        if fmt == 'table':
                            return table
                        elif fmt == 'dict':
                            return dictionary
                        elif fmt == 'pandas':
                            return table.to_pandas()
                        else:
                            return array

                else:
                    return

            else:
                print(
                'Queries must begin with a SELECT or PRAGMA statement. For database modifications use self.modify() method.')

        except IOError:
            print('Could not execute: ' + SQL)

    def references(self, criteria, publications='publications', column_name='publication_shortname', fetch=False):
        """
        Do a reverse lookup on the **publications** table. Will return every entry that matches that reference.

        Parameters
        ----------
        criteria: int or str
            The id from the PUBLICATIONS table whose data across all tables is to be printed.
        publications: str
            Name of the publications table
        column_name: str
            Name of the reference column in other tables
        fetch: bool
            Return the results.

        Returns
        -------
        data_tables: dict
             Returns a dictionary of astropy tables with the table name as the keys.

        """

        data_tables = dict()

        # If an ID is provided but the column name is publication shortname, grab the shortname
        if isinstance(criteria, type(1)) and column_name == 'publication_shortname':
            t = self.query("SELECT * FROM {} WHERE id={}".format(publications, criteria), fmt='table')
            if len(t) > 0:
                criteria = t['shortname'][0]
            else:
                print('No match found for {}'.format(criteria))
                return

        t = self.query("SELECT * FROM sqlite_master WHERE type='table'", fmt='table')
        all_tables = t['name'].tolist()
        for table in ['sources'] + [t for t in all_tables if
                                    t not in ['publications', 'sqlite_sequence', 'sources']]:

            # Get the columns, pull out redundant ones, and query the table for this source's data
            t = self.query("PRAGMA table_info({})".format(table), fmt='table')
            columns = np.array(t['name'])
            types = np.array(t['type'])

            # Only get simple data types and exclude redundant ones for nicer printing
            columns = columns[
                ((types == 'REAL') | (types == 'INTEGER') | (types == 'TEXT')) & (columns != column_name)]

            # Query the table
            try:
                data = self.query("SELECT {} FROM {} WHERE {}='{}'".format(','.join(columns), table,
                                                                         column_name, criteria), fmt='table')
            except:
                data = None

            # If there's data for this table, save it
            if data:
                if fetch:
                    data_tables[table] = self.query(
                        "SELECT {} FROM {} WHERE {}='{}'".format(
                            ','.join(columns), table, column_name, criteria), fmt='table', fetch=True)
                else:
                    data = data[list(columns)]
                    pprint(data, title=table.upper())

        if fetch: return data_tables

    def save(self, directory='tabledata'):
        """
        Dump the entire contents of the database into a folder **directory** as ascii files

        Parameters
        ==========
        directory: str
            Directory name to store individual table data
        """
        from subprocess import call

        # Create the .sql file is it doesn't exist, i.e. if the Database class called a .db file initially
        if not os.path.isfile(self.sqlpath):
            self.sqlpath = self.dbpath.replace('.db', '.sql')
            os.system('touch {}'.format(self.sqlpath))

        # # Write the data to the .sql file
        # with open(self.sqlpath, 'w') as f:
        #     for line in self.conn.iterdump():
        #         f.write('%s\n' % line)

        # Alternatively...
        # Write the schema
        os.system("echo '.output {}\n.schema' | sqlite3 {}".format(self.sqlpath, self.dbpath))

        # Write the table files to the tabledata directory
        os.system("mkdir -p {}".format(directory))
        tables = self.query("select tbl_name from sqlite_master where type='table'")['tbl_name']
        tablepaths = [self.sqlpath]
        for table in tables:
            print('Generating {}...'.format(table))
            tablepath = '{0}/{1}.sql'.format(directory, table)
            tablepaths.append(tablepath)
            with open(tablepath, 'w') as f:
                for line in self.conn.iterdump():
                    if sys.version_info.major == 2:
                        # line = line.decode('utf-8')
                        line = line.encode('utf-8').decode('utf-8')

                    line = line.strip()
                    if line.startswith('INSERT INTO "{}"'.format(table)):
                        f.write('%s\n' % line.encode('ascii', 'ignore'))

        print("Tables saved to directory {}/".format(directory))
        print("""=======================================================================================
You can now run git to commit and push these changes, if needed.
For example, if on the master branch you can do the following:
  git add {0} {1}/*.sql
  git commit -m "COMMIT MESSAGE HERE"
  git push origin master
You can then issue a pull request on GitHub to have these changes reviewed and accepted
======================================================================================="""
              .format(self.sqlpath, directory))

        # Collect name and commit message from the user and push to Github
        # if git:
        #     user = get_input('Please enter your name : ')
        #     commit = get_input('Briefly describe the changes you have made : ')
        #     if user and commit:
        #         try:  # If not on the same branch, changes are overwritten or aborted: BAD
        #             call('git checkout {}'.format(branch), shell=True)
        #             call('git pull origin {}'.format(branch), shell=True)
        #             call('git add {}'.format(' '.join(tablepaths)), shell=True)
        #             call('git commit -m "(via astrodbkit) {}"'.format(commit), shell=True)
        #             call('git push origin {}'.format(branch), shell=True)
        #         except:
        #             print('Changes written to file but not pushed to Github. Make sure you are working from a git repo.')
        #     else:
        #         print('Sorry, astrodbkit needs a username and commit message to push changes to Guthub.')

    def schema(self, table):
        """
        Print the table schema

        Parameters
        ----------
        table: str
          The table name

        """
        try:
            pprint(self.query("PRAGMA table_info({})".format(table), fmt='table'))
        except ValueError:
            print('Table {} not found'.format(table))

    def search(self, criterion, table, columns='', fetch=False, radius=1/60.):
        """
        General search method for tables. For (ra,dec) input in decimal degrees,
        i.e. (12.3456,-65.4321), returns all sources within 1 arcminute, or the specified radius.
        For string input, i.e. 'vb10', returns all sources with case-insensitive partial text
        matches in columns with 'TEXT' data type. For integer input, i.e. 123, returns all
        exact matches of columns with INTEGER data type.

        Parameters
        ----------
        criterion: (str, int, sequence, tuple)
            The text, integer, coordinate tuple, or sequence thereof to search the table with.
        table: str
            The name of the table to search
        columns: sequence
            Specific column names to search, otherwise searches all columns
        fetch: bool
            Return the results of the query as an Astropy table
        radius: float
            Radius in degrees in which to search for objects if using (ra,dec). Default: 1/60 degree

        """

        # Get list of columns to search and format properly
        t = self.query("PRAGMA table_info({})".format(table), unpack=True, fmt='table')
        all_columns = t['name'].tolist()
        types = t['type'].tolist()
        columns = columns or all_columns
        columns = np.asarray([columns] if isinstance(columns, str) else columns)

        # Separate good and bad columns and corresponding types
        badcols = columns[~np.in1d(columns, all_columns)]
        columns = columns[np.in1d(columns, all_columns)]
        columns = np.array([c for c in all_columns if c in columns])
        types = np.array([t for c, t in zip(all_columns, types) if c in columns])[np.in1d(columns, all_columns)]
        for col in badcols:
            print("'{}' is not a column in the {} table.".format(col, table.upper()))

        # Coordinate search
        if sys.version_info[0] == 2:
            str_check = (str, unicode)
        else:
            str_check = str

        results = ''

        if isinstance(criterion, (tuple, list, np.ndarray)):
            try:
                t = self.query('SELECT id,ra,dec FROM sources', fmt='table')
                df = t.to_pandas()
                df[['ra', 'dec']] = df[['ra', 'dec']].apply(pd.to_numeric)  # convert everything to floats
                mask = df['ra'].isnull()
                df = df[~mask]

                def ang_sep(row, ra1, dec1):
                    # Using Vicenty Formula (http://en.wikipedia.org/wiki/Great-circle_distance)
                    # and adapting from astropy's SkyCoord

                    factor = math.pi / 180
                    sdlon = math.sin((row['ra'] - ra1) * factor)  # RA is longitude
                    cdlon = math.cos((row['ra'] - ra1) * factor)
                    slat1 = math.sin(dec1 * factor)  # Dec is latitude
                    slat2 = math.sin(row['dec'] * factor)
                    clat1 = math.cos(dec1 * factor)
                    clat2 = math.cos(row['dec'] * factor)

                    num1 = clat2 * sdlon
                    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
                    numerator = math.sqrt(num1 ** 2 + num2 ** 2)
                    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

                    return np.arctan2(numerator, denominator) / factor

                df['theta'] = df.apply(ang_sep, axis=1, args=(criterion[0], criterion[1]))
                good = df['theta'] <= radius

                if sum(good) > 0:
                    params = ", ".join(['{}'.format(s) for s in df[good]['id'].tolist()])
                    results = self.query('SELECT * FROM sources WHERE id IN ({})'.format(params), fmt='table')
            except:
                print("Could not search {} table by coordinates {}. Try again.".format(table.upper(), criterion))

        # Text string search of columns with 'TEXT' data type
        elif isinstance(criterion, str_check) and any(columns) and 'TEXT' in types:
            try:
                q = "SELECT * FROM {} WHERE {}".format(table, ' OR '.join([r"REPLACE(" + c + r",' ','') like '%" \
                     + criterion.replace(' ', '') + r"%'" for c, t in zip(columns,types[np.in1d(columns, all_columns)]) \
                     if t == 'TEXT']))
                results = self.query(q, fmt='table')
            except:
                print("Could not search {} table by string {}. Try again.".format(table.upper(), criterion))

        # Integer search of columns with 'INTEGER' data type
        elif isinstance(criterion, int):
            try:
                q = "SELECT * FROM {} WHERE {}".format(table, ' OR '.join(['{}={}'.format(c, criterion) \
                     for c, t in zip(columns, types[np.in1d(columns, all_columns)]) if t == 'INTEGER']))
                results = self.query(q, fmt='table')
            except:
                print("Could not search {} table by id {}. Try again.".format(table.upper(), criterion))

        # Problem!
        else:
            print("Could not search {} table by '{}'. Try again.".format(table.upper(), criterion))

        # Print or return the results
        if fetch:
            return results or at.Table(names=columns, dtype=[type_dict[t] for t in types], masked=True)
        else:
            if results: 
                pprint(results, title=table.upper())
            else:
                print("No results found for {} in the {} table.".format(criterion, table.upper()))

    def snapshot(self, name_db='export.db', version=1.0):
        """
        Function to generate a snapshot of the database by version number.

        Parameters
        ----------
        name_db: string
            Name of the new database (Default: export.db)
        version: float
            Version number to export (Default: 1.0)
        """

        # Check if file exists
        if os.path.isfile(name_db):
            import datetime
            date = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M")
            print("Renaming existing file {} to {}".format(name_db, name_db.replace('.db', date + '.db')))
            os.system("mv {} {}".format(name_db, name_db.replace('.db', date + '.db')))

        # Create a new database from existing database schema
        t, = self.query("select sql from sqlite_master where type = 'table'", unpack=True)
        schema = ';\n'.join(t) + ';'
        os.system("sqlite3 {} '{}'".format(name_db, schema))

        # Attach database to newly created database
        db = Database(name_db)
        db.list('PRAGMA foreign_keys=OFF')  # Temporarily deactivate foreign keys for the snapshot
        db.list("ATTACH DATABASE '{}' AS orig".format(self.dbpath))

        # For each table in database, insert records if they match version number
        t = db.query("SELECT * FROM sqlite_master WHERE type='table'", fmt='table')
        all_tables = t['name'].tolist()
        for table in [t for t in all_tables if t not in ['sqlite_sequence', 'changelog']]:
            # Check if this is table with version column
            metadata = db.query("PRAGMA table_info({})".format(table), fmt='table')
            columns, types, required, pk = [np.array(metadata[n]) for n in ['name', 'type', 'notnull', 'pk']]
            print(table, columns)
            if 'version' not in columns:
                db.modify("INSERT INTO {0} SELECT * FROM orig.{0}".format(table))
            else:
                db.modify("INSERT INTO {0} SELECT * FROM orig.{0} WHERE orig.{0}.version<={1}".format(table, version))

        # Detach original database
        db.list('DETACH DATABASE orig')
        db.list('PRAGMA foreign_keys=ON')
        db.close()

    def table(self, table, columns, types, constraints='', pk='', new_table=False):
        """
        Rearrange, add or delete columns from database **table** with desired ordered list of **columns** and corresponding data **types**.

        Parameters
        ----------
        table: sequence
            The name of the table to modify
        columns: list
            A sequence of the columns in the order in which they are to appear in the SQL table
        types: sequence
            A sequence of the types corresponding to each column in the columns list above.
        constraints: sequence (optional)
            A sequence of the constraints for each column, e.g. '', 'UNIQUE', 'NOT NULL', etc.
        pk: string or list
            Name(s) of the primary key(s) if other than ID
        new_table: bool
            Create a new table

        """
        goodtogo = True

        # Make sure there is an integer primary key, unique, not null 'id' column
        # and the appropriate number of elements in each sequence
        if columns[0] != 'id':
            print("Column 1 must be called 'id'")
            goodtogo = False

        # if types[0].upper() != 'INTEGER PRIMARY KEY':
        #     print("'id' column type must be 'INTEGER PRIMARY KEY'")
        #     goodtogo = False

        if constraints:
            if 'UNIQUE' not in constraints[0].upper() and 'NOT NULL' not in constraints[0].upper():
                print("'id' column constraints must be 'UNIQUE NOT NULL'")
                goodtogo = False
        else:
            constraints = ['UNIQUE NOT NULL'] + ([''] * (len(columns) - 1))

        # Set UNIQUE NOT NULL constraints for the primary keys, except ID which is already has them
        if pk:
            if not isinstance(pk, type(list())):
                pk = list(pk)

            for elem in pk:
                if elem == 'id':
                    continue
                else:
                    ind, = np.where(columns == elem)
                    constraints[ind] = 'UNIQUE NOT NULL'
        else:
            pk = ['id']

        if not len(columns) == len(types) == len(constraints):
            print("Must provide equal length *columns ({}), *types ({}), and *constraints ({}) sequences." \
                  .format(len(columns), len(types), len(constraints)))
            goodtogo = False

        if goodtogo:
            t = self.query("SELECT name FROM sqlite_master", unpack=True, fmt='table')
            tables = t['name'].tolist()

            # If the table exists, modify the columns
            if table in tables and not new_table:

                # Rename the old table and create a new one
                self.list("DROP TABLE IF EXISTS TempOldTable")
                self.list("ALTER TABLE {0} RENAME TO TempOldTable".format(table))
                create_txt = "CREATE TABLE {0} ({1}".format(table, ', '.join(
                        ['{} {} {}'.format(c, t, r) for c, t, r in zip(columns, types, constraints)]))
                create_txt += ', PRIMARY KEY({}))'.format(', '.join([elem for elem in pk]))
                # print(create_txt.replace(',', ',\n'))
                self.list(create_txt)

                # Populate the new table and drop the old one
                old_columns = [c for c in self.query("PRAGMA table_info(TempOldTable)", unpack=True)[1] if c in columns]
                self.list("INSERT INTO {0} ({1}) SELECT {1} FROM TempOldTable".format(table, ','.join(old_columns)))

                # Check for and add any foreign key constraints
                t = self.query('PRAGMA foreign_key_list(TempOldTable)', fmt='table')
                if not isinstance(t, type(None)):
                    self.list("DROP TABLE TempOldTable")
                    self.add_foreign_key(table, t['table'].tolist(), t['from'].tolist(), t['to'].tolist())
                else:
                    self.list("DROP TABLE TempOldTable")

            # If the table does not exist and new_table is True, create it
            elif table not in tables and new_table:
                create_txt = "CREATE TABLE {0} ({1}".format(table, ', '.join(
                    ['{} {} {}'.format(c, t, r) for c, t, r in zip(columns, types, constraints)]))
                create_txt += ', PRIMARY KEY({}))'.format(', '.join([elem for elem in pk]))
                print(create_txt.replace(',', ',\n'))
                self.list(create_txt)

            # Otherwise the table to be modified doesn't exist or the new table to add already exists, so do nothing
            else:
                if new_table:
                    print('Table {} already exists. Set *new_table=False to modify.'.format(table.upper()))
                else:
                    print('Table {} does not exist. Could not modify. Set *new_table=True to add a new table.'.format(
                        table.upper()))

        else:
            print('The {} table has not been {}. Please make sure your table columns, \
             types, and constraints are formatted properly.'.format(table.upper(), \
                                                                    'created' if new_table else 'modified'))


class Spectrum:
    def __init__(self, data, header='', path=''):
        """
        Initialize the Spectrum object

        Parameters
        ----------
        data: sequence
            The [w,f,e] spectrum
        header: dictionary (optional)
            A dictionary of data to include in the header
        path: str (optional)
            The absolute path to the original file

        Returns
        -------
        object
            The Spectrum object

        """
        self.data = data
        self.path = path

        if header and isinstance(header, dict):
            new_header = pf.Header()
            for k, v in header.items():
                new_header[k.replace('.', '_').replace('#', '')] = v
        elif isinstance(header, pf.header.Header):
            new_header = header
        elif isinstance(header, list):
            new_header = pf.Header()
            for line in header:
                try:
                    k, v = line.split('=')[0], '\\'.join(line.split('=')[1:])
                    new_header[k.replace('\\', '').replace('.', '_').replace('#', '').strip()] = v
                except:
                    pass
        elif header:
            print('Header is {}. Must be a fits header, list, or dictionary.'.format(type(header)))
            new_header = ''
        else:
            new_header = ''

        self.header = new_header


# ==============================================================================================================================================
# ================================= Adapters and converters for special data types =============================================================
# ==============================================================================================================================================

def adapt_array(arr):
    """
    Adapts a Numpy array into an ARRAY string to put into the database.

    Parameters
    ----------
    arr: array
        The Numpy array to be adapted into an ARRAY type that can be inserted into a SQL file.

    Returns
    -------
    ARRAY
            The adapted array object

    """
    out = io.BytesIO()
    np.save(out, arr), out.seek(0)
    return buffer(out.read())  # TODO: Fix for Python 3


def convert_array(array):
    """
    Converts an ARRAY string stored in the database back into a Numpy array.

    Parameters
    ----------
    array: ARRAY
        The array object to be converted back into a Numpy array.

    Returns
    -------
    array
            The converted Numpy array.

    """
    out = io.BytesIO(array)
    out.seek(0)
    return np.load(out)


# TODO: Eliminate this, not being used anymore
def adapt_spectrum(spec):
    """
    Adapts a SPECTRUM object into a string to put into the database

    Parameters
    ----------
    spec: str, astrodbkit.astrodb.Spectrum
        The spectrum object to convert or string to put into the database

    Returns
    -------
    spec: str
            The file path to place in the database
    """
    if isinstance(spec, str):
        pass
    else:
        spec = '$BDNYC_spectra' + spec.path.split('BDNYC_spectra')[1]

    return spec


def convert_spectrum(File, verbose=False):
    """
    Converts a SPECTRUM data type stored in the database into a (W,F,E) sequence of arrays.

    Parameters
    ----------
    File: SPECTRUM
        The URL or filepath of the file to be converted into arrays.

    Returns
    -------
    sequence
            The converted spectrum.

    """
    spectrum, header = '', ''
    if isinstance(File, type(b'')):  # Decode if needed (ie, for Python 3)
        File = File.decode('utf-8')

    if isinstance(File, (str, type(u''))):

        # Convert variable path to absolute path
        if File.startswith('$'):
            abspath = os.popen('echo {}'.format(File.split('/')[0])).read()[:-1]
            if abspath: File = File.replace(File.split('/')[0], abspath)

        if File.startswith('http'):
            if verbose: print('Downloading {}'.format(File))
            downloaded_file = download_file(File, cache=True)  # download only once
        else:
            downloaded_file = File

        try:  # Try FITS files first
            # Get the data
            # try:
            spectrum, header = pf.getdata(downloaded_file, cache=True, header=True)
            # except:
            #     spectrum, header = pf.getdata(File, cache=False, header=True)

            # Check the key type
            KEY_TYPE = ['CTYPE1']
            setType = set(KEY_TYPE).intersection(set(header.keys()))
            if len(setType) == 0:
                isLinear = True
            else:
                valType = header[setType.pop()]
                isLinear = valType.strip().upper() == 'LINEAR'

            # Get wl, flux & error data from fits file
            spectrum = __get_spec(spectrum, header, File)

            # Generate wl axis when needed
            if not isinstance(spectrum[0],np.ndarray):
                spectrum[0] = __create_waxis(header, len(spectrum[1]), File)

            # If no wl axis generated, then clear out all retrieved data for object
            if not isinstance(spectrum[0],np.ndarray):
                spectrum = None

            if verbose: print('Read as FITS...')
        except (IOError, KeyError):
            # Check if the FITS file is just Numpy arrays
            try:
                spectrum, header = pf.getdata(downloaded_file, cache=True, header=True)
                if verbose: print('Read as FITS Numpy array...')
            except (IOError, KeyError):
                try: # Try ascii
                    spectrum = ii.read(downloaded_file)
                    spectrum = np.array([np.asarray(spectrum.columns[n]) for n in range(len(spectrum.columns))])
                    if verbose: print('Read as ascii...')

                    txt, header = open(downloaded_file), []
                    for i in txt:
                        if any([i.startswith(char) for char in ['#', '|', '\\']]):
                            header.append(i.replace('\n', ''))
                    txt.close()
                except:
                    pass

    if spectrum == '':
        print('Could not retrieve spectrum at {}.'.format(File))
        return File
    else:
        spectrum = Spectrum(spectrum, header, File)
        return spectrum


def __create_waxis(fitsHeader, lenData, fileName, wlog=False, verb=True):
    # Define key names in
    KEY_MIN = ['COEFF0', 'CRVAL1']  # Min wl
    KEY_DELT = ['COEFF1', 'CDELT1', 'CD1_1']  # Delta of wl
    KEY_OFF = ['LTV1']  # Offset in wl to subsection start

    # Find key names for minimum wl, delta, and wl offset in fits header
    setMin = set(KEY_MIN).intersection(set(fitsHeader.keys()))
    setDelt = set(KEY_DELT).intersection(set(fitsHeader.keys()))
    setOff = set(KEY_OFF).intersection(set(fitsHeader.keys()))

    # Get the values for minimum wl, delta, and wl offset, and generate axis
    if len(setMin) >= 1 and len(setDelt) >= 1:
        nameMin = setMin.pop()
        valMin = fitsHeader[nameMin]

        nameDelt = setDelt.pop()
        valDelt = fitsHeader[nameDelt]

        if len(setOff) == 0:
            valOff = 0
        else:
            nameOff = setOff.pop()
            valOff = fitsHeader[nameOff]

        # generate wl axis
        if nameMin == 'COEFF0' or wlog == True:
            # SDSS fits files
            wAxis = 10 ** (np.arange(lenData) * valDelt + valMin)
        else:
            wAxis = (np.arange(lenData) * valDelt) + valMin - (valOff * valDelt)

    else:
        wAxis = None
        if verb:
            print('Could not re-create wavelength axis for ' + fileName + '.')

    return wAxis


def __get_spec(fitsData, fitsHeader, fileName, verb=True):
    validData = [None] * 3

    # Identify number of data sets in fits file
    dimNum = len(fitsData)

    # Identify data sets in fits file
    fluxIdx = None
    waveIdx = None
    sigmaIdx = None

    if dimNum == 1:
        fluxIdx = 0
    elif dimNum == 2:
        if len(fitsData[0]) == 1:
            sampleData = fitsData[0][0][20]
        else:
            sampleData = fitsData[0][20]
        if sampleData < 0.0001:
            # 0-flux, 1-unknown
            fluxIdx = 0
        else:
            waveIdx = 0
            fluxIdx = 1
    elif dimNum == 3:
        waveIdx = 0
        fluxIdx = 1
        sigmaIdx = 2
    elif dimNum == 4:
        # 0-flux clean, 1-flux raw, 2-background, 3-sigma clean
        fluxIdx = 0
        sigmaIdx = 3
    elif dimNum == 5:
        # 0-flux, 1-continuum substracted flux, 2-sigma, 3-mask array, 4-unknown
        fluxIdx = 0
        sigmaIdx = 2
    elif dimNum > 10:
        # Implies that only one data set in fits file: flux
        fluxIdx = -1
        if np.isscalar(fitsData[0]):
            fluxIdx = -1
        elif len(fitsData[0]) == 2:
            # Data comes in a xxxx by 2 matrix (ascii origin)
            tmpWave = []
            tmpFlux = []
            for pair in fitsData:
                tmpWave.append(pair[0])
                tmpFlux.append(pair[1])
            fitsData = [tmpWave, tmpFlux]
            fitsData = np.array(fitsData)

            waveIdx = 0
            fluxIdx = 1
        else:
            # Indicates that data is structured in an unrecognized way
            fluxIdx = None
    else:
        fluxIdx = None

    # Fetch wave data set from fits file
    if fluxIdx is None:
        # No interpretation known for fits file data sets
        validData = None
        if verb:
            print('Unable to interpret data in ' + fileName + '.')
        return validData
    else:
        if waveIdx is not None:
            if len(fitsData[waveIdx]) == 1:
                # Data set may be a 1-item list
                validData[0] = fitsData[waveIdx][0]
            else:
                validData[0] = fitsData[waveIdx]

    # Fetch flux data set from fits file
    if fluxIdx == -1:
        validData[1] = fitsData
    else:
        if len(fitsData[fluxIdx]) == 1:
            validData[1] = fitsData[fluxIdx][0]
        else:
            validData[1] = fitsData[fluxIdx]

    # Fetch sigma data set from fits file, if requested
    if sigmaIdx is None:
        validData[2] = np.array([np.nan] * len(validData[1]))
    else:
        if len(fitsData[sigmaIdx]) == 1:
            validData[2] = fitsData[sigmaIdx][0]
        else:
            validData[2] = fitsData[sigmaIdx]

    # If all sigma values have the same value, replace them with nans
    if validData[2][10] == validData[2][11] == validData[2][12]:
        validData[2] = np.array([np.nan] * len(validData[1]))

    return validData


# Register the adapters
sqlite3.register_adapter(np.ndarray, adapt_array)
# sqlite3.register_adapter(str, adapt_spectrum)

# Register the converters
sqlite3.register_converter(str("ARRAY"), convert_array)
sqlite3.register_converter(str("SPECTRUM"), convert_spectrum)


def pprint(data, names='', title='', formats={}):
    """
    Prints tables with a bit of formatting

    Parameters
    ----------
    data: (sequence, dict, table)
        The data to print(in the table
    names: sequence
        The column names
    title: str (optional)
        The title of the table
    formats: dict
        A dictionary of column:format values

    """
    # Make the data into a table if it isn't already
    if type(data) != at.Table:
        data = at.Table(data, names=names)

    # Make a copy
    pdata = data.copy()

    # Put the title in the metadata
    try:
        title = title or pdata.meta['name']
    except:
        pass

    # Shorten the column names for slimmer data
    for old, new in zip(*[pdata.colnames, [
        i.replace('wavelength', 'wav').replace('publication', 'pub').replace('instrument', 'inst')\
        .replace('telescope','scope') for i in pdata.colnames]]):
        pdata.rename_column(old, new) if new != old else None

    # Format the columns
    formats.update({'comments': '%.15s', 'obs_date': '%.10s', 'names': '%.30s', 'description': '%.50s'})

    # print it!
    if title: print('\n' + title)
    try:
        ii.write(pdata, sys.stdout, Writer=ii.FixedWidthTwoLine, formats=formats, fill_values=[('None', '-')])
    except UnicodeDecodeError:  # Fix for Unicode characters. Print out in close approximation to ii.write()
        max_length = 50
        str_lengths = dict()
        for key in pdata.keys():
            lengths = map(lambda x: len(str(x).decode('utf-8')), pdata[key].data)
            lengths.append(len(key))
            str_lengths[key] = min(max(lengths), max_length)
        print(' '.join(key.rjust(str_lengths[key]) for key in pdata.keys()))
        print(' '.join('-' * str_lengths[key] for key in pdata.keys()))
        for i in pdata:
            print(' '.join([str(i[key]).decode('utf-8')[:max_length].rjust(str_lengths[key])
                           if i[key] else '-'.rjust(str_lengths[key]) for key in pdata.keys()]))


def clean_header(header):
    try:
        header = pf.open(File, ignore_missing_end=True)[0].header  # TODO: Fix missing File reference
        new_header = pf.Header()
        for x, y, z in header.cards: new_header[x.replace('.', '_').replace('#', '')] = (y, z)
        header = pf.PrimaryHDU(header=new_header).header
    except:
        pass

    return header


def _help():
    print(' ')
    command = '{:<33s}'.format('Command')
    result = '{:79s}'.format('Result')
    pprint(np.asarray([['<column name>', 'Display full record entry for that column without taking action'], \
                       ['k', 'Keeps both records and assigns second one new id if necessary'], \
                       ['r', 'Replaces all columns of first record with second record values and deletes second record'], \
                       ['r <column name> <column name> ...',
                        'Replaces specified columns of first record with second record values and deletes second record'], \
                       ['c', 'Complete empty columns of first record with second record values where possible and deletes second record'], \
                       ['[Enter]', 'Keep first record and delete second'], \
                       ['sql <SQLite query>', 'Execute arbitrary raw SQLite command'], \
                       ['abort', 'Abort merge of current table, undo all changes, and proceed to next table']]), \
           names=[command, result], \
           formats={command: '%-33s', result: '%-79s'})


def scrub(data, units=False):
    """
    For input data [w,f,e] or [w,f] returns the list with NaN, negative, and zero flux
    (and corresponding wavelengths and errors) removed.
    """
    units = [i.unit if hasattr(i, 'unit') else 1 for i in data]
    data = [np.asarray(i.value if hasattr(i, 'unit') else i, dtype=np.float32) for i in data if
            isinstance(i, np.ndarray)]
    data = [i[np.where(~np.isinf(data[1]))] for i in data]
    data = [i[np.where(np.logical_and(data[1] > 0, ~np.isnan(data[1])))] for i in data]
    data = [i[np.unique(data[0], return_index=True)[1]] for i in data]
    return [i[np.lexsort([data[0]])] * Q for i, Q in zip(data, units)] if units else [i[np.lexsort([data[0]])] for i in
                                                                                      data]


def _autofill_spec_record(record):
    """
    Returns an astropy table with columns auto-filled from FITS header

    Parameters
    ----------
    record: astropy.io.fits.table.table.Row
        The spectrum table row to scrape

    Returns
    -------
    record: astropy.io.fits.table.table.Row
        The spectrum table row with possible new rows inserted
    """
    try:
        record['filename'] = os.path.basename(record['spectrum'])

        if record['spectrum'].endswith('.fits'):
            header = pf.getheader(record['spectrum'])

            # Wavelength units
            if not record['wavelength_units']:
                try:
                    record['wavelength_units'] = header['XUNITS']
                except KeyError:
                    try:
                        if header['BUNIT']: record['wavelength_units'] = 'um'
                    except KeyError:
                        pass
            if 'microns' in record['wavelength_units'] or 'Microns' in record['wavelength_units'] or 'um' in record[
                'wavelength_units']:
                record['wavelength_units'] = 'um'

            # Flux units
            if not record['flux_units']:
                try:
                    record['flux_units'] = header['YUNITS'].replace(' ', '')
                except KeyError:
                    try:
                        record['flux_units'] = header['BUNIT'].replace(' ', '')
                    except KeyError:
                        pass
            if 'erg' in record['flux_units'] and 'A' in record['flux_units']:
                record['flux_units'] = 'ergs-1cm-2A-1' if 'erg' in record['flux_units'] and 'A' in record['flux_units'] \
                    else 'ergs-1cm-2um-1' if 'erg' in record['flux_units'] and 'um' in record['flux_units'] \
                    else 'Wm-2um-1' if 'W' in record['flux_units'] and 'um' in record['flux_units'] \
                    else 'Wm-2A-1' if 'W' in record['flux_units'] and 'A' in record['flux_units'] \
                    else ''

            # Observation date
            if not record['obs_date']:
                try:
                    record['obs_date'] = header['DATE_OBS']
                except KeyError:
                    try:
                        record['obs_date'] = header['DATE-OBS']
                    except KeyError:
                        try:
                            record['obs_date'] = header['DATE']
                        except KeyError:
                            pass

            # Telescope id
            if not record['telescope_id']:
                try:
                    n = header['TELESCOP'].lower() if isinstance(header['TELESCOP'], str) else ''
                    record['telescope_id'] = 5 if 'hst' in n \
                        else 6 if 'spitzer' in n \
                        else 7 if 'irtf' in n \
                        else 9 if 'keck' in n and 'ii' in n \
                        else 8 if 'keck' in n and 'i' in n \
                        else 10 if 'kp' in n and '4' in n \
                        else 11 if 'kp' in n and '2' in n \
                        else 12 if 'bok' in n \
                        else 13 if 'mmt' in n \
                        else 14 if 'ctio' in n and '1' in n \
                        else 15 if 'ctio' in n and '4' in n \
                        else 16 if 'gemini' in n and 'north' in n \
                        else 17 if 'gemini' in n and 'south' in n \
                        else 18 if ('vlt' in n and 'U2' in n) \
                        else 19 if '3.5m' in n \
                        else 20 if 'subaru' in n \
                        else 21 if ('mag' in n and 'ii' in n) or ('clay' in n) \
                        else 22 if ('mag' in n and 'i' in n) or ('baade' in n) \
                        else 23 if ('eso' in n and '1m' in n) \
                        else 24 if 'cfht' in n \
                        else 25 if 'ntt' in n \
                        else 26 if ('palomar' in n and '200-inch' in n) \
                        else 27 if 'pan-starrs' in n \
                        else 28 if ('palomar' in n and '60-inch' in n) \
                        else 29 if ('ctio' in n and '0.9m' in n) \
                        else 30 if 'soar' in n \
                        else 31 if ('vlt' in n and 'U3' in n) \
                        else 32 if ('vlt' in n and 'U4' in n) \
                        else 33 if 'gtc' in n \
                        else None
                except KeyError:
                    pass

            # Instrument id
            if not record['instrument_id']:
                try:
                    i = header['INSTRUME'].lower()
                    record[
                        'instrument_id'] = 1 if 'r-c spec' in i or 'test' in i or 'nod' in i else 2 if 'gmos-n' in i else 3 if 'gmos-s' in i else 4 if 'fors' in i else 5 if 'lris' in i else 6 if 'spex' in i else 7 if 'ldss3' in i else 8 if 'focas' in i else 9 if 'nirspec' in i else 10 if 'irs' in i else 11 if 'fire' in i else 12 if 'mage' in i else 13 if 'goldcam' in i else 14 if 'sinfoni' in i else 15 if 'osiris' in i else 16 if 'triplespec' in i else 17 if 'x-shooter' in i else 18 if 'gnirs' in i else 19 if 'wircam' in i else 20 if 'cormass' in i else 21 if 'isaac' in i else 22 if 'irac' in i else 23 if 'dis' in i else 24 if 'susi2' in i else 25 if 'ircs' in i else 26 if 'nirc' in i else 29 if 'stis' in i else 0
                except KeyError:
                    pass

    except:
        pass

    return record


type_dict = {'INTEGER': np.dtype('int64'), 'REAL': np.dtype('float64'), 'TEXT': np.dtype('object'),
             'ARRAY': np.dtype('object'), 'SPECTRUM': np.dtype('object'), 'BOOLEAN': np.dtype('bool')}
