#!/usr/bin/python
# encoding: utf-8
# Author: Joe Filippazzo, jcfilippazzo@gmail.com

import io, os, sys, itertools, sqlite3, warnings
import numpy as np, matplotlib.pyplot as plt
import astropy.io.fits as pf, astropy.io.ascii as ii, astropy.table as at
import votools
warnings.simplefilter('ignore')

def create_database(dbpath):
  """
  Create a new database at the given dbpath
  
  Parameters
  ----------
  dbpath: str
    The full path for the new database, including the filename and .db file extension.

  """
  if dbpath.endswith('.db'):
    sources_table = "CREATE TABLE sources (id INTEGER PRIMARY KEY, ra REAL, dec REAL, designation TEXT, publication_id INTEGER, shortname TEXT, names TEXT, comments TEXT)"
    os.system("sqlite3 {} '{}'".format(dbpath,sources_table))
    if os.path.isfile(dbpath):
      print "\nDatabase created! To load, run\n\ndb = astrodb.get_db('{}')\n\nThen run db.modify_table() method to create tables.".format(dbpath)
  else: print "Please provide a path and file name with a .db file extension, e.g. /Users/<username>/Desktop/test.db"

class get_db:
  def __init__(self, dbpath):
    """
    Initialize the database.
    
    Parameters
    ----------
    dbpath: str 
      The path to the database file. 
    
    Returns
    -------
    object
      The database object
         
    """
    
    if os.path.isfile(dbpath):
    
      # Create connection
      con = sqlite3.connect(dbpath, isolation_level=None, detect_types=sqlite3.PARSE_DECLTYPES)
      con.text_factory = str
      self.conn = con
      self.list = con.cursor().execute
            
      # Make dictionary
      def dict_factory(cursor, row):
        d = {}
        for idx,col in enumerate(cursor.description): d[col[0]] = row[idx]
        return d
        
      self.dict = con.cursor()
      self.dict.row_factory = dict_factory
      self.dict = self.dict.execute
    
      # Make sure the ignore table exists
      self.list("CREATE TABLE IF NOT EXISTS ignore (id INTEGER PRIMARY KEY, id1 INTEGER, id2 INTEGER, tablename TEXT)")
    
    else: print "Sorry, no such file '{}'".format(dbpath)

  def add_data(self, data, table, delimiter='|', bands=''):
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
    
    """ 
    # Store raw entry
    entry = data
    
    # Digest the ascii file into table
    if isinstance(data,str) and os.path.isfile(data): 
      data = ii.read(data)
    
    # Or add the sequence of data elements into a table
    elif isinstance(data,(list,tuple,np.ndarray)):
      data = at.Table(list(np.asarray(data[1:]).T), names=data[0], dtype=[type(i) for i in data[1]])
    
    else: data = None
    
    if data:

      # Get list of all columns and make an empty table for new records
      metadata = self.query("PRAGMA table_info({})".format(table), fmt='table')
      columns, types, required = [np.array(metadata[n]) for n in ['name','type','notnull']]
      new_records = at.Table(names=columns, dtype=[type_dict[t] for t in types])
    
      # If a row contains photometry for multiple bands, use the *multiband argument and execute this
      if bands and table.lower()=='photometry':    
      
        # Pull out columns that are band names
        for b in list(set(bands)&set(data.colnames)):
          try:
            # Get the repeated data plus the band data and rename the columns
            band = data[list(set(columns)&set(data.colnames))+[b,b+'_unc']]
            for suf in ['','_unc']: band.rename_column(b+suf,'magnitude'+suf)
            band.add_column(at.Column([b]*len(band), name='band'))

            # Add the band data to the list of new_records
            new_records = at.vstack([new_records,band])
          except IOError: pass
            
      else:      
        # Inject data into full database table format
        new_records = at.vstack([new_records,data])[new_records.colnames]
    
      # Reject rows that fail column requirements, e.g. NOT NULL fields like 'source_id'
      for r in columns[np.where(np.logical_and(required,columns!='id'))]: 
        new_records = new_records[np.where(new_records[r])]
    
      # For spectra, try to populate the table by reading the FITS header
      if table.lower()=='spectra':
        del_records = []
        for n,new_rec in enumerate(new_records):
          
          # Convert relative path to absolute path
          relpath = new_rec['spectrum']
          if relpath.startswith('$'):
            abspath = os.popen('echo {}'.format(relpath.split('/')[0])).read()[:-1]
            if abspath: new_rec['spectrum'] = relpath.replace(relpath.split('/')[0],abspath)
          
          # Test if the file exists and try to pull metadata from the FITS header
          if os.path.isfile(new_rec['spectrum']):
            new_records[n]['spectrum'] = relpath
            new_records[n] = _autofill_spec_record(new_rec)
          else: 
            print 'Error adding the spectrum at {}'.format(new_rec['spectrum'])
            del_records.append(n)
        
        # Remove bad records from the table
        new_records.remove_rows(del_records)

      # Add the new records
      for new_rec in new_records:
        new_rec = list(new_rec)
        for n,col in enumerate(new_rec): 
          if type(col)==np.ma.core.MaskedConstant: new_rec[n] = None
        self.modify("INSERT INTO {} VALUES({})".format(table, ','.join('?'*len(columns))), new_rec)
    
      # Print a table of the new records or bad news
      if new_records: 
        pprint(new_records, names=columns, title="{} new records added to the {} table.".format(len(new_records),table.upper()))
      else: 
        print 'No new records added to the {} table. Please check your input: {}'.format(table,entry)
    
      # Run table clean up
      try: self.clean_up(table)
      except: print 'Could not run clean_up() method.' 
    
    else: print 'Please check your input: {}'.format(entry)
 
  def clean_up(self, table):
    """
    Removes exact duplicates, blank records or data without a *source_id* from the specified **table**. 
    Then finds possible duplicates and prompts for conflict resolution.
    
    Parameters
    ----------
    table: str
      The name of the table to remove duplicates, blanks, and data without source attributions.
    
    """
    # Get the table info and all the records
    metadata = self.query("PRAGMA table_info({})".format(table), fmt='table')
    columns, types, required = [np.array(metadata[n]) for n in ['name','type','notnull']]
    records = self.query("SELECT * FROM {}".format(table), fmt='table')
    ignore = self.query("SELECT * FROM ignore WHERE tablename LIKE ?", (table,))
    duplicate, command = [1], ''

    # Remove records with missing required values
    req_keys = columns[np.where(required)]
    try:
      self.modify("DELETE FROM {} WHERE {}".format(table, ' OR '.join([i+' IS NULL' for i in req_keys])), verbose=False)
      self.modify("DELETE FROM {} WHERE {}".format(table, ' OR '.join([i+" IN ('null','None','')" for i in req_keys])), verbose=False)
    except: pass
    
    # Remove exact duplicates
    self.modify("DELETE FROM {0} WHERE id NOT IN (SELECT min(id) FROM {0} GROUP BY {1})".format(table,', '.join(columns[1:])), verbose=False)

    # Check for records with identical required values but different ids.            
    if table.lower()!='sources': req_keys = columns[np.where(np.logical_and(required,columns!='id'))]

    # List of old and new pairs to ignore
    ignore, new_ignore = ignore or [], []

    while any(duplicate):      
      # Pull out duplicates one by one
      duplicate = self.query("SELECT t1.id, t2.id FROM {0} t1 JOIN {0} t2 ON t1.source_id=t2.source_id WHERE t1.id!=t2.id AND {1}{2}{3}"\
                              .format(table, ' AND '.join(['t1.{0}=t2.{0}'.format(i) for i in req_keys]), (' AND '\
                              +' AND '.join(["(t1.id NOT IN ({0}) and t2.id NOT IN ({0}))".format(','.join(map(str,[id1,id2]))) for id1,id2 \
                              in zip(ignore['id1'],ignore['id2'])])) if ignore else '', (' AND '\
                              +' AND '.join(["(t1.id NOT IN ({0}) and t2.id NOT IN ({0}))".format(','.join(map(str,ni))) for ni \
                              in new_ignore])) if new_ignore else ''), fetch='one')

      # Compare potential duplicates and prompt user for action on each
      try:
        
        # Run record matches through comparison and return the command 
        command = self._compare_records(table, duplicate, delete=True)
        
        # Add acceptible duplicates to ignore list or abort
        if command=='keep': 
          new_ignore.append([duplicate[0],duplicate[1]])
          self.list("INSERT INTO ignore VALUES(?,?,?,?)", (None,duplicate[0],duplicate[1],table.lower()))
        elif command=='undo': pass # Add this functionality!
        elif command=='abort': break
        else: pass
      except: break

    # Finish or abort table clean up
    if command=='abort': 
      print '\nAborted clean up of {} table.'.format(table.upper())
      return 'abort'
    else: print '\nFinished clean up on {} table.'.format(table.upper())

  def _compare_records(self, table, duplicate, options=['r','c','k','sql'], delete=False):
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
    delete: bool
      Delete the record with the higher id.
    
    """
    # Print the old and new records suspectred of being duplicates
    data = self.query("SELECT * FROM {} WHERE id IN ({})".format(table,','.join(map(str,duplicate))), fmt='table', verbose=True)
    columns = data.colnames[1:]
    old, new = [[data[n][k] for k in columns[1:]] for n in [0,1]]

    # Print the command key
    replace = raw_input("\nKeep both records [k]? Or replace [r], complete [c], or keep only [Press *Enter*] record {}? (Type column name to inspect or 'help' for options): ".format(duplicate[0])).lower()
    while replace in columns or replace=='help':
      if replace in columns: pprint(np.asarray([[i for idx,i in enumerate(old) if idx in [0,columns.index(replace)]],[i for idx,i in enumerate(new) if idx in [0,columns.index(replace)]]]), names=['id',replace])    
      elif replace=='help': _help()
      replace = raw_input("\nKeep both records [k]? Or replace [r], complete [c], or keep only [Press *Enter*] record {}? (Type column name to inspect or 'help' for options): ".format(duplicate[0]))

    if replace in options:
      
      # Replace the old record with the new record
      if replace=='r':
        sure = raw_input('Are you sure you want to replace record {} with record {}? [y/n] : '.format(*duplicate))
        if sure.lower()=='y':
          self.modify("DELETE FROM {} WHERE id={}".format(table, duplicate[1]), verbose=False)
          self.modify("UPDATE {} SET id={} WHERE id={}".format(table, duplicate[0], duplicate[1]), verbose=False)
      
      # Replace specific records
      elif replace.startswith('r') and all([i in list(columns)+options for i in replace.split()]):
        empty_cols, new_vals = zip(*[['{}=?'.format(e),n] for e,n in zip(columns,new) if e in replace])
        if empty_cols:
          if delete: self.modify("DELETE FROM {} WHERE id={}".format(table, duplicate[1]), verbose=False)
          self.modify("UPDATE {} SET {} WHERE id={}".format(table, ','.join(empty_cols), duplicate[0]), tuple(new_vals), verbose=False)
      
      # Complete the old record with any missing data provided in the new record, then delete the new record
      elif replace=='c':
        try:
          empty_cols, new_vals = zip(*[['{}=?'.format(e),n] for e,o,n in zip(columns,old,new) if repr(o).lower() in ['','none','null'] and repr(n).lower() not in ['','none','null']])
          self.modify("DELETE FROM {} WHERE id={}".format(table, duplicate[1]), verbose=False)
          self.modify("UPDATE {} SET {} WHERE id={}".format(table, ','.join(empty_cols), duplicate[0]), tuple(new_vals), verbose=False)
        except:
          pass
      
      # Keep both records
      elif replace=='k': return 'keep'
      
      # Execute raw SQL
      elif replace.startswith('sql ') and 'sql' in options: 
        try: self.modify(replace[4:], verbose=False)
        except: print 'Invalid SQLite statement: {}'.format(replace[4:])
    
    # Abort the current database clean up
    elif replace=='abort': return 'abort'
    
    # Delete the current record if *delete=True
    elif not replace:
      self.modify("DELETE FROM {} WHERE id={}".format(table, duplicate[1]), verbose=False)
          
  def inventory(self, source_id, plot=False, fetch=False, fmt='table'):
    """
    Prints a summary of all objects in the database. Input string or list of strings in **ID** or **unum** for specific objects.
    
    Parameters
    ----------
    source_id: int
      The id from the SOURCES table whose data across all tables is to be printed.
    plot: bool
      Plots all spectra for the object.
    fetch: bool
      Return the results.
      
    Returns
    -------
    data_tables: dict
      Returns a dictionary of astropy tables with the table name as the keys.
    
    """
    data_tables = {}

    for table in ['sources']+[t for t in zip(*self.list("SELECT * FROM sqlite_master WHERE type='table'"))[1] if t not in ['sources','sqlite_sequence']]:
        
      try:
          
        # Get the columns, pull out redundant ones, and query the table for this source's data
        columns, types = self.query("PRAGMA table_info({})".format(table), unpack=True)[1:3]
        
        if table=='sources' or 'source_id' in columns:
        
          # Only get simple data types and exclude redundant 'source_id' for nicer printing
          columns = columns[((types=='REAL')|(types=='INTEGER')|(types=='TEXT'))&(columns!='source_id')]

          # Query the table
          try: 
            id = 'id' if table.lower()=='sources' else 'source_id'
            data = self.query("SELECT {} FROM {} WHERE {}={}".format(','.join(columns),table,id,source_id), fmt='table')
            
            if not data and table.lower()=='sources':
              print 'No source with id {}. Try db.search() to search the database for a source_id.'.format(source_id)
            
          except: 
            data = None
        
          # If there's data for this table, save it
          if data: 
            if fetch: 
              data_tables[table] = self.query("SELECT {} FROM {} WHERE {}={}".format(','.join(columns),table,id,source_id), \
                                              fetch=True, fmt=fmt)
            else:
              data = data[list(columns)]
              pprint(data, title=table.upper())          
        
          # Plot all the spectra
          if plot and table.lower()=='spectra' and data:
            for i in self.query("SELECT id FROM spectra WHERE source_id={}".format(source_id), unpack=True)[0]: 
              self.plot_spectrum(i)
              
        else: pass
    
      except:
        print 'Could not retrieve data from {} table.'.format(table.upper())
    
    if fetch: return data_tables

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
      con, master, reassign = get_db(conflicted), self.list("PRAGMA database_list").fetchall()[0][2], {}
      con.modify("ATTACH DATABASE '{}' AS m".format(master), verbose=False)
      self.modify("ATTACH DATABASE '{}' AS c".format(conflicted), verbose=False)
      con.modify("ATTACH DATABASE '{}' AS c".format(conflicted), verbose=False)
      self.modify("ATTACH DATABASE '{}' AS m".format(master), verbose=False)
      
      # Drop any backup tables from failed merges
      for table in tables: self.modify("DROP TABLE IF EXISTS Backup_{0}".format(table), verbose=False)
      
      # Gather user data to add to CHANGELOG table
      import socket, datetime
      if not diff_only: user = raw_input('Please enter your name : ')
      machine_name = socket.gethostname()
      date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
      modified_tables = []
      
      # Merge table by table, starting with SOURCES
      tables = tables or ['sources']+[t for t in zip(*self.list("SELECT * FROM sqlite_master WHERE name NOT LIKE '%Backup%' AND name!='sqlite_sequence' AND type='table'{}".format(" AND name IN ({})".format("'"+"','".join(tables)+"'") if tables else '')))[1] if t!='sources']
      for table in tables:
        # Get column names and data types from master table and column names from conflicted table
        columns, types, constraints = self.query("PRAGMA table_info({})".format(table), unpack=True)[1:4]
        conflicted_cols = con.query("PRAGMA table_info({})".format(table), unpack=True)[1]
                
        if any([i not in columns for i in conflicted_cols]):
          # Abort table merge if conflicted has new columns not present in master. New columns must be added to the master database first via db.edit_columns().
          print "\nMerge of {0} table aborted since conflicted copy has columns {1} not present in master.\nAdd new columns to master with astrodb.table() method and try again.\n".format(table.upper(),[i for i in conflicted_cols if i not in columns])
        
        else:
          # Add new columns from master table to conflicted table if necessary
          if any([i not in conflicted_cols for i in columns]): 
            con.modify("DROP TABLE IF EXISTS Conflicted_{0}".format(table))
            con.modify("ALTER TABLE {0} RENAME TO Conflicted_{0}".format(table))
            con.modify("CREATE TABLE {0} ({1})".format(table, ', '.join(['{} {} {}{}'.format(c,t,r,' UNIQUE' if c=='id' else '') for c,t,r in zip(columns,types,constraints*['NOT NULL'])])))
            con.modify("INSERT INTO {0} ({1}) SELECT {1} FROM Conflicted_{0}".format(table, ','.join(conflicted_cols)))
            con.modify("DROP TABLE Conflicted_{0}".format(table))
        
          # Pull unique records from conflicted table
          data = map(list, con.list("SELECT * FROM (SELECT 1 AS db, {0} FROM m.{2} UNION ALL SELECT 2 AS db, {0} FROM c.{2}) GROUP BY {1} HAVING COUNT(*)=1 AND db=2".format(','.join(columns), ','.join(columns[1:]), table)).fetchall())
          
          if data:
            
            # Just print the table differences
            if diff_only:
              pprint(zip(*data)[1:], names=columns, title='New {} records'.format(table.upper()))
            
            # Add new records to the master and then clean up tables
            else:
              # Make temporary table copy so changes can be undone at any time
              self.modify("DROP TABLE IF EXISTS Backup_{0}".format(table), verbose=False)
              self.modify("ALTER TABLE {0} RENAME TO Backup_{0}".format(table), verbose=False)
              self.modify("CREATE TABLE {0} ({1})".format(table, ', '.join(['{} {} {}{}'.format(c,t,r,' UNIQUE' if c=='id' else '') for c,t,r in zip(columns,types,constraints*['NOT NULL'])])), verbose=False)
              self.modify("INSERT INTO {0} ({1}) SELECT {1} FROM Backup_{0}".format(table, ','.join(columns)), verbose=False)

              # Create a dictionary of any reassigned ids from merged SOURCES tables and replace applicable source_ids in other tables.
              print "\nMerging {} tables.\n".format(table.upper())
              try: count = self.query("SELECT MAX(id) FROM {}".format(table), fetch='one')[0]+1
              except TypeError: count = 1
              for n,i in enumerate([d[1:] for d in data]):
                if table=='sources': reassign[i[0]] = count
                elif 'source_id' in columns and i[1] in reassign.keys(): i[1] = reassign[i[1]]
                else: pass
                i[0] = count
                data[n] = i
                count += 1
            
              # Insert unique records into master
              for d in data: self.modify("INSERT INTO {} VALUES({})".format(table, ','.join(['?' for c in columns])), d, verbose=False)
              pprint(zip(*data), names=columns, title="{} records added to {} table at '{}':".format(len(data), table, master))
              
              # Run clean_up on the table to check for conflicts
              abort = self.clean_up(table)
          
              # Undo all changes to table if merge is aborted. Otherwise, push table changes to master.
              if abort: 
                self.modify("DROP TABLE {0}".format(table), verbose=False)
                self.modify("ALTER TABLE Backup_{0} RENAME TO {0}".format(table), verbose=False)
              else: 
                self.modify("DROP TABLE Backup_{0}".format(table), verbose=False)
                modified_tables.append(table.upper())
          
          else: print "\n{} tables identical.".format(table.upper())
      
      # Add data to CHANGELOG table
      if not diff_only:
        user_description = raw_input('\nPlease describe the changes made in this merge: ')
        self.list("INSERT INTO changelog VALUES(?, ?, ?, ?, ?, ?, ?)", \
                    (None, date, str(user), machine_name, ', '.join(modified_tables), user_description, os.path.basename(conflicted)))
      
      # Finish up and detach
      if diff_only:
        print "\nDiff complete. No changes made to either database. Set `diff_only=False' to apply merge."
      else:
        print "\nMerge complete!"
        
      con.modify("DETACH DATABASE c", verbose=False)
      self.modify("DETACH DATABASE c", verbose=False)
      con.modify("DETACH DATABASE m", verbose=False)
      self.modify("DETACH DATABASE m", verbose=False)
    else: print "File '{}' not found!".format(conflicted)

  def modify(self, SQL, params='', verbose=True):
    """
    Wrapper for CRUD operations to make them distinct from queries and automatically pass commit() method to cursor.
    
    Parameters
    ----------
    SQL: str
      The SQL query to execute
    params: sequence
      Mimicks the native parameter substitution of sqlite3
    verbose: bool
      Prints the number of modified records
    """
    try:
      
      # Make sure the database isn't locked
      self.conn.commit()
      
      if SQL.lower().startswith('select'):
        print 'Use self.query method for queries.'
      else:
        self.list(SQL, params)
        self.conn.commit()
        if verbose:
          print 'Number of records modified: {}'.format(self.list("SELECT changes()").fetchone()[0] or '0')
    except:
      print "Could not execute: "+SQL
    
  def output_spectrum(self, spectrum_id, filepath, original=False):
    """
    Prints a file of the spectrum with id **spectrum_id** to an ascii file with specified **filepath**.
    
    Parameters
    ----------
    spectrum_id: int
      The id from the SPECTRA table of the spectrum to print to file.
    filepath: str
      The path of the file to print the data to.
    original: bool
      Return the file in the original uploaded form
    
    """
    # Add option to just copy original file
    
    data = self.query("SELECT * FROM spectra WHERE id={}".format(spectrum_id), fetch='one', fmt='dict')
    if data:
      fn = '{}{}.txt'.format(filepath, data['filename'] or spectrum_id)

      # Write the header
      header = np.asarray(data['header'].cards)
      for h in header: h[0] = '# '+h[0]
      if data['header']: ii.write(header, fn, delimiter='\t', format='no_header')
      
      # Write the data
      with open(fn, mode='a') as f: ii.write([np.asarray(i, dtype=np.float64) for i in data['spectrum']], f, names=['# wavelength [{}]'.format(data['wavelength_units']),'flux [{}]'.format(data['flux_units']),'unc [{}]'.format(data['flux_units'])], delimiter='\t')
      
    else: print "No spectrum found with id {}".format(spectrum_id)
  
  def plot_spectrum(self, spectrum_id, overplot=False, color='b'):
    """
    Plots spectrum **ID** from SPECTRA table.
    
    Parameters
    ----------
    spectrum_id: int
      The id from the SPECTRA table of the spectrum to plot.
    overplot: bool
      Overplot the spectrum
    color: str
      The color used for the data
      
    """
    i = self.query("SELECT * FROM spectra WHERE id={}".format(spectrum_id), fetch='one', fmt='dict')
    if i:
      try:
        spec = i['spectrum']
        
        # Draw the axes and add the metadata
        if not overplot: 
          fig, ax = plt.subplots()
          plt.rc('text', usetex=False)
          ax.set_yscale('log', nonposy='clip'), plt.title('source_id = {}'.format(i['source_id']))
          plt.figtext(0.15,0.88, '{}\n{}\n{}\n{}'.format(i['filename'],self.query("SELECT name FROM telescopes WHERE id={}".format(i['telescope_id']), fetch='one')[0] if i['telescope_id'] else '',self.query("SELECT name FROM instruments WHERE id={}".format(i['instrument_id']), fetch='one')[0] if i['instrument_id'] else '',i['obs_date']), verticalalignment='top')
          ax.set_xlabel('[{}]'.format(i['wavelength_units'])), ax.set_ylabel('[{}]'.format(i['flux_units'])), ax.legend(loc=8, frameon=False)
        else: ax = plt.gca()
        
        # Plot the data
        ax.loglog(spec.data[0], spec.data[1], c=color, label='spec_id: {}'.format(i['id']))
        X, Y = plt.xlim(), plt.ylim()
        try: ax.fill_between(spec.data[0], spec.data[1]-spec.data[2], spec.data[1]+spec.data[2], color=color, alpha=0.3), ax.set_xlim(X), ax.set_ylim(Y)
        except: print 'No uncertainty array for spectrum {}'.format(spectrum_id)
      except: print "Could not plot spectrum {}".format(spectrum_id); plt.close()
    else: print "No spectrum {} in the SPECTRA table.".format(spectrum_id)

  def query(self, SQL, params='', fmt='array', fetch='all', unpack=False, export='', verbose=False):
    """
    Wrapper for cursors so data can be retrieved as a list or dictionary from same method
    
    Parameters
    ----------
    SQL: str
      The SQL query to execute
    params: sequence
      Mimicks the native parameter substitution of sqlite3
    fmt: str
      Returns the data as a dictionary, array, or astropy.table given 'dict', 'array', or 'table'
    unpack: bool
      Returns the transpose of the data
    export: str
      The file path of the ascii file to which the data should be exported
    verbose: bool
      Print the data also
      
    Returns
    -------
    result: (array,dict,table)
      The result of the database query
    """
    try:
      # Restricy queries to SELECT and PRAGMA statements
      if SQL.lower().startswith('select') or SQL.lower().startswith('pragma'):

        # Make the query explicit so that column and table names are preserved
        SQL, columns = self._explicit_query(SQL)
        
        # Get the data as a dictionary
        dictionary = self.dict(SQL, params).fetchall()
        
        if any(dictionary):
        
          # Fetch one
          if fetch=='one': dictionary = [dictionary.pop(0)]

          # Make an Astropy table
          table = at.Table(dictionary)
                
          # Reorder the columns
          table = table[columns]
          
          # Make an array
          array = np.asarray(table)

          # Unpack the results if necessary (data types are not preserved)
          if unpack: array = np.array(zip(*array))
        
          # Print on screen
          if verbose: pprint(table)         
        
          # Print the results to file
          if export:
            # If .vot or .xml, assume VOTable export with votools
            if export.lower().endswith('.xml') or export.lower().endswith('.vot'): votools.dict_tovot(dictionary, export)
          
            # Otherwise print as ascii
            else: ii.write(table, export, Writer=ii.FixedWidthTwoLine, fill_values=[('None', '-')])
        
          # Or return the results
          else: 
            if fetch=='one': dictionary, array = dictionary[0], array if unpack else np.array(list(array[0]))
            return table if fmt=='table' else dictionary if fmt=='dict' else array
          
        else: return
          
      else:
        print 'Queries must begin with a SELECT or PRAGMA statement. For database modifications use self.modify() method.'  
    
    except IOError:
      print 'Could not execute: '+SQL

  def search(self, criterion, table, columns='', fetch=False):
    """
    General search method for tables. For (ra,dec) input in decimal degrees, 
    i.e. '(12.3456,-65.4321)', returns all sources within 1 arcminute.
    For string input, i.e. 'vb10', returns all sources with case-insensitive partial text 
    matches in columns with 'TEXT' data type. For integer input, i.e. 123, returns all 
    exact matches of columns with INTEGER data type.
    
    Parameters
    ----------
    criterion: (str, int, tuple)
      The text, integer, or coordinate tuple to search the table with.
    table: str
      The name of the table to search
    columns: sequence
      Specific column names to search, otherwise searches all columns
    fetch: bool
      Return the results of the query as an Astropy table
      
    """
    results = ''
    
    # Get list of columns to search and format properly
    all_columns, types = self.query("PRAGMA table_info({})".format(table), unpack=True)[1:3]
    columns = columns or all_columns
    columns = np.asarray([columns] if isinstance(columns,str) else columns)
    
    # Separate good and bad columns and corresponding types
    badcols = columns[~np.in1d(columns,all_columns)]
    columns = columns[np.in1d(columns,all_columns)]
    columns = np.array([c for c in all_columns if c in columns])
    types = np.array([t for c,t in zip(all_columns,types) if c in columns])[np.in1d(columns,all_columns)]
    for col in badcols: print "'{}' is not a column in the {} table.".format(col,table.upper())
            
    # Coordinate search
    if isinstance(criterion,(tuple,list,np.ndarray)):
      try:
        q = "SELECT * FROM {} WHERE ra BETWEEN ".format(table)\
          +str(criterion[0]-0.01667)+" AND "\
          +str(criterion[0]+0.01667)+" AND dec BETWEEN "\
          +str(criterion[1]-0.01667)+" AND "\
          +str(criterion[1]+0.01667)
        results = self.query(q, fmt='table')
      except:
        print "Could not search {} table by coordinates {}. Try again.".format(table.upper(),criterion)
    
    # Text string search of columns with 'TEXT' data type
    elif isinstance(criterion, (str,unicode)) and any(columns) and 'TEXT' in types:
      try:
        q = "SELECT * FROM {} WHERE {}".format(table,' OR '.join([r"REPLACE("+c+r",' ','') like '%"\
          +criterion.replace(' ','')+r"%'" for c,t in zip(columns,types[np.in1d(columns,all_columns)]) if t=='TEXT']))
        results = self.query(q, fmt='table')
      except:
        print "Could not search {} table by string {}. Try again.".format(table.upper(),criterion)
    
    # Integer search of columns with 'INTEGER' data type
    elif isinstance(criterion, int):
      try:
        q = "SELECT * FROM {} WHERE {}".format(table,' OR '.join(['{}={}'.format(c,criterion) \
            for c,t in zip(columns,types[np.in1d(columns,all_columns)]) if t=='INTEGER']))
        results = self.query(q, fmt='table')
      except:
        print "Could not search {} table by id {}. Try again.".format(table.upper(),criterion)        
    
    # Problem!
    else: print "Could not search {} table by '{}'. Try again.".format(table.upper(),criterion)
    
    # Print or return the results
    if results: 
      if fetch: return results
      else: pprint(results, title=table.upper())
    else: print "No results found for {} in {} the table.".format(criterion,table.upper())

  def table(self, table, columns, types, constraints='', new_table=False):
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
    new_table: bool
      Create a new table
    
    """
    goodtogo = True
    
    # Make sure there is an integer primary key, unique, not null 'id' column 
    # and the appropriate number of elements in each sequence
    if columns[0]!='id':
      print "Column 1 must be called 'id'"; goodtogo = False
    if types[0].upper()!='INTEGER PRIMARY KEY': 
      print "'id' column type must be 'INTEGER PRIMARY KEY'"; goodtogo = False
    if constraints:
      if 'UNIQUE' not in constraints[0].upper() and 'NOT NULL' not in constraints[0].upper(): 
        print "'id' column constraints must be 'UNIQUE NOT NULL'";  goodtogo = False
    else:
      constraints = ['UNIQUE NOT NULL']+(['']*(len(columns)-1))
    if not len(columns)==len(types)==len(constraints):
      print "Must provide equal length *columns ({}), *types ({}), and *constraints ({}) sequences."\
            .format(len(columns),len(types),len(constraints));  goodtogo = False
    
    if goodtogo:
      tables = self.query("SELECT name FROM sqlite_master", unpack=True)[0]
      
      # If the table exists, modify the columns
      if table in tables and not new_table:
        
        # Rename the old table and create a new one
        self.list("ALTER TABLE {0} RENAME TO TempOldTable".format(table))
        self.list("CREATE TABLE {0} ({1})".format(table, ', '.join(['{} {} {}'.format(c,t,r) for c,t,r in zip(columns,types,constraints)])))
        
        # Populate the new table and drop the old one
        old_columns = [c for c in self.query("PRAGMA table_info(TempOldTable)", unpack=True)[1] if c in columns]
        self.list("INSERT INTO {0} ({1}) SELECT {1} FROM TempOldTable".format(table, ','.join(old_columns)))
        self.list("DROP TABLE TempOldTable")
    
      # If the table does not exist and new_table is True, create it
      elif table not in tables and new_table:
        self.list("CREATE TABLE {0} ({1})".format(table, ', '.join(['{} {} {}'.format(c,t,r) for c,t,r in zip(columns,types,constraints)])))
    
      # Otherwise the table to be modified doesn't exist or the new table to add already exists, so do nothing
      else:
        if new_table: print 'Table {} already exists. Set *new_table=False to modify.'.format(table.upper())
        else: print 'Table {} does not exist. Could not modify. Set *new_table=True to add a new table.'.format(table.upper())
        
    else: print 'The {} table has not been {}. Please make sure your table columns, \
                types, and constraints are formatted properly.'.format(table.upper(),\
                'created' if new_table else 'modified')

  def _explicit_query(self, SQL):
    """
    Sorts the column names so they are returned in the same order they are queried. Also turns 
    ambiguous SELECT statements into explicit SQLite language in case column names are not unique.
    
    Parameters
    ----------
    SQL: str
      The SQLite query to parse
    
    Returns
    -------
    (SQL, columns): (str, sequence)
      The new SQLite string to use in the query and the ordered column names
    
    """
    # If field names are given, sort so that they come out in the same order they are fetched
    if 'select' in SQL.lower() and 'from' in SQL.lower():
                      
      # Make a dictionary of the table aliases
      tdict = {}
      from_clause = SQL.lower().split('from ')[-1].split(' where')[0]
      tables = [j for k in [i.split(' on ') for i in from_clause.split(' join ')] for j in k if '=' not in j]
      
      for t in tables:
        t = t.replace('as','')
        try:
          name, alias = t.split()
          tdict[alias] = name
        except:
          tdict[t] = t
      
      # Get all the column names
      columns = SQL.replace(' ','').lower().split('distinct' if 'distinct' in SQL.lower() else 'select')[1].split('from')[0].split(',')
  
      # Replace * with the field names
      for n,col in enumerate(columns):
        if '.' in col:
          t, col = col.split('.')
        else:
          t = tables[0]
  
        if '*' in col: 
          col = np.array(self.list("PRAGMA table_info({})".format(tdict.get(t))).fetchall()).T[1]
        else: 
          col = [col]
          
        columns[n] = ["{}.{}".format(t,c) if len(tables)>1 else c for c in col]
  
      # Flatten the list of columns
      columns = [j for k in columns for j in k]
  
      # Reconstruct SQL query
      SQL = "SELECT {}".format('DISTINCT ' if 'distinct' in SQL.lower() else '')\
          +','.join(["{0} AS '{0}'".format(col) for col in columns])\
          +' FROM '\
          +SQL.replace('from','FROM').split('FROM')[-1]
  
    elif 'pragma' in SQL.lower():
      columns = ['cid','name','type','notnull','dflt_value','pk']
    
    return SQL, columns

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
  return buffer(out.read())

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
  if isinstance(spec,str): pass
  else: spec = '$BDNYC_spectra'+spec.path.split('BDNYC_spectra')[1]
  
  return spec
  
  
def convert_spectrum(File):
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
  
  class Spectrum:
    def __init__(self, data, header, path):
        self.data = data
        self.header = header
        self.path = path
  
  if isinstance(File,str):
    
    # Convert variable path to absolute path
    if File.startswith('$'):
      abspath = os.popen('echo {}'.format(File.split('/')[0])).read()[:-1]
      if abspath: File = File.replace(File.split('/')[0],abspath)
    
    # For FITS files
    if File.endswith('.fits'):
      try: 
        # Get the data 
        spectrum, header = pf.getdata(File, cache=True, header=True)
      
        # Check the key type
        KEY_TYPE = ['CTYPE1']
        setType  = set(KEY_TYPE).intersection(set(header.keys()))
        if len(setType) == 0: isLinear = True
        else:
            valType = header[setType.pop()]
            isLinear = valType.strip().upper()=='LINEAR'
        
        # Get wl, flux & error data from fits file
        spectrum = __get_spec(spectrum, header, File)
    
        # Generate wl axis when needed
        if not spectrum[0]: spectrum[0] = __create_waxis(header, len(spectrum[1]), File)
        
        # If no wl axis generated, then clear out all retrieved data for object
        if not spectrum[0]: spectrum = None
      
      except: pass
    
    # For .txt files
    if File.endswith('.txt'): 
      try: 
        spectrum = ii.read(File)
        spectrum = np.array([np.asarray(spectrum.columns[n]) for n in [0,1,2]])
        try: 
          txt, header = open(File), []
          for i in txt: 
            if any([i.startswith(char) for char in ['#','|','\\']]): header.append(i.replace('\n',''))
          txt.close()
        except: pass
      except: pass

  if spectrum=='': print 'Could not retrieve spectrum at {}.'.format(File); return File
  else: 
    spectrum = Spectrum(spectrum, header, File)
    return spectrum

def __create_waxis(fitsHeader, lenData, fileName, wlog=False):
    # Define key names in
    KEY_MIN  = ['COEFF0','CRVAL1']         # Min wl
    KEY_DELT = ['COEFF1','CDELT1','CD1_1'] # Delta of wl
    KEY_OFF  = ['LTV1']                    # Offset in wl to subsection start
    
    # Find key names for minimum wl, delta, and wl offset in fits header
    setMin  = set(KEY_MIN).intersection(set(fitsHeader.keys()))
    setDelt = set(KEY_DELT).intersection(set(fitsHeader.keys()))
    setOff  = set(KEY_OFF).intersection(set(fitsHeader.keys()))
    
    # Get the values for minimum wl, delta, and wl offset, and generate axis
    if len(setMin) >= 1 and len (setDelt) >= 1:
        nameMin = setMin.pop()
        valMin  = fitsHeader[nameMin]
        
        nameDelt = setDelt.pop()
        valDelt  = fitsHeader[nameDelt]
        
        if len(setOff) == 0:
            valOff = 0
        else:
            nameOff = setOff.pop()
            valOff  = fitsHeader[nameOff]
        
        # generate wl axis
        if nameMin == 'COEFF0' or wlog==True:
            # SDSS fits files
            wAxis = 10 ** (np.arange(lenData) * valDelt + valMin)
        else:
            wAxis = (np.arange(lenData) * valDelt) + valMin - (valOff * valDelt)
        
    else:
        wAxis = None
        if verb:
            print 'Could not re-create wavelength axis for ' + fileName + '.'
    
    return wAxis

def __get_spec(fitsData, fitsHeader, fileName):
    validData = [None] * 3
    
    # Identify number of data sets in fits file
    dimNum = len(fitsData)
    
    # Identify data sets in fits file
    fluxIdx  = None
    waveIdx  = None
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
            fluxIdx  = 0
        else:
            waveIdx = 0
            fluxIdx = 1
    elif dimNum == 3:
        waveIdx  = 0
        fluxIdx  = 1
        sigmaIdx = 2
    elif dimNum == 4:
    # 0-flux clean, 1-flux raw, 2-background, 3-sigma clean
        fluxIdx  = 0
        sigmaIdx = 3
    elif dimNum == 5:
    # 0-flux, 1-continuum substracted flux, 2-sigma, 3-mask array, 4-unknown
        fluxIdx  = 0
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
            fitsData = [tmpWave,tmpFlux]
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
            print 'Unable to interpret data in ' + fileName + '.'
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
sqlite3.register_converter("ARRAY", convert_array)
sqlite3.register_converter("SPECTRUM", convert_spectrum)

def pprint(data, names='', title=''):
  """
  Prints tables with a little bit 'o formatting
  
  Parameters
  ----------
  data: (sequence, dict, table)
    The data to print in the table
  names: sequence
    The column names
  title: str (optional)
    The title of the table
    
  """
  # Make the data into a table if it isn't already
  if type(data)!=at.Table: data = at.Table(data, names=names)
  
  # Put the title in the metadata
  try: title = title or data.meta['name']
  except: pass
  
  # Shorten the column names for slimmer data
  for old,new in zip(*[data.colnames,[i.replace('wavelength','wav').replace('publication','pub').replace('instrument','inst').replace('telescope','scope') for i in data.colnames]]): 
    data.rename_column(old,new) if new!=old else None
  
  # Print it!
  if title: print '\n'+title
  ii.write(data, sys.stdout, Writer=ii.FixedWidthTwoLine, formats={'comments': '%.15s', 'obs_date': '%.10s', 'names': '%.20s', 'description': '%.50s'}, fill_values=[('None', '-')])

def clean_header(header):
  try:
    header = pf.open(File, ignore_missing_end=True)[0].header
    new_header = pf.Header()
    for x,y,z in header.cards: new_header[x.replace('.','_').replace('#','')] = (y,z)
    header = pf.PrimaryHDU(header=new_header).header
  except: pass

def _help():
  print ' '
  pprint(np.asarray([['[column name]','Display full record entry for that column without taking action'], \
                   ['k','Keeps both records and assigns second one new id if necessary'], \
                   ['r','Replaces all columns of first record with second record values'], \
                   ['r [column name] [column name]...','Replaces specified columns of first record with second record values'], \
                   ['c','Complete empty columns of first record with second record values where possible'], \
                   ['[Enter]','Keep first record and delete second'],\
                   ['abort','Abort merge of current table, undo all changes, and proceed to next table']]), \
                   names=['Command','Result'])

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
  if record['spectrum'].endswith('.fits'):
    header = pf.getheader(record['spectrum'])

    # Wavelength units
    if not record['wavelength_units']:
      try: 
        record['wavelength_units'] = header['XUNITS'] 
      except KeyError:
        try:
           if header['BUNIT']: record['wavelength_units'] = 'um'
        except KeyError: pass
    if 'microns' in record['wavelength_units'] or 'Microns' in record['wavelength_units'] or 'um' in record['wavelength_units']: record['wavelength_units'] = 'um'

    # Flux units
    if not record['flux_units']:
      try: record['flux_units'] = header['YUNITS'].replace(' ','')
      except KeyError:
        try: record['flux_units'] = header['BUNIT'].replace(' ','')
        except KeyError: pass
    if 'erg' in record['flux_units'] and 'A' in record['flux_units']: record['flux_units'] = 'ergs-1cm-2A-1' if 'erg' in record['flux_units'] and 'A' in record['flux_units'] else 'ergs-1cm-2um-1' if 'erg' in record['flux_units'] and 'um' in record['flux_units'] else 'Wm-2um-1' if 'W' in record['flux_units'] and 'um' in record['flux_units'] else 'Wm-2A-1' if 'W' in record['flux_units'] and 'A' in record['flux_units'] else ''

    # Observation date
    if not record['obs_date']:
      try: record['obs_date'] = header['DATE_OBS']
      except KeyError:
        try: record['obs_date'] = header['DATE-OBS']
        except KeyError:
          try: record['obs_date'] = header['DATE']
          except KeyError: pass

    # Telescope id
    if not record['telescope_id']:
      try:
        n = header['TELESCOP'].lower() if isinstance(header['TELESCOP'],str) else ''
        record['telescope_id'] = 5 if 'hst' in n else 6 if 'spitzer' in n else 7 if 'irtf' in n else 9 if 'keck' in n and 'ii' in n else 8 if 'keck' in n and 'i' in n else 10 if 'kp' in n and '4' in n else 11 if 'kp' in n and '2' in n else 12 if 'bok' in n else 13 if 'mmt' in n else 14 if 'ctio' in n and '1' in n else 15 if 'ctio' in n and '4' in n else 16 if 'gemini' in n and 'north' in n else 17 if 'gemini' in n and 'south' in n else 18 if ('vlt' in n and 'U2' in n) else 19 if '3.5m' in n else 20 if 'subaru' in n else 21 if ('mag' in n and 'ii' in n) or ('clay' in n) else 22 if ('mag' in n and 'i' in n) or ('baade' in n) else 23 if ('eso' in n and '1m' in n) else 24 if 'cfht' in n else 25 if 'ntt' in n else 26 if ('palomar' in n and '200-inch' in n) else 27 if 'pan-starrs' in n else 28 if ('palomar' in n and '60-inch' in n) else 29 if ('ctio' in n and '0.9m' in n) else 30 if 'soar' in n else 31 if ('vlt' in n and 'U3' in n) else 32 if ('vlt' in n and 'U4' in n) else 33 if 'gtc' in n else None
      except KeyError: pass

    # Instrument id
    if not record['instrument_id']:
      try: 
        i = header['INSTRUME'].lower()
        record['instrument_id'] = 1 if 'r-c spec' in i or 'test' in i or 'nod' in i else 2 if 'gmos-n' in i else 3 if 'gmos-s' in i else 4 if 'fors' in i else 5 if 'lris' in i else 6 if 'spex' in i else 7 if 'ldss3' in i else 8 if 'focas' in i else 9 if 'nirspec' in i else 10 if 'irs' in i else 11 if 'fire' in i else 12 if 'mage' in i else 13 if 'goldcam' in i else 14 if 'sinfoni' in i else 15 if 'osiris' in i else 16 if 'triplespec' in i else 17 if 'x-shooter' in i else 18 if 'gnirs' in i else 19 if 'wircam' in i else 20 if 'cormass' in i else 21 if 'isaac' in i else 22 if 'irac' in i else 23 if 'dis' in i else 24 if 'susi2' in i else 25 if 'ircs' in i else 26 if 'nirc' in i else 29 if 'stis' in i else 0
      except KeyError: pass

  return record
  
type_dict = {'INTEGER':np.dtype('int64'), 'REAL':np.dtype('float64'), 'TEXT':np.dtype('S64'), 'ARRAY':np.dtype('object'), 'SPECTRUM':np.dtype('S128')}
