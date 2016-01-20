#!/usr/bin/python
# Author: Joe Filippazzo, jcfilippazzo@gmail.com

import io, os, sys, itertools, urllib2, sqlite3 as sql, numpy as np, matplotlib.pyplot as plt, astropy.io.fits as pf, astropy.io.ascii as ii, astropy.table as at
package_path = os.path.dirname(os.path.abspath(__file__))

def clean_header(fitsPath):
  """
  Fix broken fits headers by cleaning illegal characters from keywords, inserting an END card, and rewriting the header.
  
  Parameters
  ----------
  fitsPath: str
    The path of the FITS file to be repaired.
  
  Returns
  -------
  fits.header
    The FITS header with an END card inserted and illegal characters removed from keywords.
    
  """
  header = pf.open(fitsPath, ignore_missing_end=True)[0].header
  new_header = pf.Header()
  for x,y,z in header.cards: new_header[x.replace('.','_').replace('#','')] = (y,z)
  return pf.PrimaryHDU(header=new_header).header

def create_database(dbpath):
  """
  Create a new database at the given *dbpath.
  
  Parameters
  ----------
  dbpath: str
    The full path for the new database, including the filename and .db file extension.

  """
  if dbpath.endswith('.db'):
    sources_table = "CREATE TABLE sources (id INTEGER PRIMARY KEY, ra REAL, dec REAL, designation TEXT, publication_id INTEGER, comments TEXT, shortname TEXT, names TEXT)"
    os.system("sqlite3 {} '{}'".format(dbpath,sources_table))
    if os.path.isfile(dbpath):
      print "\nDatabase created! To load, run\n\ndb = astrodb.get_db('{}')\n\nThen run db.modify_table() method to create tables.".format(dbpath)
  else: print "Please provide a path and file name with a .db file extension, e.g. /Users/Me/Desktop/test.db"

class get_db:
  def __init__(self, dbpath=package_path+'/bdnyc198.db'):
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
      def dict_factory(cursor, row):
        d = {}
        for idx,col in enumerate(cursor.description): d[col[0]] = row[idx]
        return d
    
      con = sql.connect(dbpath, isolation_level=None, detect_types=sql.PARSE_DECLTYPES)
      con.text_factory = str
      self.conn = con
      self.list = con.cursor().execute
      self.dict = con.cursor()
      self.dict.row_factory = dict_factory
      self.dict = self.dict.execute
    
    else: print "Sorry, no such file '{}'".format(dbpath)

  def add_data(self, ascii, table, delimiter='|', multiband=False):
    """
    Adds data in **ascii** file to the specified database **table**. Note column names (row 1 of ascii file) must match table fields to insert, however order and completeness don't matter.
    
    Parameters
    ----------
    ascii: str
      The path to the ascii file to be read in.
    table: str
      The name of the table into which the data should be inserted
    delimiter: str
      The string to use as the delimiter when parsing the ascii file
    multiband: bool
      Digest columns of multiple photometric measurements (e.g. J, H, Ks) into individual rows of data for database insertion

    Returns
    -------
    None
    
    """
    # Digest the ascii file into Python 
    data, insert, update = np.genfromtxt(ascii, delimiter=delimiter, dtype=object).tolist(), [], []
    
    # Get the column names and data types from the table
    columns, types = zip(*self.list("PRAGMA table_info({})".format(table)).fetchall())[1:3]

    # Grab the column names from the first line of the input ascii file 
    data_columns = data.pop(0)
    
    # If a row contains photometry for multiple bands, use the *multiband argument and execute this
    if multiband and table=='photometry':    
      Tdata, new_data = zip(*data), []
      
      # Recognized filters for multiband photometry upload. Update this to pull from managed list such as SVO Filter Profile Service
      bands = ['HST_F336W', 'GALEX_FUV', 'DES_g', 'HST_F775W', 'HST_F190N', 'DES_r', 'MIPS_[24]', 'HST_F090M', 'GAIA_BP', 'HST_F215N', 'GAIA_RP', 'IRAC_[8]', 'HST_F140W', 'GALEX_NUV', 'DENIS_J', 'DENIS_I', 'HST_F673N', 'HST_F164N', 'SDSS_r', 'HST_F170M', 'DES_Y', 'WISE_W3', 'HST_F475W', 'MKO_J', 'JOHNSON_B', 'DENIS_Ks', 'GAIA_G', 'COUSINS_I', '2MASS_H', '2MASS_J', 'HST_F555W', '2MASS_Ks', 'IRAC_[3.6]', 'COUSINS_R', 'JOHNSON_U', 'JOHNSON_V', 'HST_F110W', 'HST_F625W', 'MKO_K', 'IRAC_[5.8]', 'MKO_Y', 'SDSS_g', 'DES_u', 'SDSS_i', 'HST_F656N', 'DES_z', 'HST_F850LP', 'FS_J1', 'FS_J2', 'FS_J3', 'HST_F390N', 'MKO_M', 'WISE_W4', 'SDSS_u', 'WISE_W2', 'WISE_W1', 'MKO_H', 'DES_i', 'IRAC_[4.5]', 'SDSS_z', 'MKO_L']
      
      # Get the columns that are not magnitudes
      repeat_data = [Tdata[idx] for idx,c in enumerate(data_columns) if c.replace('_unc','') not in bands]
      repeat_cols = [c for c in data_columns if c.replace('_unc','') not in bands]+['band','magnitude','magnitude_unc']
      
      # For each band, make a new data row
      for c in data_columns:
        if c in bands:
          new_data += zip(*repeat_data+[[c]*len(data),Tdata[data_columns.index(c)],Tdata[data_columns.index(c+'_unc')]])
    
      # Transpose the data columns into rows and identify the new data columns to insert
      data, data_columns = new_data, repeat_cols
                
    # Insert data strictly matching SQL table field names to ascii input column names
    for row in data:
      values = [None for i in columns]
      for col in columns: values[columns.index(col)] = row[data_columns.index(col)] if col in data_columns and row[data_columns.index(col)] else None
      # values[0] = sorted(list(set(range(1,self.list("SELECT max(id) FROM {}".format(table)).fetchone()[0]+2))-set(zip(*self.list("SELECT id FROM {}".format(table)).fetchall())[0])))[0]
      if isinstance(values[1],int) or (values[1] or 'None').isdigit(): update.append(tuple(values)) if values[columns.index('id')] else insert.append(tuple(values))

    # Reject rows that fail column constraints, e.g. NOT NULL fields like 'source_id'
    constraints, keep = np.asarray(zip(*self.query("pragma table_info({})".format(table)))[3],dtype=bool), []
    for n,d in enumerate(insert):
      if all(np.asarray(d,dtype=bool)[constraints][1:]): keep.append(n)
    try: insert = list(np.array(insert)[np.array(keep)])
    except IndexError: insert = ''

    # If they are unique records (i.e. don't have an 'id' column value specified in the ascii file), add them as new records
    if insert:
      print at.Table(np.asarray(insert), names=columns)
      for i in insert: self.modify("INSERT INTO {} VALUES({})".format(table, ','.join('?'*len(columns))), i, verbose=False)
      print "{} new records added to the {} table.".format(len(insert),table.upper())

    # If they do have a specified 'id', update the fields for that record with the supplied information
    if update:
      print at.Table(np.asarray([['-'*30,'-'*100],['[column name]','Display full record entry for that column without taking action'],['k','Keeps both records and assigns second one new id if necessary'],['r','Replaces all columns of first record with second record values'],['r [column name] [column name]...','Replaces specified columns of first record with second record values'],['c','Complete empty columns of first record with second record values where possible'],['[Enter]','Keep first record and delete second'],['abort','Abort merge of current table, undo all changes, and proceed to next table']]), names=['Command','Result'])
      for item in update:
        record = self.list("SELECT * FROM {} WHERE id={}".format(table, item[columns.index('id')])).fetchone()
        if record: self.compare_records(self, table, columns, record, item)
    
    # Run table clean up
    self.clean_up(table)

  def add_spectrum(self, filepath, source_id, wavelength_units='', flux_units='', publication_id='', obs_date='', wavelength_order='', regime='', instrument_id='', telescope_id='', mode_id='', comment='', header_chars=['#']):
    """
    Adds an ascii or .fits file to the SQL database
    
    Parameters
    ----------
    filepath: str
      The path to the FITS or ascii file
    source_id: int
      The id from the SOURCES table of the object to which the spectrum will be associated.
    wavelength_units: str (optional)
      The wavelength units of the spectrum, e.g. 'um', 'A', 'nm'
    flux_units: str (optional)
      The flux units of the spectrum, e.g. 'erg/s/cm2/A', 'W m-2 um-1'
    publication_id: int (optional)
      The id from the PUBLICATIONS table that indicates the reference for the spectrum
    obs_date: str (optional)
      The date of the observation
    wavelength_order: int (optional)
      For high resolution spectra, the order of the observation, e.g. '65' for NIRSPEC order 65. Leave blank for low and medium resolution spectra.
    instrument_id: int (optional)
      The id from the INSTRUMENTS table used to take the spectrum
    telescope_id: int (optional)
      The id from the TELESCOPES table used to take the spectrum
    mode_id: int (optional)
      The id from the MODES table used to take the spectrum, if applicable
    regime: str (optional)
      The regime of the spectrum, e.g. 'OPT','NIR' or 'MIR'
    comment: str (optional)
      Any comments about the spectrum that might be useful to future users
    header_chars: sequence 
      If a line begins with any character in this list, that text will be used to generate a FITS header for the spectrum.
    """
    filename = os.path.basename(filepath)
    
    # If it's a FITS file, try to get the missing metadata from the header
    if filepath.endswith('.fits'):
      header = clean_header(filepath)

      # Wavelength units
      if not wavelength_units:
        try:
          wavelength_units = header['XUNITS'] 
          if 'microns' in wavelength_units or 'Microns' in wavelength_units or 'um' in wavelength_units: wavelength_units = 'um'
        except KeyError:
          try:
             if header['BUNIT']: wavelength_units = 'um'
          except KeyError: wavelength_units = ''
      
      # Flux units
      if not flux_units:
        try: flux_units = header['YUNITS'].replace(' ','')
        except KeyError:
          try: flux_units = header['BUNIT'].replace(' ','')
          except KeyError: flux_units = ''
      if 'erg' in flux_units and 'A' in flux_units: flux_units = 'ergs-1cm-2A-1' if 'erg' in flux_units and 'A' in flux_units else 'ergs-1cm-2um-1' if 'erg' in flux_units and 'um' in flux_units else 'Wm-2um-1' if 'W' in flux_units and 'um' in flux_units else 'Wm-2A-1' if 'W' in flux_units and 'A' in flux_units else ''

      # Observation date
      if not obs_date:
        try: obs_date = header['DATE_OBS']
        except KeyError:
          try: obs_date = header['DATE-OBS']
          except KeyError:
            try: obs_date = header['DATE']
            except KeyError: obs_date = ''
      
      # Telescope id
      if not telescope_id:
        try:
          n = header['TELESCOP'].lower() if isinstance(header['TELESCOP'],str) else ''
          telescope_id = 5 if 'hst' in n else 6 if 'spitzer' in n else 7 if 'irtf' in n else 9 if 'keck' in n and 'ii' in n else 8 if 'keck' in n and 'i' in n else 10 if 'kp' in n and '4' in n else 11 if 'kp' in n and '2' in n else 12 if 'bok' in n else 13 if 'mmt' in n else 14 if 'ctio' in n and '1' in n else 15 if 'ctio' in n and '4' in n else 16 if 'gemini' in n and 'north' in n else 17 if 'gemini' in n and 'south' in n else 18 if ('vlt' in n and 'U2' in n) else 19 if '3.5m' in n else 20 if 'subaru' in n else 21 if ('mag' in n and 'ii' in n) or ('clay' in n) else 22 if ('mag' in n and 'i' in n) or ('baade' in n) else 23 if ('eso' in n and '1m' in n) else 24 if 'cfht' in n else 25 if 'ntt' in n else 26 if ('palomar' in n and '200-inch' in n) else 27 if 'pan-starrs' in n else 28 if ('palomar' in n and '60-inch' in n) else 29 if ('ctio' in n and '0.9m' in n) else 30 if 'soar' in n else 31 if ('vlt' in n and 'U3' in n) else 32 if ('vlt' in n and 'U4' in n) else 33 if 'gtc' in n else None
        except KeyError: telescope_id = ''
      
      # Instrument id
      if not instrument_id:
        try: 
          i = header['INSTRUME'].lower()
          instrument_id = 1 if 'r-c spec' in i or 'test' in i or 'nod' in i else 2 if 'gmos-n' in i else 3 if 'gmos-s' in i else 4 if 'fors' in i else 5 if 'lris' in i else 6 if 'spex' in i else 7 if 'ldss3' in i else 8 if 'focas' in i else 9 if 'nirspec' in i else 10 if 'irs' in i else 11 if 'fire' in i else 12 if 'mage' in i else 13 if 'goldcam' in i else 14 if 'sinfoni' in i else 15 if 'osiris' in i else 16 if 'triplespec' in i else 17 if 'x-shooter' in i else 18 if 'gnirs' in i else 19 if 'wircam' in i else 20 if 'cormass' in i else 21 if 'isaac' in i else 22 if 'irac' in i else 23 if 'dis' in i else 24 if 'susi2' in i else 25 if 'ircs' in i else 26 if 'nirc' in i else 29 if 'stis' in i else 0
        except KeyError: instrument_id = ''

    # Pull comments out of text file (lines which begin with one of the specified *header_chars*) and create FITS header for database insertion
    if filepath.endswith('.txt'):
      h = [i.strip() for i in open(filepath) if any([i.startswith(char) for char in header_chars])]
      if h:
        header = pf.Header()
        for i in h: header['COMMENT'] = i
        header = pf.PrimaryHDU(header=header).header
      else: header = None

    try:
      # Upload the file to AWS S3
      import boto3
      s3 = boto3.resource('s3')
      data = open(filepath, 'rb')
      s3.Bucket('bdnyc').put_object(Key=filename, Body=data)

      # Return the URL for database insertion
      spectrum = 'https://s3.amazonaws.com/bdnyc/{}'.format(filename)
        
      # Add data to the database
      spec_id = sorted(list(set(range(1,self.list("SELECT max(id) FROM spectra").fetchone()[0]+2))-set(zip(*self.list("SELECT id FROM spectra").fetchall())[0])))[0]
      self.modify("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (spec_id, source_id, spectrum, wavelength_units, flux_units, wavelength_order, regime, publication_id, obs_date, instrument_id, telescope_id, mode_id, filename, comment, header))
      print at.Table(np.asarray([[spec_id, source_id, wavelength_units, flux_units, regime, publication_id, obs_date, instrument_id, telescope_id, mode_id, filename, comment]]), names=['spec_id','source_id','w_unit','f_unit','regime','pub_id','obs_date', 'inst_id', 'scope_id', 'mode_id', 'filename', 'comment'])
    except: print "Couldn't add spectrum {} to database.".format(filepath)
 
  def clean_up(self, table):
    """
    Removes exact duplicates, blank records or data without a *source_id* from the specified **table**. Then finds possible duplicates and prompts for conflict resolution.
    
    Parameters
    ----------
    table: str
      The name of the table to remove duplicates, blanks, and data without source attributions.
    
    """
    (columns, types), dup, ignore, I = self.query("PRAGMA table_info({})".format(table), unpack=True)[1:3], 1, [], ''
    
    # Delete blank records, exact duplicates, or data without a source_id
    self.modify("DELETE FROM {0} WHERE ({1})".format(table, columns[1]+' IS NULL' if len(columns)==2 else (' IS NULL AND '.join(columns[1:])+' IS NULL')))
    self.modify("DELETE FROM {0} WHERE id NOT IN (SELECT min(id) FROM {0} GROUP BY {1})".format(table,', '.join(columns[1:])))
    if 'source_id' in columns: self.list("DELETE FROM {0} WHERE source_id IS NULL OR source_id IN ('null','None','')".format(table))
    
    if table in ['sources','spectra','photometry','spectral_types','radial_velocities','parallaxes','proper_motions']:
      # Define columns to check for uniqueness
      primary, secondary, ignore = ['flux','magnitude','parallax','spectral_type','proper_motion_ra','proper_motion_dec','radial_velocity'], ['band','regime'], []
      
      # Find non-unique records and run astrodb.compare_records()
      while dup:
        # Check SOURCES table for records with similar RA and Dec but different ids that are not companions or system components.
        if table=='sources': dup = self.query(  "SELECT t1.*, t2.* FROM sources t1 JOIN sources t2 WHERE t1.id!=t2.id \
                                                AND (t1.ra BETWEEN t2.ra-0.0007 AND t2.ra+0.0007) AND (t1.dec BETWEEN t2.dec-0.00077 AND t2.dec+0.0007)\
                                                {}".format(' AND '+' AND '.join(['(t1.id NOT IN ({0}) AND t2.id NOT IN ({0}))'.format(','.join(map(str,i))) for i in ignore]) if ignore else ''), fetch='one')

        # Check all other tables for records with identical primary and secondary column values but different ids.        
        else: dup = self.query("SELECT t1.*, t2.* FROM {0} AS t1 JOIN {0} AS t2 ON t1.source_id=t2.source_id WHERE t1.id!=t2.id AND t1.{1}=t2.{1}{2}{3}".format(table, [c for c in columns if c in primary].pop(), ' AND t1.{0}=t2.{0}'.format([c for c in columns if c in secondary].pop()) if [c for c in columns if c in secondary] else '', ' AND '+' AND '.join(['(t1.id NOT IN ({0}) AND t2.id NOT IN ({0}))'.format(','.join(map(str,i))) for i in ignore]) if ignore else ''), fetch='one')        
        
        # Compare potential duplicates and prompt user for action on each
        if dup and dup[:len(dup)/2][0]!=dup[len(dup)/2:][0]:
          I = self.compare_records(self, table, columns, dup[:len(dup)/2], dup[len(dup)/2:], delete=True)
          if isinstance(I,list): ignore.append(I)
          elif I=='undo': pass # Add this functionality!
          elif I=='abort': break
          else: pass
    
    # Finish or abort table merge
    if I=='abort': 
      print '\nAborted merge of {} table. Undoing all changes.\n'.format(table.upper())
      return 'abort'
    else: print 'Finished clean up on {} table.'.format(table.upper())

  def compare_records(self, table, columns, old, new, options=['r','c','k','sql'], delete=False):
    """
    Compares similar records and prompts the user to make decisions about keeping, updating, or modifying records in question.
  
    Parameters
    ----------
    table: str
      The name of the table whose records are being compared.
    columns: list
      The list of columns across which the comparison should be made.
    old: (str, int, float, blob)
      The value of the record with the lower id.
    new: (str, int, float, blob)
      The value of the record with the higher id.
    options: list
      The allowed options: 'r' for replace, 'c' for complete, 'k' for keep, 'sql' for raw SQL input.
    delete: bool
      Delete the record with the higher id.
    
    """
    # Print the old and new records suspectred of being duplicates
    pold, pnew = ['{:.3g}...{:.3g}'.format(i[0],i[-1]) if isinstance(i, np.ndarray) else i for i in old], ['{:.3g}...{:.3g}'.format(i[0],i[-1]) if isinstance(i, np.ndarray) else i for i in new]
    print at.Table(np.asarray([pold,pnew]), names=columns)

    # Print the command key
    replace = raw_input("Keep both records [k]? Or replace [r], complete [c], or keep only [Press *Enter*] record {}? (Type column name to inspect or 'help' for options): ".format(old[0]))
    while replace.lower() in columns or replace.lower()=='help':
      if replace.lower() in columns: print at.Table(np.asarray([[i for idx,i in enumerate(pold) if idx in [0,columns.index(replace.lower())]],[i for idx,i in enumerate(pnew) if idx in [0,columns.index(replace.lower())]]]), names=['id',replace.lower()])    
      elif replace.lower()=='help': print at.Table(np.asarray([['-'*30,'-'*100],['[column name]','Display full record entry for that column without taking action'],['k','Keep both records and assign second one new id if necessary'],['r','Replace all columns of first record with second record values'],['r [column name] [column name]...','Replace specified columns of first record with second record values'],['c','Complete empty columns of first record with second record values where possible'],['[Enter]','Keep first record and delete second'],['abort','Abort merge of current table, undo all changes, and proceed to next table']]), names=['Command','Result'])
      replace = raw_input("Keep both records [k]? Or replace [r], complete [c], or keep only [Press *Enter*] record {}? (Type column name to inspect or 'help' for options): ".format(old[0]))

    if replace and (all([i in list(columns)+options for i in replace.lower().split()]) or replace.lower().startswith('sql')):
      
      # Replace the old record with the new record
      if replace.lower().startswith('r') and 'r' in options:
        if replace.lower()=='r':
          sure = raw_input('Are you sure you want to replace record {} with record {}? [y/n] : '.format(old[0],new[0]))
          if sure.lower()=='y':
            empty_cols, new_vals = zip(*[['{}=?'.format(e),n] for e,n in zip(columns[1:],new[1:])])
            if delete: self.modify("DELETE FROM {} WHERE id={}".format(table, new[0]))
            self.modify("UPDATE {} SET {} WHERE id={}".format(table, ','.join(empty_cols), old[0]), tuple(new_vals))
        elif all([i in list(columns)+options for i in replace.lower().split()]):
          empty_cols, new_vals = zip(*[['{}=?'.format(e),n] for e,n in zip(columns[1:],new[1:]) if e in replace])
          if empty_cols:
            if delete: self.modify("DELETE FROM {} WHERE id={}".format(table, new[0]))
            self.modify("UPDATE {} SET {} WHERE id={}".format(table, ','.join(empty_cols), old[0]), tuple(new_vals))
      
      # Complete the old record with any missing data provided in the new record, then delete the new record
      elif replace.lower()=='c' and 'c' in options:
        try:
          empty_cols, new_vals = zip(*[['{}=?'.format(e),n] for e,o,n in zip(columns[1:],old[1:],new[1:]) if repr(o).lower() in ['','none','null'] and repr(n).lower() not in ['','none','null']])
          if delete: self.modify("DELETE FROM {} WHERE id={}".format(table, new[0]))
          self.modify("UPDATE {} SET {} WHERE id={}".format(table, ','.join(empty_cols), old[0]), tuple(new_vals))
        except:
          if delete: self.modify("DELETE FROM {} WHERE id={}".format(table, new[0]))
          else: pass
      elif replace.lower()=='k' and 'k' in options: return [old[0],new[0]]
      
      # Execute raw SQL
      elif replace.lower().startswith('sql ') and 'sql' in options: 
        try: self.modify(replace[4:])
        except: pass 
    
    # Abort the current database clean up
    elif replace.lower()=='abort': return 'abort'
    
    # Delete the current record if *delete=True
    elif not replace:
      if delete: self.modify("DELETE FROM {} WHERE id={}".format(table, new[0]))
      else: pass
    
  def header(self, spectrum_id_or_path):
    """
    Prints the header information for the given **spectrum_id_or_path**.
    
    Parameters
    ----------
    spectrum_id_or_path: (int, str)
      The id from the SPECTRA table or the path to the FITS file of the spectrum header to print.
    
    """
    if isinstance(spectrum_id_or_path,int):
      try: 
        H = self.dict("SELECT * FROM spectra WHERE id={}".format(spectrum_id_or_path)).fetchone()['header']
        if H: return H
        else: print 'No header for spectrum {}'.format(spectrum_id_or_path)
      except TypeError: print 'No spectrum with id {}'.format(spectrum_id_or_path)
    elif os.path.isfile(spectrum_id_or_path):
      if spectrum_id_or_path.endswith('.fits'):
        return clean_header(spectrum_id_or_path)
      else:
        txt, H = open(spectrum_id_or_path), []
        for i in txt: 
          if i.startswith('#'): H.append(i)
        txt.close()
        print ''.join(H) if H else 'No header for spectrum {}'.format(spectrum_id_or_path)
    else: print 'No such file {}'.format(spectrum_id_or_path)

  def identify(self, search):
    """
    For **search** input of (ra,dec) decimal degree tuple, i.e. '(12.3456,-65.4321)', returns all sources within 1 arcminute.
    For **search** input of text string, i.e. 'vb10', returns all sources with case-insensitive partial text matches in *names* or *designation* columns.
    
    Parameters
    ----------
    search: (str, tuple)
      The text or coordinate tuple to search the SOURCES table with.
      
    """
    try: q = "SELECT id,ra,dec,designation,unum,shortname,names FROM sources WHERE ra BETWEEN "+str(search[0]-0.01667)+" AND "+str(search[0]+0.01667)+" AND dec BETWEEN "+str(search[1]-0.01667)+" AND "+str(search[1]+0.01667)
    except TypeError: q = "SELECT id,ra,dec,designation,unum,shortname,components,companions,names FROM sources WHERE REPLACE(names,' ','') like '%"+search.replace(' ','')+"%' or designation like '%"+search+"%' or unum like '%"+search+"%' or shortname like '%"+search+"%'"
    results = self.query(q)
    if results: 
      if len(results)==1: self.inventory(int(results[0][0]))
      else: print at.Table(np.asarray(results), names=['id','ra','dec','designation','unum','short','components','companions','names'])
    else: print "No objects found by {}".format(search)
      
  def inventory(self, source_id, plot=False, return_data=False):
    """
    Prints a summary of all objects in the database. Input string or list of strings in **ID** or **unum** for specific objects.
    
    Parameters
    ----------
    source_id: int
      The id from the SOURCES table whose data across all tables is to be printed.
    plot: bool
      Plots all spectra for the object.
    return_data: bool
      Return the results.
      
    Returns
    -------
    data_tables: dict
      Returns a dictionary of astropy tables with the table name as the keys.
    
    """
    data_tables = {}
    try:
      for table in ['sources']+[t for t in self.query("SELECT * FROM sqlite_master WHERE type='table'", unpack=True)[1] if t!='sources']:
        
        # Get the columns, pull out redundant ones, and query the table for this source's data
        columns, types = self.query("PRAGMA table_info({})".format(table), unpack=True)[1:3]
        columns, types = [z[np.where(np.logical_and(columns!='spectrum',columns!='source_id'))] for z in [columns,types]]
        try: data = self.query("SELECT {} FROM {} WHERE {}".format(','.join(columns),table,'id={}'.format(source_id) if table=='sources' else 'source_id={}'.format(source_id)))
        except: data = np.array([])
        
        # If there's data for this table, save it
        if data.any(): 
          data_table = at.Table(np.asarray(data), names=columns)
          data_tables[table] = data_table
          print ''; print table.upper(); print data_table
        
      if plot:
        for i in self.query("SELECT id FROM spectra WHERE source_id={}".format(source_id), unpack=True)[0]: self.plot_spectrum(i)
    
    except: print 'No source with id {}. Try db.identify() to search the database for a source_id.'.format(source_id)
    
    if return_data: return data_tables

  def lookup(self, table, ids=None, concatenate='', delim='/'):
    """
    Quickly look up records from the specified *table* and list *ids* to limit results. Specify column values to *concatenate* into a string.
    
    Parameters
    ----------
    table: str
      The name of the table to perform a lookup on.
    ids: (int, list)
      An id or list of ids to lookup.
    concatenate: str
      The name of the column whose values should be joined by delim and returned.
    delim: str
      If concatenate, the delimiter to be used to join the results.
      
    Returns
    -------
    (list, str)
      A list of the complete record(s) or a concatenated string from the desired table with the given ids. 
    
    """
    if ids=='-' or ids==['-']: return '-'
    elif type(ids)==str and table.lower()=='publications' and ',' not in ids:
      return self.list("SELECT * FROM publications WHERE shortname LIKE '%{0}%' OR description LIKE '%{0}%' OR bibtex LIKE '%{0}%'".format(ids)).fetchall()
    else:
      if type(ids)==int: ids = [ids]
      if type(ids)==str and ',' in ids: ids = ids.split(',')
      if type(ids)==list and not ids: ids = ''
      try: results = self.list("SELECT {} FROM {} WHERE id IN ({})".format(concatenate or '*',table, ','.join(map(str,ids)))).fetchall() if ids else self.list("SELECT * FROM {}".format(table)).fetchall()
      except: results = ''
      return '' if not results else delim.join(map(str,zip(*results)[0])) if concatenate else results

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
      con.list("ATTACH DATABASE '{}' AS m".format(master)), self.list("ATTACH DATABASE '{}' AS c".format(conflicted)), con.list("ATTACH DATABASE '{}' AS c".format(conflicted)), self.list("ATTACH DATABASE '{}' AS m".format(master))
      
      # Drop any backup tables from failed merges
      for table in tables: self.list("DROP TABLE IF EXISTS Backup_{0}".format(table))
      
      # Gather user data to add to CHANGELOG table
      import socket, datetime
      user, machine_name, date, modified_tables = raw_input('Please enter your name : '), socket.gethostname(), datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), []
      
      # Print instructions for user
      print at.Table(np.asarray([['-'*30,'-'*100],['[column name]','Display full record entry for that column without taking action'],['k','Keeps both records and assigns second one new id if necessary'],['r','Replaces all columns of first record with second record values'],['r [column name] [column name]...','Replaces specified columns of first record with second record values'],['c','Complete empty columns of first record with second record values where possible'],['[Enter]','Keep first record and delete second'],['abort','Abort merge of current table, undo all changes, and proceed to next table']]), names=['Command','Result'])
      
      # Merge table by table, starting with SOURCES
      tables = tables or ['sources']+[t for t in zip(*self.list("SELECT * FROM sqlite_master WHERE name NOT LIKE '%Backup%' AND type='table'{}".format(" AND name IN ({})".format("'"+"','".join(tables)+"'") if tables else '')).fetchall())[1] if t!='sources']
      for table in tables:
        # Get column names and data types from master table and column names from conflicted table
        (columns, types), conflicted_cols = zip(*self.list("PRAGMA table_info({})".format(table)).fetchall())[1:3], zip(*con.list("PRAGMA table_info({})".format(table)).fetchall())[1]
                
        if any([i not in columns for i in conflicted_cols]):
          # Abort table merge if conflicted has new columns not present in master. New columns must be added to the master database first via db.edit_columns().
          print "\nMerge of {0} table aborted since conflicted copy has columns {1} not present in master.\nAdd new columns to master with BDdb.edit_columns() and try again.\n".format(table.upper(),[i for i in conflicted_cols if i not in columns])
        
        else:
          # Add new columns from master table to conflicted table if necessary
          if any([i not in conflicted_cols for i in columns]): con.modify("DROP TABLE IF EXISTS Conflicted_{0}".format(table)), con.modify("ALTER TABLE {0} RENAME TO Conflicted_{0}".format(table)), con.modify("CREATE TABLE {0} ({1})".format(table, ', '.join(['{} {}'.format(c,t) for c,t in zip(columns,types)]))), con.modify("INSERT INTO {0} ({1}) SELECT {1} FROM Conflicted_{0}".format(table, ','.join(conflicted_cols))), con.modify("DROP TABLE Conflicted_{0}".format(table))
        
          # Pull unique records from conflicted table
          data = map(list, con.list("SELECT * FROM (SELECT 1 AS db, {0} FROM m.{2} UNION ALL SELECT 2 AS db, {0} FROM c.{2}) GROUP BY {1} HAVING COUNT(*)=1 AND db=2".format(','.join(columns),','.join(columns[1:]),table)).fetchall())

          if data:
            if diff_only:
              print at.Table(np.asarray([[repr(i) for i in d] for d in data]), names=columns)
            else:
              # Make temporary table copy so changes can be undone at any time
              self.modify("DROP TABLE IF EXISTS Backup_{0}".format(table))
              self.modify("ALTER TABLE {0} RENAME TO Backup_{0}".format(table))
              self.modify("CREATE TABLE {0} ({1})".format(table, ', '.join(['{} {}'.format(c,t) for c,t in zip(columns,types)])))
              self.modify("INSERT INTO {0} ({1}) SELECT {1} FROM Backup_{0}".format(table, ','.join(columns)))

              # Create a dictionary of any reassigned ids from merged SOURCES tables and replace applicable source_ids in other tables.
              print "\nMerging {} tables.\n".format(table.upper())
              try: count = self.list("SELECT MAX(id) FROM {}".format(table)).fetchone()[0]+1
              except TypeError: count = 1
              for n,i in enumerate([d[1:] for d in data]):
                if table=='sources': reassign[i[0]] = count
                elif 'source_id' in columns and i[1] in reassign.keys(): i[1] = reassign[i[1]]
                else: pass
                i[0] = count
                data[n] = i
                count += 1
            
              # Insert unique conflicted records into master and run BDdb.clean_up()
              self.modify("INSERT INTO {} VALUES({})".format(table, ','.join(['?' for c in columns])), data)
              print "{} records added to {} table at '{}':".format(len(data), table, master)
              print at.Table(np.asarray([[repr(i) for i in d] for d in data]), names=columns)
              abort = self.clean_up(table)
          
              # Undo all changes to table if merge is aborted. Otherwise, push table changes to master.
              if abort: self.modify("DROP TABLE {0}".format(table)), self.modify("ALTER TABLE Backup_{0} RENAME TO {0}".format(table))
              else: self.modify("DROP TABLE Backup_{0}".format(table)), modified_tables.append(table.upper())
          
          else: print "{} tables identical.".format(table.upper())
      
      # Add data to CHANGELOG table
      if not diff_only:
        user_description = raw_input('\nPlease describe the changes made in this merge : ')
        self.modify("INSERT INTO changelog VALUES(?, ?, ?, ?, ?, ?, ?)", (None, date, user, machine_name, ', '.join(modified_tables), user_description, os.path.basename(conflicted)))
      
      # Finish up and detach
      print "\nMerge complete!" if not diff_only else "\nDiff complete. No changes made to either database."
      con.modify("DETACH DATABASE c"), self.modify("DETACH DATABASE c"), con.modify("DETACH DATABASE m"), self.modify("DETACH DATABASE m"), con.modify.close()
    else: print "File '{}' not found!".format(conflicted)

  def modify(self, SQL, params='', verbose=False):
    """
    Wrapper for CRUD operations to make them distinct from queries and automatically pass commit() method to cursor.
    
    Parameters
    ----------
    SQL: str
      The SQL query to execute
    params: sequence
      Mimicks the native parameter substitution of sqlite3
    verbose: bool
      Print the number of modified records
    """
    try:
      if SQL.lower().startswith('select'):
        print 'Use self.query method for queries.'
      else:
        self.list(SQL, params)
        self.conn.commit()
        if verbose: print 'Number of records modified: {}'.format(self.list("SELECT changes()").fetchone()[0])
    except IOError:
      print 'Could not execute! Please check the query syntax and parameters format.'

  def modify_table(self, table, columns, types, new_table=False):
    """
    Rearrange, add or delete columns from database **table** with desired ordered list of **columns** and corresponding data **types**.
    
    Parameters
    ----------
    table: str
      The name of the table to modify
    columns: list
      A list of the columns in the order in which they are to appear in the SQL table
    types: list
      A list of the types corresponding to each column in the columns list above.
    new_table: bool
      Create a new table
    
    """
    # Make sure there is an 'id' and a 'source_id' column
    goodtogo = True
    if columns[0]!='id' or types[0].upper()!='INTEGER PRIMARY KEY': print "Column 1 must be 'id' with type 'INTEGER PRIMARY KEY'"; goodtogo = False
    if table!='sources' and (columns[1]!='source_id' or types[1].upper()!='INTEGER'): print "Column 2 must be 'source_id' with type 'INTEGER'"; goodtogo = False
    
    if goodtogo:
      # If the table exists, modify the columns
      if table in zip(*self.list("SELECT name FROM sqlite_master").fetchall())[0] and not new_table:
        self.list("ALTER TABLE {0} RENAME TO TempOldTable".format(table))
        self.list("CREATE TABLE {0} ({1})".format(table, ', '.join(['{} {}'.format(c,t) for c,t in zip(columns,types)])))
        self.list("INSERT INTO {0} ({1}) SELECT {1} FROM TempOldTable".format(table, ','.join([c for c in list(zip(*self.list("PRAGMA table_info(TempOldTable)").fetchall())[1]) if c in columns])))
        self.list("DROP TABLE TempOldTable")
    
      # If the table does not exist and new_tabe is True, create it
      elif table not in zip(*self.list("SELECT name FROM sqlite_master").fetchall())[0] and new_table:
        self.list("CREATE TABLE {0} ({1})".format(table, ', '.join(['{} {}'.format(c,t) for c,t in zip(columns,types)])))
    
      # Otherwise the table to be modified doesn't exist or the new table to add already exists, so do nothing
      else:
        if new_table: print 'Table {} already exists. Set *new_table=False to modify.'.format(table.upper())
        else: print 'Table {} does not exist. Could not modify. Set *new_table=True to add a new table.'.format(table.upper())
        
    else: print 'The {} table has not been modified. Please make sure your table columns and types are modified properly.'.format(table.upper())
    
  def output_spectrum(self, spectrum_id, filepath):
    """
    Prints a file of the spectrum with id **spectrum_id** to an ascii file with specified **filepath**.
    
    Parameters
    ----------
    spectrum_id: int
      The id from the SPECTRA table of the spectrum to print to file.
    filepath: str
      The path of the file to print the data to.
    
    """
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
    i = self.dict("SELECT * FROM spectra WHERE id={}".format(spectrum_id)).fetchone()
    if i:
      try:
        if not overplot: 
          plt.figure()
          plt.rc('text', usetex=False)
          plt.grid(True)
          plt.yscale('log', nonposy='clip')
          plt.title('source_id = {}'.format(i['source_id']))
          plt.figtext(0.15,0.88, '{}\n{}\n{}\n{}'.format(i['filename'],self.list("SELECT name FROM telescopes WHERE id={}".format(i['telescope_id'])).fetchone()[0] if i['telescope_id'] else '',self.list("SELECT name FROM instruments WHERE id={}".format(i['instrument_id'])).fetchone()[0] if i['instrument_id'] else '',i['obs_date']), verticalalignment='top')
          plt.xlabel('[{}]'.format(i['wavelength_units']))
          plt.ylabel('[{}]'.format(i['flux_units']))
          plt.legend(loc=8, frameon=False)
        wav, flx, err = i['spectrum']
        plt.loglog(wav, flx, c=color, label='spec_id: {}'.format(i['id']))
        X, Y = plt.xlim(), plt.ylim()
        try: plt.fill_between(wav, flx-err, flx+err, color='b', alpha=0.3), plt.xlim(X), plt.ylim(Y)
        except TypeError: print 'No uncertainty array for spectrum {}'.format(spectrum_id)
      except: print "Could not plot spectrum {}".format(spectrum_id)
    else: print "No spectrum {} in the SPECTRA table.".format(ID)

  def query(self, SQL, params='', fmt='array', fetch='all', unpack=False):
    """
    Wrapper for cursors so data can be retrieved as a list or dictionary from same method
    
    Parameters
    ----------
    SQL: str
      The SQL query to execute
    params: sequence
      Mimicks the native parameter substitution of sqlite3
    fmt: str
      Returns the data as a dictionary, list, array, or astropy.table given 'dict', 'list', 'array', 'table'
    unpack: bool
      Returns the transpose of the data
      
    Returns
    -------
    result: (array,dict,list,table)
      The result of the database query
    """
    try:
      # Restricy queries to SELECT and PRAGMA statements
      if SQL.lower().startswith('select') or SQL.lower().startswith('pragma'):
        
        # Return a dictionary or astropy table
        if fmt.lower() in ['dict','table']: 
          # Get the results as a dictionary
          result = self.dict(SQL, params).fetchone() if fetch=='one' else self.dict(SQL, params).fetchall()
          
          # Make an astropy table if necessary
          if fmt.lower()=='table': result = at.Table(result)
        
        # Return an array or list
        else: 
          # Get the results as an array by default
          result = np.asarray(self.list(SQL, params).fetchone() if fetch=='one' else self.list(SQL, params).fetchall())
          
          # Unpack the results if necessary
          if unpack: result = result.T
          
          # Turn back into list if necessary
          if fmt.lower()=='list': result = list(result)
        
        return result
          
      else:
        print 'Queries must begin with a SELECT or PRAGMA statement. For database modifications use self.modify() method.'  
    except IOError:
      print 'Could not execute! Please check the query syntax and parameters format.'

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

def adapt_header(header): 
  """
  Adapts a FITS header into a HEADER string to put into the database.
  
  Parameters
  ----------
  header: fits.header
    The FITS header to be adapted into a HEADER type that can be inserted into a SQL file.
    
  Returns
  -------
  HEADER
    The adapted HEADER object
    
  """
  return header.tostring(sep='\n')

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

def convert_header(header):
  """
  Converts a HEADER string stored in the database back into a FITS header.
  
  Parameters
  ----------
  header: HEADER
    The header object to be converted back into a FITS header.
    
  Returns
  -------
  fits.header
    The converted FITS header.
    
  """
  return pf.Header().fromstring(header, sep='\n') if header else None

def convert_spectrum(url):
  """
  Converts a URL string stored in the database into a (W,F,E) sequence of arrays.
  
  Parameters
  ----------
  url: SPECTRUM
    The CUNY Academic Works or AWS S3 URL of the file to be converted into arrays.
    
  Returns
  -------
  sequence
    The converted spectrum.
    
  """
  spectrum = ''
  
  # If it is a url, retrieve the data
  if isinstance(url,str):
    if url.endswith('.fits'):
      try: spectrum = pf.open(url, cache=True)[0].data
      except: pass
    if url.endswith('.txt'): 
      try: spectrum = np.genfromtxt(urllib2.urlopen(url), unpack=True)
      except: pass
    else:
      pass
  
  # THIS BIT IS JUST IN CASE THE DATA IS NOT UPLOADED TO A SERVER
  # If it is a [W,F,E] spectrum, just return it
  elif isinstance(url, (list,np.ndarray,tuple)): 
    try: spectrum = url
    except: pass

  if spectrum!='': return spectrum
  else: print 'Could not retrieve spectrum at {}. Are you connected to the Internet?'.format(url)

# Register the adapters
sql.register_adapter(np.ndarray, adapt_array)
sql.register_adapter(pf.header.Header, adapt_header)

# Register the converters
sql.register_converter("ARRAY", convert_array)
sql.register_converter("HEADER", convert_header)
sql.register_converter("SPECTRUM", convert_spectrum)
