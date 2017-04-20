#!/usr/bin/python

# David Rodriguez
# Created: November 4, 2015
# Read in a table from an SQL query of the BDNYC Database and convert to a VOTable file for use elsewhere

from astropy.table import Table, Column
from astropy.io.votable import from_table


def dict_tovot(tabdata, tabname='votable.xml', phot=False, binary=True):
    """
    Converts dictionary table **tabdata** to a VOTable with name **tabname**

    Parameters
    ----------
    tabdata: list
      SQL query dictionary list from running query_dict.execute()
    tabname: str
      The name of the VOTable to be created
    phot: bool
      Parameter specifying if the table contains photometry to be merged
    binary: bool
      Parameter specifying if the VOTable should be saved as a binary.
      This is necessary for tables with lots of text columns.

    """

    # Check if input is a dictionary
    if not isinstance(tabdata[0], dict):
        raise TypeError('Table must be a dictionary. Call the SQL query with query_dict.execute()')

    # Create an empty table to store the data
    t = Table()

    colnames = tabdata[0].keys()

    # If this is a photometry table, parse it and make sure to have the full list of columns
    if phot:
        tabdata = photparse(tabdata)

        colnames = tabdata[0].keys()

        for i in range(len(tabdata)):
            tmpcol = tabdata[i].keys()
            for elem in tmpcol:
                if elem not in colnames:
                    colnames.append(elem)

        # No need for band column any more
        try:
            colnames.remove('band')
        except ValueError:
            pass

    # Run through all the columns and create them
    for elem in colnames:
        table_add(t, tabdata, elem)

    # Output to a file
    print('Creating table...')
    votable = from_table(t)

    # Required in some cases (ie, for lots of text columns)
    if binary:
        votable.set_all_tables_format('binary')

    votable.to_xml(tabname)

    print('Table created: {}'.format(tabname))


def photaddline(tab, sourceid):
    """
    Loop through the dictionary list **tab** creating a line for the source specified in **sourceid**

    Parameters
    ----------
    tab:
      Dictionary list of all the photometry data
    sourceid:
      ID of source in the photometry table (source_id)

    Returns
    -------
    tmpdict: dict
        Dictionary with all the data for the specified source

    """

    colnames = tab[0].keys()
    tmpdict = dict()
    for i in range(len(tab)):

        # If not working on the same source, continue
        if tab[i]['source_id'] != sourceid:
            continue

        # Check column names and create new ones for band-specific ones
        for elem in colnames:
            if elem not in ['comments', 'epoch', 'instrument_id', 'magnitude', 'magnitude_unc', 'publication_id',
                            'system', 'telescope_id']:
                tmpdict[elem] = tab[i][elem]
            elif elem == 'band':
                continue
            else:
                tmpstr = tab[i]['band']+'.'+elem
                tmpdict[tmpstr] = tab[i][elem]

    return tmpdict


def photparse(tab):
    """
    Parse through a photometry table to group by source_id

    Parameters
    ----------
    tab: list
      SQL query dictionary list from running query_dict.execute()

    Returns
    -------
    newtab: list
      Dictionary list after parsing to group together sources

    """

    # Check that source_id column is present
    if 'source_id' not in tab[0].keys():
        raise KeyError('phot=TRUE requires the source_id columb be included')

    # Loop through the table and grab unique band names and source IDs
    uniqueid = []
    for i in range(len(tab)):
        tmpid = tab[i]['source_id']

        if tmpid not in uniqueid:
            uniqueid.append(tmpid)

    # Loop over unique id and create a new table for each element in it
    newtab = []
    for sourceid in uniqueid:
        tmpdict = photaddline(tab, sourceid)
        newtab.append(tmpdict)

    return newtab


def table_add(tab, data, col):
    """
    Function to parse dictionary list **data** and add the data to table **tab** for column **col**

    Parameters
    ----------
    tab: Table class
      Table to store values
    data: list
      Dictionary list from the SQL query
    col: str
      Column name (ie, dictionary key) for the column to add

    """

    x = []
    for i in range(len(data)):

        # If the particular key is not present, use a place-holder value (used for photometry tables)
        if col not in data[i]:
            temp = ''
        else:
            temp = data[i][col]

        # Fix up None elements
        if temp is None: temp = ''

        x.append(temp)

    print('Adding column {}'.format(col))
    tab.add_column(Column(x, name=col))
