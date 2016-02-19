.. astrodbkit documentation master file, created by
   sphinx-quickstart on Tue Jan 19 10:54:25 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to astrodbkit's documentation!
======================================

This documentation describes a tool kit of classes, methods and functions useful for CRUD operations and analysis of data from an SQL database.

Getting Started
===============

To install, do::

    pip install astrodbkit

Creating a Database
===================

To create a database from scratch, do::

    dbpath = '/desired/path/to/my_new_database.db'
    astrodb.create_database(dbpath)
    
Alternatively, you can `download and use the BDNYC Database`_, which contains the astrometry, photometry and spectra for the 198 objects in the `Filippazzo et al. (2015)`_ sample.

.. _download and use the BDNYC Database: https://github.com/BDNYC/BDNYCdb
.. _Filippazzo et al. (2015): http://adslabs.org/adsabs/abs/2015ApJ...810..158F/

.. note:: For access to the full dataset, an email request must be made to a BDNYC group admin.

Accessing the Database
======================

To start using the database, launch iPython, import the module, then initialize the database with the :class:`astrodb.get_db()` class like so::

    db = astrodb.get_db(dbpath)

Voila! 

Querying the Database
=====================

Now that you have the database at your fingertips, you’ll want to get some info out of it. 

You can see an inventory of all data for a specific source by passing an integer id to the :py:meth:`~astrodb.get_db.inventory` method::

    data = db.inventory(86)

This will retrieve the data across all tables with the specified source_id for visual inspection. Setting *fetch=True* will return the data as a dictionary of Astropy tables so that table and column keys can be used to access the results. For example::

    data['photometry'][['band','magnitude','magnitude_unc']]

will return a table of the band, magnitude and uncertainty for all records in the sources table with that source_id.

You can search any table in the database with the :py:meth:`~astrodb.get_db.identify` method by supplying a string, integer, or (ra,dec) coordinates along with the table to search. For example, if I want to find all the records in the SOURCES table in the HR 8799 system::

    db.search('8799', 'sources')
    
Or all the papers published by Joe Filippazzo::

    db.search('Fili', 'publications')

You can also pass SQL queries wrapped in double-quotes (") to the :py:meth:`~astrodb.get_db.query` method::

    data = db.query( "SQL_query_goes_here" )

`Here is a detailed post about how to write a SQL query`_.

.. _Here is a detailed post about how to write a SQL query: http://www.bdnyc.org/?p=898

Contents
========

.. toctree::
   :maxdepth: 2

   astrodb
   votools

Indices and tables

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

