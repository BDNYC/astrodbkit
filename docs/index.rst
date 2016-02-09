.. astrodbkit documentation master file, created by
   sphinx-quickstart on Tue Jan 19 10:54:25 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to astrodbkit's documentation!
======================================

The BDNYC Database is an advanced SQL relational database of published spectra, photometry and astrometry for over 1300 very low-mass stars, brown dwarfs and planetary mass objects.

This documentation describes a tool kit of classes, methods and functions useful for CRUD operations and analysis of data from The BDNYC Data Archive.

Getting Started
===============

To install, just do::

    pip install astrodbkit

This package includes the initial release of the BDNYC Database, which contains the astrometry, photometry and spectra for the 198 objects in the `Filippazzo et al. (2015)`_ sample.

.. _Filippazzo et al. (2015): http://adslabs.org/adsabs/abs/2015ApJ...810..158F/

.. note:: For access to the full dataset, an email request must be made to a BDNYC group admin.

Accessing the Database
======================

To start using the database, launch iPython, import the module, then initialize the database with the :class:`astrodb.get_db()` class like so::

    from astrodbkit import astrodb
    db = astrodb.get_db()


Voila! You can see an inventory of all data for a specific source by passing a *source_id* to the :py:meth:`~astrodb.get_db.inventory` method::

    db.inventory(86)

This will also plot all available spectra for that source for visual inspection if you set **plot=True**.

Querying the Database
=====================

Now that you have the database at your fingertips, you’ll want to get some info out of it. To do this, you can pass SQL queries wrapped in double-quotes (") to the :py:meth:`~astrodb.get_db.query` method::

    data = db.query( "SQL_query_goes_here" )

`Here is a detailed post about how to write a SQL query`_.

.. _Here is a detailed post about how to write a SQL query: http://www.bdnyc.org/?p=898

Contents
========

.. toctree::
   :maxdepth: 2

   astrodbkit.astrodb
   astrodbkit.votools

Indices and tables

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

