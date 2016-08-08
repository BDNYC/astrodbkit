import pytest
import tempfile
import os
from astropy.utils.data import download_file
from .. import astrodb
from sqlite3 import IntegrityError


def setup_module(module):
    try:
        db_path = download_file("http://github.com/BDNYC/BDNYCdb/raw/master/bdnyc_database.db")
    except:
        db_path = download_file("http://github.com/BDNYC/BDNYCdb/raw/master/BDNYCv1.0.db")
    module.bdnyc_db = astrodb.Database(db_path)
    filename = os.path.join(tempfile.mkdtemp(), 'empty_db.db')
    astrodb.create_database(filename)
    module.empty_db = astrodb.Database(filename)


def test_load_bdnyc():
    print(bdnyc_db)
    assert not isinstance(bdnyc_db, type(None))


def test_load_empty():
    print(empty_db)
    assert not isinstance(empty_db, type(None))


def test_search(): 
    bdnyc_db.search('1234','sources')


def test_inventory():
    bdnyc_db.inventory(825)


def test_sqlquery():
    bdnyc_db.query("SELECT s.id, s.ra, s.dec, s.shortname, p.source_id, p.band, p.magnitude "
        "FROM sources as s JOIN photometry as p ON s.id=p.source_id "
        "WHERE s.dec<=-10 AND p.band=='W1'")


def test_schema():
    bdnyc_db.schema('sources')


def test_table():
    columns = ['id', 'ra', 'dec', 'shortname', 'source_id']
    types = ['INTEGER', 'REAL', 'REAL', 'TEXT', 'INTEGER']
    constraints = ['NOT NULL UNIQUE', '', '', '', '']
    empty_db.table('new_sources', columns, types, constraints, new_table=True)


def test_add_data():
    data = list()
    data.append(['ra', 'dec', 'shortname'])
    data.append([12, -12, 'fakesource'])
    empty_db.add_data(data, 'sources')


def test_add_foreign_key():
    empty_db.add_foreign_key('new_sources', ['sources'], ['source_id'], ['id'])


def test_foreign_key_support():
    data = list()
    data.append(['ra', 'dec', 'shortname', 'source_id'])
    data.append([12, -12, 'fakesource', 9999])
    with pytest.raises(IntegrityError):
        empty_db.add_data(data, 'new_sources')  # Foreign key error


def test_lookup():
    bdnyc_db.lookup([1, '2MASS'], 'sources')


@pytest.mark.xfail
def test_clean_up():
    assert False

@pytest.mark.xfail
def test_merge():
    assert False

@pytest.mark.xfail
def test_output_spectrum():
    assert False

@pytest.mark.xfail
def test_plot_spectrum():
    assert False
