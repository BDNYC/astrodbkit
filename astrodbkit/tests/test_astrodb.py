import pytest
import tempfile
import os
from astropy.utils.data import download_file
from .. import astrodb


def setup_module(module):
    db_path = download_file("http://github.com/BDNYC/BDNYCdb/raw/master/BDNYCv1.0.db")
    module.bdnyc_db = astrodb.Database(db_path)
    filename = os.path.join(tempfile.mkdtemp(), 'empty_db.db')
    module.empty_db = astrodb.create_database(filename)

def test_loaddb():
    print(bdnyc_db)
    print(empty_db)

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
    
@pytest.mark.xfail
def test_add_data():
    assert False
    
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

@pytest.mark.xfail
def test_table():
    assert False