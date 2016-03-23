from astropy.utils.data import download_file
from .. import astrodb

def setup_module(module):
    db_path = download_file("http://github.com/BDNYC/BDNYCdb/raw/master/BDNYCv1.0.db")
    module.bdnyc_db = astrodb.Database(db_path)
    filename = tmpdir.join('testdb.db').strpath
    empty_db = astrodb.create_database(filename)

def test_newdb(tmpdir):
    print(empty_db)

def test_loaddb():
    print(bdnyc_db)

def test_search(): 
    db.search('1234','sources')

def test_inventory():
    db.inventory(825)

def test_sqlquery():
    db.query("SELECT s.id, s.ra, s.dec, s.shortname, p.source_id, p.band, p.magnitude "
        "FROM sources as s JOIN photometry as p ON s.id=p.source_id "
        "WHERE s.dec<=-10 AND p.band=='W1'")

def test_schema():
    db.schema('sources')
    
def test_add_data():
    assert False

def test_clean_up():
    assert False

def test_merge():
    assert False

def test_output_spectrum():
    assert False

def test_plot_spectrum():
    assert False

def test_table():
    assert False