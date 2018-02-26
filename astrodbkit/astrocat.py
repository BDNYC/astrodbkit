"""
Module to ingest arbitrary VizieR catalogs into a custom cross-matched database

Authors: Joe Filippazzo, Andrea Lin
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import astropy.units as q
import astropy.table as at
import astropy.coordinates as coord
import datetime
from sklearn.cluster import DBSCAN
from collections import Counter
from scipy.stats import norm
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch
from sklearn.externals import joblib
from astropy.coordinates import SkyCoord

Vizier.ROW_LIMIT = -1

class Catalog(object):
    
    def __init__(self, name='Test'):
        """
        Initialize a catalog object
        
        Parameters
        ----------
        name: str
            The name of the database
        """
        self.name = name
        self.catalog = pd.DataFrame(columns=('source_id','ra','dec','flag','cat_name','catID'))
        self.n_sources = 0
        self.history = "{}: Database created".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        self.catalogs = {}
        self.xmatch_radius = 0.0001
        
    @property
    def info(self):
        """
        Print the history
        """
        print(self.history)
        
    def ingest_data(self, data, cat_name, id_col, ra_col='_RAJ2000', dec_col='_DEJ2000', count=-1):
        """
        Ingest a data file and regroup sources
        
        Parameters
        ----------
        data: str, pandas.DataFrame, astropy.table.Table
            The path to the exported VizieR data or the data table
        cat_name: str
            The name of the added catalog
        id_col: str
            The name of the column containing the unique ids
        ra_col: str
            The name of the RA column
        dec_col: str
            The name of the DEC column
        count: int
            The number of table rows to add
            (This is mainly for testing purposes)
        """
        # Check if the catalog is already ingested
        if cat_name in self.catalogs:
            
            print('Catalog {} already ingested.'.format(cat_name))
            
        else:
            
            if isinstance(data, str):
                path = data
                data = pd.read_csv(data, sep='\t', comment='#', engine='python')[:count]
                
            elif isinstance(data, pd.core.frame.DataFrame):
                path = type(data)
                
            elif isinstance(data, (at.QTable, at.Table)):
                path = type(data)
                data = pd.DataFrame(list(data), columns=data.colnames)
                
            else:
                print("Sorry, but I cannot read that data. Try an ascii file path, astropy table, or pandas data frame.")
                return
                
            # Make sure ra and dec are decimal degrees
            if isinstance(data[ra_col][0], str):
                
                crds = coord.SkyCoord(ra=data[ra_col], dec=data[dec_col], unit=(q.hour, q.deg), frame='icrs')
                data.insert(0,'dec', crds.dec)
                data.insert(0,'ra', crds.ra)
                
            elif isinstance(data[ra_col][0], float):
                
                data.rename(columns={ra_col:'ra', dec_col:'dec'}, inplace=True)
            
            else:
                print("I can't read the RA and DEC of the input data. Please try again.")
                return
                
            # Change some names
            data.insert(0,'catID', ['{}_{}'.format(cat_name,n+1) for n in range(len(data))])
            data.insert(0,'dec_corr', data['dec'])
            data.insert(0,'ra_corr', data['ra'])
            data.insert(0,'source_id', np.nan)
            
            print('Ingesting {} rows from {} catalog...'.format(len(data),cat_name))
            
            # Save the raw data as an attribute
            setattr(self, cat_name, data)
                
            # Update the history
            self.history += "\n{}: Catalog {} ingested.".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),cat_name)
            self.catalogs.update({cat_name:(path,id_col)})
            
    def inventory(self, source_id):
        """
        Look at the inventory for a given source
        
        Parameters
        ----------
        source_id: int
            The id of the source to inspect
        """
        if self.n_sources==0:
            print('Please run group_sources() to create the catalog first.')
        
        else:
            
            if source_id>self.n_sources or source_id<1 or not isinstance(source_id, int):
                print('Please enter an integer between 1 and',self.n_sources)
            
            else:
            
                print('Source:')
                print(at.Table.from_pandas(self.catalog[self.catalog['id']==source_id]).pprint())
                for cat_name in self.catalogs:
                    cat = getattr(self, cat_name)
                    rows = cat[cat['source_id']==source_id]
                    if not rows.empty:
                        print('\n{}:'.format(cat_name))
                        at.Table.from_pandas(rows).pprint()
                        
    def _catalog_check(self, cat_name):
        """
        Check to see if the name of the ingested catalog is valid
        
        Parameters
        ----------
        cat_name: str
            The name of the catalog in the Catalog object
        
        Returns
        -------
        str
            The catalog name in astroquery format
        """
        good = True
        
        # Make sure the attribute name is good
        if cat_name[0].isdigit():
            print("No names beginning with numbers please!")
            good = False
            
        # Make sure catalog is unique
        if cat_name in self.catalogs:
            print('Catalog {} already ingested.'.format(cat_name))
            good = False
        
        return good
    
    def Vizier_query(self, viz_cat, cat_name, ra, dec, radius, ra_col='RAJ2000', dec_col='DEJ2000', group=True):
        """
        Use astroquery to search a catalog for sources within a search cone
        
        Parameters
        ----------
        viz_cat: str
            The catalog string from Vizier (e.g. 'II/246' for 2MASS PSC)
        cat_name: str
            A name for the imported catalog (e.g. '2MASS')
        ra: astropy.units.quantity.Quantity
            The RA of the center of the cone search
        dec: astropy.units.quantity.Quantity
            The Dec of the center of the cone search
        radius: astropy.units.quantity.Quantity
            The radius of the cone search
        """
        # Verify the cat_name
        if self._catalog_check(cat_name):
            
            # Prep the current catalog as an astropy.QTable
            tab = at.Table.from_pandas(self.catalog)
            
            # Cone search Vizier
            print("Searching {} for sources withiin {} of ({}, {}). Please be patient...".format(viz_cat, radius, ra, dec))
            crds = coord.SkyCoord(ra=ra, dec=dec, frame='icrs')
            data = Vizier.query_region(crds, radius=radius, catalog=viz_cat)[0]
            
            # Ingest the data
            self.ingest_data(data, cat_name, 'id', ra_col=ra_col, dec_col=dec_col)
            
            # Regroup
            if len(self.catalogs)>1 and group:
                self.group_sources(self.xmatch_radius)
            
    def Vizier_xmatch(self, viz_cat, cat_name, ra_col='_RAJ2000', dec_col='_DEJ2000', radius='', group=True):
        """
        Use astroquery to pull in and cross match a catalog with sources in self.catalog
        
        Parameters
        ----------
        viz_cat: str
            The catalog string from Vizier (e.g. 'II/246' for 2MASS PSC)
        cat_name: str
            A name for the imported catalog (e.g. '2MASS')
        radius: astropy.units.quantity.Quantity
            The matching radius
        """
        # Make sure sources have been grouped
        if self.catalog.empty:
            print('Please run group_sources() before cross matching.')
            return
            
        if self._catalog_check(cat_name):
            
            # Verify the cat_name
            viz_cat = "vizier:{}".format(viz_cat)
            
            # Prep the current catalog as an astropy.QTable
            tab = at.Table.from_pandas(self.catalog)
            
            # Crossmatch with Vizier
            print("Cross matching {} sources with {} catalog. Please be patient...".format(len(tab), viz_cat))
            data = XMatch.query(cat1=tab, cat2=viz_cat, max_distance=radius or self.xmatch_radius*q.deg, colRA1='ra', colDec1='dec', colRA2=ra_col, colDec2=dec_col)
            
            # Ingest the data
            self.ingest_data(data, cat_name, 'id', ra_col=ra_col, dec_col=dec_col)
            
            # Regroup
            if group:
                self.group_sources(self.xmatch_radius)
    
    def group_sources(self, radius='', plot=False):
        """
        Calculate the centers of the point clusters given the
        radius and minimum number of points

        Parameters
        ----------
        coords: array-like
            The list of (x,y) coordinates of all clicks
        radius: int
            The distance threshold in degrees for cluster membership
            [default of 0.36 arcseconds]

        Returns
        -------
        np.ndarray
            An array of the cluster centers
        """
        if len(self.catalogs)==0:
            print("No catalogs to start grouping! Add one with the ingest_data() method first.")
            
        else:
            
            # Gather the catalogs
            cats = pd.concat([getattr(self, cat_name) for cat_name in self.catalogs])
            
            # Clear the source grouping
            cats['oncID'] = np.nan
            cats['oncflag'] = ''
            self.xmatch_radius = radius or self.xmatch_radius
            
            # Make a list of the coordinates of each catalog row
            coords = cats[['ra','dec']].values
            
            # Perform DBSCAN to find clusters
            db = DBSCAN(eps=radius, min_samples=1, n_jobs=-1).fit(coords)
            
            # Group the sources
            core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
            core_samples_mask[db.core_sample_indices_] = True
            source_ids = db.labels_+1
            unique_source_ids = list(set(source_ids))
            self.n_sources = len(unique_source_ids)
            
            # Get the average coordinates of all clusters
            unique_coords = np.asarray([np.mean(coords[source_ids==id], axis=0) for id in list(set(source_ids))])
            
            # Generate a source catalog
            self.catalog = pd.DataFrame(columns=('id','ra','dec','flag'))
            self.catalog['id'] = unique_source_ids
            self.catalog[['ra','dec']] = unique_coords
            # self.catalog['flag'] = ['d{}'.format(i) if i>1 else '' for i in Counter(source_ids).values()]
            self.catalog['datasets'] = Counter(source_ids).values()
            
            # Update history
            self.history += "\n{}: Catalog grouped with radius {} arcsec.".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), self.xmatch_radius)
            
            # Update the source_ids in each catalog
            cats['source_id'] = source_ids
            for cat_name in self.catalogs:
                
                # Get the source_ids for the catalog
                cat_source_ids = cats.loc[cats['catID'].str.startswith(cat_name)]['source_id']
                
                # Get the catalog
                cat = getattr(self, cat_name)
                
                # Update the source_ids and put it back
                cat['source_id'] = cat_source_ids
                setattr(self, cat_name, cat)
                
                del cat, cat_source_ids
                
            del cats
            
            # Plot it
            if plot:
                plt.figure()
                plt.title('{} clusters for {} sources'.format(self.n_sources,len(coords)))
                
                colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, n_sources)]
                for k, col in zip(unique_source_ids, colors):
                    
                    class_member_mask = (source_ids == k)
                    xy = coords[class_member_mask & core_samples_mask]
                    
                    marker = 'o'
                    if len(xy)==1:
                        col = [0,0,0,1]
                        marker = '+'
                        
                    plt.plot(xy[:, 0], xy[:, 1], color=tuple(col), marker=marker, markerfacecolor=tuple(col))
                    
    def find_outliers(self):
        """
        Find pairwise distance mean of each cluster and flag outliers
        """
        pass
    
    def drop_catalog(self, cat_name):
        """
        Remove an imported catalog from the Dataset object
        
        Parameters
        ----------
        cat_name: str
            The name given to the catalog
        """
        # Delete the name and data
        self.catalogs.pop(cat_name)
        delattr(self, cat_name)
        
        # Update history
        print("Deleted {} catalog.".format(cat_name))
        self.history += "\n{}: Deleted {} catalog.".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), cat_name)
        
    def load(self, path):
        """
        Load the catalog from file
        
        Parameters
        ----------
        path: str
            The path to the file
        """
        # Get the object
        DB = joblib.load(path)
        
        # Load the attributes
        self.catalog   = DB.catalog
        self.n_sources = DB.n_sources
        self.name      = DB.name
        self.history   = DB.history
        
        del DB
        
    def save(self, path):
        """
        Save the catalog to file for faster loading next time
        
        Parameters
        ----------
        path: str
            The path to the file
        """
        joblib.dump(self, path)
        
    def correct_offsets(self, cat_name, truth='ACS'):
        """
        Function to determine systematic, linear offsets between catalogs
        
        FUTURE -- do this with TweakReg, which also accounts for rotation/scaling
        See thread at https://github.com/spacetelescope/drizzlepac/issues/77
        
        Parameters
        ----------
        cat_name: str
            Name of catalog to correct
        truth: str
            The catalog to measure against
        """
        # Must be grouped!
        if not self.xmatch_radius:
            
            print("Please run group_sources() before running correct_offsets().")
            
        else:
            
            # First, remove any previous catalog correction
            self.catalog.loc[self.catalog['cat_name']==cat_name, 'ra_corr'] = self.catalog.loc[self.catalog['cat_name']==cat_name, '_RAJ2000']
            self.catalog.loc[self.catalog['cat_name']==cat_name, 'dec_corr'] = self.catalog.loc[self.catalog['cat_name']==cat_name, '_DEJ2000']
            
            # Copy the catalog
            onc_gr = self.catalog.copy()
            
            # restrict to one-to-one matches, sort by oncID so that matches are paired
            o2o_new = onc_gr.loc[(onc_gr['oncflag'].str.contains('o')) & (onc_gr['cat_name'] == cat_name) ,:].sort_values('oncID')
            o2o_old = onc_gr.loc[(onc_gr['oncID'].isin(o2o_new['oncID']) & (onc_gr['cat_name'] == truth)), :].sort_values('oncID')
            
            # get coords
            c_o2o_new = SkyCoord(o2o_new.loc[o2o_new['cat_name'] == cat_name, 'ra_corr'],\
                                 o2o_new.loc[o2o_new['cat_name'] == cat_name, 'dec_corr'], unit='degree')
            c_o2o_old = SkyCoord(o2o_old.loc[o2o_old['cat_name'] == truth, 'ra_corr'],\
                                 o2o_old.loc[o2o_old['cat_name'] == truth, 'dec_corr'], unit='degree')
                             
            print(len(c_o2o_old), 'one-to-one matches found!')
            
            if len(c_o2o_old)>0:
                
                delta_ra = []
                delta_dec = []
                
                for i in range(len(c_o2o_old)):
                    # offsets FROM ACS TO new catalog
                    ri, di = c_o2o_old[i].spherical_offsets_to(c_o2o_new[i])
                    
                    delta_ra.append(ri.arcsecond)
                    delta_dec.append(di.arcsecond)
                    
                    progress_meter((i+1)*100./len(c_o2o_old))
                    
                delta_ra = np.array(delta_ra)
                delta_dec = np.array(delta_dec)
                
                print('\n')
                
                # fit a gaussian
                mu_ra, std_ra = norm.fit(delta_ra)
                mu_dec, std_dec = norm.fit(delta_dec)
                
                # Fix precision
                mu_ra = round(mu_ra, 6)
                mu_dec = round(mu_dec, 6)
                
                # Update the coordinates of the appropriate sources
                print('Shifting {} sources by {}" in RA and {}" in Dec...'.format(cat_name,mu_ra,mu_dec))
                self.catalog.loc[self.catalog['cat_name']==cat_name, 'ra_corr'] += mu_ra
                self.catalog.loc[self.catalog['cat_name']==cat_name, 'dec_corr'] += mu_dec
                
                # Update history
                now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.history += "\n{}: {} sources shifted by {} deg in RA and {} deg in Declination.".format(now, cat_name, mu_ra, mu_dec)
                
                # Regroup the sources since many have moved
                self.group_sources(self.xmatch_radius)
                
            else:
                
                print('Cannot correct offsets in {} sources.'.format(cat_name))

def progress_meter(progress):
    """
    Print nice progress update
    
    Parameters
    ----------
    progress: float
        Some fraction of the completed job
    """
    sys.stdout.write("\rloading... %.1f%%" % progress)
    sys.stdout.flush()
