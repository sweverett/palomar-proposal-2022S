from abc import abstractmethod
import numpy as np
from astropy.table import Table

class ClusterCatalog(object):
    '''
    Base class for a cluster catalog
    '''

    def __init__(self, filename, cols=None, ext=1):
        '''
        can pass a list of colnames to load in only a subset
        '''

        self.filename = filename

        self.data = Table(fitsio.read(filename, columns=cols, ext=ext))

        self.

        # depends on the catalog type
        self.ra_col = None
        self.dec_col = None
        self.scale_col = None
        self.z_col = None

        return

    @abstractmethod
    def set_colnames(self):
        pass

class RedmapperCatalog(ClusterCatalog):

    def set_colnames(self):

        self.ra_col = 'RA'
        self.dec_col = 'DEC'
        self.scale_col = 'R_LAMBDA'
        self.z_col = 'Z_LAMBDA'

        return
