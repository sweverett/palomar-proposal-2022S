import numpy as np
from astropy.table import Table, hstack
import esutil.htm as htm
import os
from time import time
from argparse import ArgumentParser
import matplotlib.pyplot as plt

parser = ArgumentParser()

parser.add_argument('cluster_file', type=str,
                    help='Filepath for cluster file to add cols to')
parser.add_argument('source_file', type=str,
                    help='Filepath for source file to add cols to')
parser.add_argument('outfile', type=str, default=None,
                    help='Filepath for output file w/ added cols')
parser.add_argument('--overwrite', action='store_true', default=False,
                    help='Set to overwrite output files')
parser.add_argument('--plot', action='store_true', default=False,
                    help='Set to plot separation distribution')

class MatchedCatalog(object):

    def __init__(self, cat1_file, cat2_file,
                 cat1_ratag='ra', cat1_dectag='dec',
                 cat2_ratag='ra', cat2_dectag='dec',
                 cat1_hdu=1, cat2_hdu=1,
                 match_radius=30/3600., depth=14, table_names=None):
        '''
        match_radius is in deg, same as htm
        '''

        self.cat1_file = cat1_file
        self.cat2_file = cat2_file

        self.cat1_ratag  = cat1_ratag
        self.cat2_ratag  = cat2_ratag
        self.cat1_dectag = cat1_dectag
        self.cat2_dectag = cat2_dectag

        self.match_radius = match_radius
        self.depth = depth

        self.cat2 = None
        self.cat1 = None
        self.cat = None # matched catalog
        self.Nobjs = 0

        self.table_names = None

        self._match()

        return

    def _match(self):
        cat1_cat, cat2_cat = self._load_cats()

        h = htm.HTM(self.depth)

        self.matcher = htm.Matcher(
            depth=self.depth,
            ra=cat1_cat[self.cat1_ratag],
            dec=cat1_cat[self.cat1_dectag]
            )

        id_m, id_t, dist = self.matcher.match(
            ra=cat2_cat[self.cat2_ratag],
            dec=cat2_cat[self.cat2_dectag],
            radius=self.match_radius
            )

        self.cat1 = cat1_cat[id_t]
        self.cat2 = cat2_cat[id_m]
        self.cat2['separation'] = dist
        self.dist = dist

        assert len(self.cat1) == len(self.cat2)

        self.cat = hstack([self.cat1, self.cat2], table_names=self.table_names)

        self.Nobjs = len(self.cat)

        return

    def _load_cats(self):
        cat1_cat = Table.read(self.cat1_file)
        cat2_cat = Table.read(self.cat2_file)

        self.Ncat1 = len(cat1_cat)
        self.Ncat2 = len(cat2_cat)

        return cat1_cat, cat2_cat

    def write(self, outfile, overwrite=False):
        self.cat.write(outfile, overwrite=overwrite)

        return

def match_clusters2sources(source_file, cluster_file, outfile,
                           overwrite=False, plot=False):

    start = time()
    matched = MatchedCatalog(cluster_file, source_file,
    # matched = MatchedCatalog(source_file, cluster_file,
                             cat1_ratag='RA', cat1_dectag='DEC',
                             cat2_ratag='RA', cat2_dectag='DEC',
                             # table_names=['source', 'cluster']
                             table_names=['cluster', 'source']
                             )

    T = time() - start

    print(f'matching took {T:.1f}s for {matched.Ncat1} reference sources')

    matched.write(outfile, overwrite=overwrite)

    return matched

def plot_separations(matched, size=(9,5)):

    sep = matched.cat['separation'] * 60. # arcmin
    N = len(sep)

    plt.hist(sep, bins=30, ec='k')
    plt.xlabel('Angular separation (arcmin)')
    plt.ylabel('Counts')
    plt.yscale('log')
    plt.title(f'{N} matches')

    plt.gcf().set_size_inches(size)

    plt.show()

    return

def main(args):

    cluster_file = args.cluster_file
    source_file = args.source_file
    outfile = args.outfile
    overwrite = args.overwrite
    plot = args.plot

    print(f'Matching clusters to sources...')
    matched = match_clusters2sources(
        source_file, cluster_file, outfile, overwrite=overwrite
        )

    if plot is True:
        plot_separations(matched)

    return 0

if __name__ == '__main__':
    args = parser.parse_args()
    rc = main(args)

    if rc == 0:
        print('Matching completed succesfully!')
    else:
        print(f'Script failed with a rc of {rc}')