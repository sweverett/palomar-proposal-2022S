import numpy as np
from astropy.table import Table

def apply_cuts(catalog, cuts):
    '''
    cluster_cuts: dict
        A dictionary in the format of `colname`: {`min`: val, `max`:val}.
        Values must either be a float or None
    '''

    N = len(catalog)
    bad = np.zeros(N, dtype=bool)

    for col, val in cuts.items():
        if val['min'] is not None:
            bad[catalog[col] < value['min']] = True

        if val['max'] is not None:
            bad[catalog[col] > val['max']] = True

    cut_catalog = catalog[~bad]

    return cut_catalog

def apply_pre_cluster_cuts(cluster_file, cluster_cuts, outfile,
                           overwrite=False, return_cat=False):
    '''
    cluster_file: str
        Filename of cluster catalog to apply cuts to.
        NOTE: This should be a file that has undergone cluster
        preprocessing
    cluster_cuts: dict
        A dictionary in the format of `colname`: {`min`: val, `max`:val}.
        Values must either be a float or None
    '''

    cat = Table.read(cluster_file)

    cut_cat = apply_cuts(cat, cluster_cuts)

    cut_cat.write(outfile, overwrite=overwrite)

    if return_cat is True:
        return cut_cat
    else:
        return

def apply_pre_source_cuts(source_file, source_cuts, outfile,
                          overwrite=False):
    '''
    source_file: str
        Filename of source catalog to apply cuts to.
        NOTE: This should be a file that has undergone source
        preprocessing
    source_cuts: dict
        A dictionary in the format of `colname`: {`min`: val, `max`:val}.
        Values must either be a float or None
    '''

    cat = Table.read(source_file)

    cut_cat = apply_cuts(cat, source_cuts)

    cut_cat.write(outfile, overwrite=overwrite)

    if return_cat is True:
        return cut_cat
    else:
        return
