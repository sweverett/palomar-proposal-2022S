import numpy as np
import os
from astropy.table import Table
from emission_lines import emission_line_index
import utils
import pudb

def apply_cuts(catalog, cuts, lines=None):
    '''
    cuts: dict
        A dictionary in the format of `colname`: {`min`: val, `max`:val}.
        Values must either be a float or None
    '''

    N = len(catalog)
    bad = np.zeros(N, dtype=bool)

    for col, val in cuts.items():
        if col == 'FIT_WARNING':
            # have to deal w/ this separately as it is an array
            assert lines is not None

            for line in lines:
                indx = emission_line_index[line]
                try:
                    bad[catalog[col][:,indx] <= val['min']] = True
                except KeyError:
                    pass

                try:
                    bad[catalog[col][:,indx] >= val['max']] = True
                except KeyError:
                    pass

                try:
                    bad[catalog[col][:,indx] != val['equal']] = True
                except KeyError:
                    pass

        if col in ['FRACDEV', 'AB_EXP', 'THETA_EXP']:
            indx = 2 # r-band; ugriz
            try:
                # pudb.set_trace()
                bad[catalog[col][:,indx] <= val['min']] = True
            except KeyError:
                pass

            try:
                bad[catalog[col][:,indx] >= val['max']] = True
            except KeyError:
                pass

            try:
                bad[catalog[col][:,indx] != val['equal']] = True
            except KeyError:
                pass
        else:
            try:
                bad[catalog[col] <= val['min']] = True
            except KeyError:
                pass

            try:
                bad[catalog[col] >= val['max']] = True
            except KeyError:
                pass

            try:
                bad[catalog[col] != val['equal']] = True
            except KeyError:
                pass

    cut_catalog = catalog[~bad]

    return cut_catalog

def apply_cluster_pre_cuts(cluster_file, cluster_cuts, outfile, overwrite=False,
                       return_cat=False):
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

    utils.make_dir(os.path.dirname(outfile))
    cut_cat.write(outfile, overwrite=overwrite)

    if return_cat is True:
        return cut_cat
    else:
        return

def apply_source_pre_cuts(source_file, source_cuts, lines, outfile, overwrite=False,
                          return_cat=False):
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

    cut_cat = apply_cuts(cat, source_cuts, lines=lines)

    utils.make_dir(os.path.dirname(outfile))
    cut_cat.write(outfile, overwrite=overwrite)

    if return_cat is True:
        return cut_cat
    else:
        return

def apply_post_cuts(targets, cuts, outfile, overwrite=False, return_cat=False):

    cut_cat = apply_cuts(targets, cuts)

    cut_cat.write(outfile, overwrite=overwrite)

    if return_cat is True:
        return cut_cat
    else:
        return
