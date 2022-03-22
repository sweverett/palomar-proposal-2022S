import numpy as np
from astropy.table import Table
from emission_lines import emission_line_index
import utils

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
                    bad[catalog[col][indx] < value['min']] = True
                except:
                    pass

                try:
                    bad[catalog[col][indx] > val['max']] = True
                except:
                    pass

                try:
                    bad[catlog[col][indx] != val['equal']] = True
                except:
                    pass

        else:
            try:
                bad[catalog[col] < value['min']] = True
            except:
                pass

            try:
                bad[catalog[col] > val['max']] = True
            except:
                pass

            try:
                bad[catlog[col] != val['equal']] = True
            except:
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

    utils.make_dir(os.dirname(outfile))
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

    utils.make_dir(os.dirname(outfile))
    cut_cat.write(outfile, overwrite=overwrite)

    if return_cat is True:
        return cut_cat
    else:
        return

def apply_post_cuts(matched_file, cuts, outfile, overwrite=False, return_cat=False):

    cat = Table.read(matched_file)

    cut_cat = apply_cuts(cat, cuts)

    utils.make_dir(os.dirname(outfile))
    cut_cat.write(outfile, overwrite=overwrite)

    if return_cat is True:
        return cut_cat
    else:
        return
