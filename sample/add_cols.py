import os
from argparse import ArgumentParser
from astropy.table import Table
from airmass import add_good_nights_col
from emission_lines import compute_visible_lines, compute_line_sb
from match import match_clusters2sources, match_source_catalogs
from cuts import apply_cluster_pre_cuts, apply_source_pre_cuts, apply_post_cuts
from compute_shear import compute_shear
from ellipticity import compute_ellipticity
import cuts
import utils
import pudb

parser = ArgumentParser()

parser.add_argument('config_file', type=str,
                    help='Filepath for run config')

def run_cluster_preprocessing(config):

    cluster_file = config['cluster_file']
    cluster_outfile = config['cluster_outfile']
    overwrite = config['overwrite']
    plot = config['plot']

    print(f'Adding `good_nights` to cluster file {cluster_file}...')

    start_date = config['airmass']['start_date']
    end_date = config['airmass']['end_date']
    utc_offset = config['airmass']['utc_offset']
    min_airmass = config['airmass']['min_airmass']
    add_good_nights_col(
        cluster_file, start_date, end_date, utc_offset,
        min_airmass=min_airmass, overwrite=overwrite,
        outfile=cluster_outfile, plot=plot
        )

    return

def run_source_preprocessing(config):
    '''
    Source preprocess funcs return cat instead of writing to a file
    as there are multiple steps
    '''

    source_file = config['source_file']
    source_outfile = config['source_outfile']
    sdss_file = config['sdss_file']
    targets_outfile = config['targets_outfile']
    overwrite = config['overwrite']
    plot = config['plot']

    print(f'Matching emission sources to SDSS photometry in {sdss_file}')

    match_outfile = source_file.replace('.fits', 'matched_sdss.fits')
    match_radius = eval(config['sdss']['match_radius'])
    matched = match_source_catalogs(
        sdss_file, source_file, match_radius=match_radius, outfile=match_outfile,
        overwrite=overwrite, plot=plot
        )
    sources = matched.cat

    print(f'Adding `visible_lines` to source file {source_outfile}')

    blue_lim = config['visible']['blue_lim']
    red_lim = config['visible']['red_lim']
    lines = config['visible']['lines']
    sources = compute_visible_lines(
        sources, lines, blue_lim, red_lim, plot=plot
        )

    print(f'Adding `line_sb` to source file {source_outfile}')

    fiber_diam = config['line_sb']['fiber_diam']
    sources = compute_line_sb(
        sources, lines, fiber_diam, plot=plot
        )

    print(f'Adding `e` (ellipticity) to source file {source_outfile}')
    sources = compute_ellipticity(
        sources, plot=plot
        )

    sources.write(source_outfile, overwrite=overwrite)

    return

def main(args):

    config = utils.read_yaml(args.config_file)

    cluster_file = config['cluster_file']
    source_file = config['source_file']

    # pre-processed filenames
    cluster_outfile = config['cluster_outfile']
    source_outfile = config['source_outfile']

    # after cuts filenames
    cluster_cut_outfile = config['cluster_cut_outfile']
    source_cut_outfile = config['source_cut_outfile']
    match_outfile = config['match_outfile']

    targets_outfile = config['targets_outfile']

    try:
        plot = config['plot']
    except:
        plot = False

    try:
        overwrite = config['overwrite']
    except:
        overwrite = False

    # make sure plot dir is present
    plot_dir = utils.get_plot_dir()
    out_dir = utils.get_out_dir()
    print('plot_dir:', plot_dir)
    print('out_dir:', out_dir)
    utils.make_dir(plot_dir)
    utils.make_dir(out_dir)
    assert os.path.exists(plot_dir)
    assert os.path.exists(out_dir)

    #-----------------------------------------------------------------
    if 'cluster_preprocess' in config['run']:
        print('Running cluster catalog preprocessing...')
        run_cluster_preprocessing(config)

    #-----------------------------------------------------------------
    if 'source_preprocess' in config['run']:
        print('Running source catalog preprocessing...')
        run_source_preprocessing(config)

    #-----------------------------------------------------------------
    if 'cluster_cuts' in config['run']:
        cluster_cuts = config['cluster_cuts']
        print(f'Applying initial cluster cuts on {cluster_outfile}...')
        apply_cluster_pre_cuts(
            cluster_outfile, cluster_cuts, outfile=cluster_cut_outfile,
            overwrite=overwrite
            )

    #-----------------------------------------------------------------
    if 'source_cuts' in config['run']:
        print(f'Applying initial source cuts on {source_outfile}...')
        source_cuts = config['source_cuts']
        lines = config['visible']['lines']
        apply_source_pre_cuts(
            source_outfile, source_cuts, lines, outfile=source_cut_outfile,
            overwrite=overwrite
            )

    #-----------------------------------------------------------------
    if 'match' in config['run']:
        print(f'Matching clusters from {cluster_cut_outfile} to ' +\
              f'sources in {source_cut_outfile}...')
        match_radius = eval(config['match']['match_radius'])
        matched = match_clusters2sources(
            source_cut_outfile, cluster_cut_outfile, outfile=match_outfile,
            match_radius=match_radius, overwrite=overwrite, plot=plot
            )
        matches = matched.cat
    else:
        matches = Table.read(match_outfile)

    #-----------------------------------------------------------------
    # if 'shear' in config['run']:
    print(f'Matching clusters from {cluster_file} to ' +\
            f'sources in {source_file}...')
    targets = compute_shear(matches)

    if 'cuts' in config['run']:
        print(f'Applying post-cuts & saving to {targets_outfile}...')

        cat_cuts = config['cuts']
        targets = apply_post_cuts(
            targets, cat_cuts, outfile=targets_outfile, return_cat=True,
            overwrite=overwrite
        )

        Ntargets = len(targets)
        print(f'Success! Final target list has {Ntargets} sources')

        if Ntargets > 0:
            print(':)')
        else:
            print(':(')

    return 0

if __name__ == '__main__':
    args = parser.parse_args()
    rc = main(args)
