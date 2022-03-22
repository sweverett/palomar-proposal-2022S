from argparse import ArgumentParser
from airmass import add_good_nights_col
from emission_lines import compute_visible_lines, compute_line_sb
from match import match_clusters2sources
import cuts
import pudb

parser = ArgumentParser()

parser.add_argument('config_file', type=str,
                    help='Filepath for run config')

def run_cluster_preprocessing(cluster_outfile, config):

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
        outfile=outfile, plot=plot
        )

    return

def run_source_preprocessing(config):
    '''
    Source preprocess funcs return cat instead of writing to a file
    as there are multiple steps
    '''

    source_file = config['source_file']
    source_outfile = config['source_outfile']
    overwrite = config['overwrite']
    plot = config['plot']

    sources = Table.read(source_file)

    print(f'Adding `visible_lines` to source file {source_file}')

    blue_lim = config['visible']['blue_lim']
    red_lim = config['visible']['red_lim']
    lines = config['visible']['lines']
    sources = compute_visible_lines(
        sources, lines, blue_lim, red_lim,
        outfile=source_outfile, overwrite=overwrite, plot=plot
        )

    print(f'Adding `line_sb` to source file {source_file}')

    fiber_diam = config['sb']['fiber_diam']
    sources = compute_line_sb(
        sources, lines, fiber_diam, plot=plot
        )

    sources.write(source_outfile, overwrite=overwrite)

    return

def main(args):

    config = args.config_file

    cluster_file = config['cluster_file']
    source_file = config['source_file']

    # pre-processed filenames
    cluster_outfile = config['cluster_outfile']
    source_outfile = config['source_outfile']

    # after cuts filenames
    cluster_cut_otufile = config['cluster_cut_outfile']
    source_cut_otufile = config['source_cut_outfile']
    match_outfile = confg['match_outfile']

    try:
        plot = config['plot']
    except:
        plot = False

    try:
        overwrite = config['overwrite']
    except:
        overwrite = False

    #-----------------------------------------------------------------
    if 'cluster_preprocess' in config['run']:
        print('Running cluster catalog preprocessing...')
        run_cluster_preprocessing(cluster_outfile, config)

    #-----------------------------------------------------------------
    if 'source_preprocess' in config['run']:
        print('Running source catalog preprocessing...')
        run_source_preprocessing(source_outfile, config)

    #-----------------------------------------------------------------
    if 'cluster_cuts' in config['run']:
        cluster_cuts = config['cluster_cuts']
        print(f'Applying initial cluster cuts on {cluster_outfile}...')
        cuts.apply_pre_cluster_cuts(
            cluster_outfile, cluster_cuts, outfile=cluster_cut_outfile
            )

    #-----------------------------------------------------------------
    if 'source_cuts' in config['run']:
        print(f'Applying initial source cuts on {source_outfile}...')
        source_cuts = config['source_cuts']
        cuts.apply_pre_source_cuts(
            source_outfile, source_cuts, outfile=source_cut_outfile
            )

    #-----------------------------------------------------------------
    if 'match' in config['run']:
        print(f'Matching clusters from {cluster_file} to ' +\
              f'sources in {source_file}...')
        match_radius = config['match']['match_radius']
        matched = match_clusters2sources(
            source_outfile, cluster_outfile, outfile=match_outfile,
            match_radius=match_radius, overwrite=overwrite, plot=plot
            )

    #-----------------------------------------------------------------
    if 'match' in config['run']:
        print(f'Matching clusters from {cluster_file} to ' +\
              f'sources in {source_file}...')

    return 0

if __name__ == '__main__':
    args = parser.parse_args()
    rc = main(args)

    if rc == 0:
        print('All columns added succesfully!')
    else:
        print(f'Script failed with a rc of {rc}')
