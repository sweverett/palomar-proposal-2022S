from argparse import ArgumentParser
from airmass import add_good_nights_col
from emission_lines import compute_visible_lines
import pudb

parser = ArgumentParser()

parser.add_argument('config_file', type=str,
                    help='Filepath for run config')
parser.add_argument('-outfile', type=str, default=None,
                    help='Filepath for output file w/ added cols')

def main(args):

    config = args.config_file
    outfile = args.outfile

    cluster_file = config['cluster_file']
    source_file = config['source_file']

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

        print(f'Adding `good_nights` to cluster file {cluster_file}...')

        cluster_outfile = config['cluster_outfile']

        start_date = config['airmass']['start_date']
        end_date = config['airmass']['end_date']
        utc_offset = config['airmass']['utc_offset']
        min_airmass = config['airmass']['min_airmass']
        add_good_nights_col(
            cluster_file, start_date, end_date, utc_offset,
            min_airmass=min_airmass, overwrite=overwrite,
            outfile=outfile, plot=plot
            )

    #-----------------------------------------------------------------
    if 'source_preprocess' in config['run']:
        print('Running source catalog preprocessing...')

        print(f'Adding `visible_lines` to source file {source_file}')
        blue_lim = config['visible']['blue_lim']
        red_lim = config['visible']['red_lim']
        lines = config['visible']['lines']
        compute_visible_lines(
            source_file, lines, blue_lim, red_lim,
            outfile=source_outfile, overwrite=overwrite, plot=plot
            )

        print(f'Adding `line_sb` to source file {source_file}')
        compute_line_sb(
            source_file, lines, blue_lim, red_lim, plot=plot
            )

    #-----------------------------------------------------------------
    if 'match' in config['run']:
        print(f'Matching clusters to sources in {source_file}...')

    return 0

if __name__ == '__main__':
    args = parser.parse_args()
    rc = main(args)

    if rc == 0:
        print('All columns added succesfully!')
    else:
        print(f'Script failed with a rc of {rc}')
