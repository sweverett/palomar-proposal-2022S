from argparse import ArgumentParser
from airmass import add_good_nights_col
import pudb

parser = ArgumentParser()

parser.add_argument('config_file', type=str,
                    help='Filepath for run config')
parser.add_argument('-outfile', type=str, default=None,
                    help='Filepath for output file w/ added cols')
parser.add_argument('--overwrite', action='store_true', default=False,
                    help='Set to overwrite output files')

def main(args):

    config = args.config_file
    outfile = args.outfile
    overwrite = args.overwrite

    cluster_file = config['cluster_file']
    source_file = config['source_file']

    #-----------------------------------------------------------------
    if 'airmass' in config['run']:
        print(f'Adding cols to cluster file {cluster_file}...')

        start_date = config['airmass']['start_date']
        end_date = config['airmass']['end_date']
        utc_offset = config['airmass']['utc_offset']
        min_airmass = config['airmass']['min_airmass']
        add_good_nights_col(
            cluster_file, start_date, end_date, utc_offset,
            min_airmass=min_airmass, overwrite=overwrite
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
