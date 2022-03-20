from argparse import ArgumentParser
from airmass import add_good_nights_col
import pudb

parser = ArgumentParser()

parser.add_argument('cluster_file', type=str,
                    help='Filepath for cluster file to add cols to')
parser.add_argument('source_file', type=str,
                    help='Filepath for source file to add cols to')
parser.add_argument('-start_date', type=str, default='2022-08-01',
                    help='Starting date for search period')
parser.add_argument('-end_date', type=str, default='2023-01-31',
                    help='Ending date for search period')
parser.add_argument('-utc_offset', type=int, default=-7,
                    help='UTC offset for observatory (currently doesnt) ' +\
                    'account for time zones')
parser.add_argument('-min_airmass', type=float, default=1.5,
                    help='Minimum required airmass for a target exposure')
parser.add_argument('-outfile', type=str, default=None,
                    help='Filepath for output file w/ added cols')
parser.add_argument('--overwrite', action='store_true', default=False,
                    help='Set to overwrite output files')

def main(args):

    cluster_file = args.cluster_file
    outfile = args.outfile
    start_date = args.start_date
    end_date = args.end_date
    utc_offset = args.utc_offset
    min_airmass = args.min_airmass
    overwrite = args.overwrite

    print(f'Adding cols to cluster file {cluster_file}...')
    add_good_nights_col(
        cluster_file, start_date, end_date, utc_offset,
        min_airmass=min_airmass, overwrite=overwrite
        )

    return 0

if __name__ == '__main__':
    args = parser.parse_args()
    rc = main(args)

    if rc == 0:
        print('All columns added succesfully!')
    else:
        print(f'Script failed with a rc of {rc}')
