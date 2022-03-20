from astropy.table import Table
from argparse import ArgumentParser
import matplotlib.pyplot as plt

parser = ArgumentParser()

parser.add_argument('cluster_file', type=str,
                    help='Cluster filename to have Rlambda added')
parser.add_argument('-outfile', type=str, default=None,
                    help='Filepath for output file w/ added cols')
parser.add_argument('--overwrite', action='store_true', default=False,
                    help='Set to overwrite output files')
parser.add_argument('--plot', action='store_true', default=False,
                    help='Set to plot separation distribution')

def add_rlambda(cat):

    cat['R_LAMBDA'] = (cat['LAMBDA'] / 100.) ** 0.2 # h^-1 Mpc

    return cat

def main(args):
    cluster_file = args.cluster_file
    outfile = args.outfile
    overwrite = args.overwrite
    plot = args.plot

    cat = Table.read(cluster_file)

    cat = add_rlambda(cat)

    if outfile is None:
        outfile = cluster_file
    cat.write(outfile, overwrite=overwrite)

    if plot is True:
        plt.hist(cat['R_LAMBDA'], bins=40, ec='k')
        plt.xlabel('R_LAMBDA (h^-1 Mpc)')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.gcf().set_size_inches(9,6)

        plt.show()

    return

if __name__ == '__main__':
    args = parser.parse_args()

    main(args)
