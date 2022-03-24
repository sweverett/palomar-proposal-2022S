import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.table import Table
from argparse import ArgumentParser
from ellipticity import compute_ellipticity
import utils
import pudb

parser = ArgumentParser()

parser.add_argument('target_file', type=str,
                    help='filename of target list')
parser.add_argument('--show', action='store_true',
                    help='Set to show plots')

def plot_histograms(targets, outdir, show=False):

    N = len(targets)

    cols = ['good_nights_1.3',
            'visible_lines',
            'halpha_sb_flux',
            'OII_sb_flux',
            'min_line_sb_flux',
            'max_line_sb_mag',
            'FRACDEV',
            'e',
            'THETA_EXP',
            'Z',
            'gtan'
            ]

    # bin_lims = [(0, )]
    bins = [30,
            30,
            np.linspace(-1, 10, 20),
            np.linspace(-1, 10, 20),
            np.linspace(-1, 10, 20),
            30,
            30,
            30,
            30,
            30,
            30
            ]
    log = [False, False, True, True, True, False, False, False, False, False, False]
    band = [False, False, False, False, False, False, True, True, True, False, False]

    for i, col in enumerate(cols):
        x = targets[col]
        if band[i] is True:
            indx = 2
            x = x[:,indx]
        if log[i] is True:
            x = np.log10(x)
            p = f'log10({col})'
        else:
            p = col
        plt.hist(x, ec='k', bins=bins[i])
        plt.xlabel(p)
        plt.ylabel('counts')
        plt.yscale('log')
        plt.title(f'{N} targets passed all cuts')
        plt.gcf().set_size_inches(7,4)

        outname = f'hist_{col}.png'
        outfile = os.path.join(outdir, outname)
        plt.savefig(outfile, bbox_inches='tight', dpi=300)

        if show is True:
            plt.show()

        plt.close()

    return

def main(args):
    target_file = args.target_file
    show = args.show

    targets = Table.read(target_file)

    # this got messed up...
    targets = compute_ellipticity(targets)

    outdir = utils.get_plot_dir()

    plot_histograms(targets, outdir, show=show)

    return

if __name__ == '__main__':
    args = parser.parse_args()

    main(args)
