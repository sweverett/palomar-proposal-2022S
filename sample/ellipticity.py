import numpy as np
import os
import matplotlib.pyplot as plt
import utils

def compute_ellipticity(sources, indx=2, plot=False):
    '''
    indx is band indx (ugriz)
    '''

    AB_EXP = sources['AB_EXP']

    e = (1 - (AB_EXP)**2 ) / (1 + (AB_EXP)**2)

    sources['e'] = e

    if plot is True:
        print('Making ellipticity plots...')

        plt.hist(e, ec='k', bins=50)
        plt.xlabel('Ellipticity (1-(A/B)^2) / (1+(A/B)^2)')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.gcf().set_size_inches(8,5)

        plotfile = os.path.join(utils.get_plot_dir(), 'ellipticity_all.png')
        plt.savefig(plotfile, bbox_inches='tight', dpi=300)
        plt.close()

        plt.hist(e, ec='k', bins=np.linspace(0,1,50))
        plt.xlabel('Ellipticity (1-(A/B)^2) / (1+(A/B)^2)')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.gcf().set_size_inches(8,5)

        plotfile = os.path.join(utils.get_plot_dir(), 'ellipticity_bounds.png')
        plt.savefig(plotfile, bbox_inches='tight', dpi=300)
        plt.close()

    return sources
