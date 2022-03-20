import numpy as np
from astropy.table import Table

def compute_line_sb(source_file, fiber_area):
    pass

def compute_visible_lines(source_file, lines, blue_lim, red_lim, outfile=None,
                          overwrite=False, plot=False):
    '''
    This function computes how many emission lines are visible
    in the given catalog and filter boundaries and saves to an
    output catalog

    lines: dict
        A dictionary in the format of `line_name`: wavelength
        (in nm)
    blue/red_lim: float
        The filter limits in nm
    '''

    if outfile is None:
        if overwrite is False:
            raise Exception('Must either pass an outfile or set '\
                            'overwrite=True')
        else:
            outfile = source_file

    cat = Table.read(source_file)

    N = len(cat)

    cat['visible_lines'] = np.zeros(N, dtype=np.int)

    zsource = cat['Z']

    for line, wavelength in lines.items():

        # visible = np.zeros(N, dtype=bool)
        shifted = wavelength * (1. + zsource)

        visible = (shifted > blue_lim) and (shifted < red_limit):

        cat[f'{line}_visible'] = visible

        cat['visible_lines'][visible] += 1

    cat.write(outfile, overwrite=overwrite)

    if plot is True:
        plt.hist(cat['visible_lines'])
        plt.xlabel('Number of visible emission lines')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.title(f'Visible lines between {blue_lim} & {red_lim} nm')
        plt.gcf().set_size_inches(9,6)
        plt.show()

    return
