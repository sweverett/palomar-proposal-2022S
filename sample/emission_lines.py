import numpy as np
from astropy.table import Table
import astropy.units as u

def compute_line_sb(sources, lines, fiber_diam, flux_col='FLUX' plot=False):
    '''
    This function computes the surface brightness of each desired emission line
    using the fiber diam used for measurements in the source catalog

    From SDSS website: flux units in 10-17 erg cm^-2 s^-1 Å^-1
    '''

    N = len(sources)

    flux_unit = (1e-17) * u.erg / (u.cm**2 * u.s * u.A)

    # see https://www.sdss.org/dr12/algorithms/fluxcal/
    sb_zp = 3631 * u.Jy

    fiber_rad = (fiber_diam*u.arcsec) / 2.
    fiber_area = np.pi * fiber_rad**2

    min_sb = None
    for line, wavelength in lines.items():
        indices = emission_line_index[line]

        sb = 0
        for indx in indices:
            sb += sources[flux_col][:,indx] / fiber_area # flux_unit / arcsec^2

        sources[f'{line_sb_flux}'] = sb
        sources[f'{line_sb_mag}'] = -2.5*np.log10(sb / sb_zp)

        if min_sb is None:
            min_sb = sb
        else:
            min_sb = np.min([min_sb, sb], axis=0)

    # keep track of minimum line sb
    sources['min_line_sb'] = min_sb

    return sources

def compute_visible_lines(sources, lines, blue_lim, red_lim, plot=False):
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

    sources = Table.read(source_file)

    N = len(sources)

    sources['visible_lines'] = np.zeros(N, dtype=np.int)

    zsource = sources['Z']

    for line, wavelength in lines.items():

        # to handle doublets, etc.
        try:
            wavelength = np.mean(wavelength)

        # visible = np.zeros(N, dtype=bool)
        shifted = wavelength * (1. + zsource)

        visible = (shifted > blue_lim) and (shifted < red_limit):

        sources[f'{line}_visible'] = visible

        sources['visible_lines'][visible] += 1

    if plot is True:
        plt.hist(sources['visible_lines'])
        plt.xlabel('Number of visible emission lines')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.title(f'Visible lines between {blue_lim} & {red_lim} nm')
        plt.gcf().set_size_inches(9,6)
        plt.show()

    return sources

emission_line_index = {
    'halpha': 24,
    'OII': (3, 4)
}
