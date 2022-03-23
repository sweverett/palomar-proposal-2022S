import numpy as np
from astropy.table import Table
import astropy.units as u
import matplotlib.pyplot as plt

def compute_line_sb(sources, lines, fiber_diam, amp_col='AMPLITUDE', plot=False):
    '''
    This function computes the surface brightness of each desired emission line
    using the fiber diam used for measurements in the source catalog

    From SDSS website: flux units in 10-17 erg cm^-2 s^-1
    From SDSS website: amplitude units in 10-17 erg cm^-2 s^-1 Å^-1
    '''

    N = len(sources)

    amp_unit = (1e-17) * u.erg / (u.cm**2 * u.s * u.AA)

    # see https://www.sdss.org/dr12/algorithms/fluxcal/
    # sb_zp = 3631 * u.Jy # (don't actually need, astropy will do this)

    fiber_rad = (fiber_diam*u.arcsec) / 2.
    fiber_area = np.pi * fiber_rad**2

    min_sb = None
    for line, wavelength in lines.items():
        indices = emission_line_index[line]

        if isinstance(indices, int):
            indices = [indices]

        flux_density = 0
        for indx in indices:
            flux_density += sources[amp_col][:,indx] * amp_unit

        # surface brightness
        sb = flux_density / fiber_area # flux_unit / arcsec^2

        # convert to AB mag
        z = sources['Z']
        wav = wavelength * (1. + z) * u.nm # (config wavelength is in nm)
        mag = flux_density.to(u.ABmag, u.spectral_density(wav))

        # account for fiber area for sb
        mag = mag + 2.5*np.log10(fiber_area)

        sources[f'{line}_sb_flux'] = sb
        sources[f'{line}_sb_mag'] = mag

        if min_sb is None:
            min_sb = sb
            min_sb_mag = None
        else:
            min_sb = np.min([min_sb, sb], axis=0)
            min_sb_mag = np.min([min_sb_mag, mag], axis=0)

    # keep track of minimum line sb
    sources['min_line_sb_flux'] = min_sb
    sources['min_line_sb_mag'] = min_sb_mag

    if plot is True:
        print('Making sb line plots...')
        for line in lines.keys():
            plt.hist(sources[f'{line}_sb_flux'], ec='k', label=line)
            plt.xlabel('Flux')
            plt.ylabel('Counts')
            plt.xscale('log')
            plt.yscale('log')
            plt.gcf().set_size_inches(8,5)

        plotfile = os.path.join(utils.plot_dir(), 'sb_flux.png')
        plt.savefig(plotfile, bbox_inches='tight', dpi=300)
        plt.close()

        for line in lines.keys():
            plt.hist(sources[f'{line}_sb_mag'], ec='k', label=line)
            plt.xlabel('AB Mag')
            plt.ylabel('Counts')
            plt.xscale('log')
            plt.yscale('log')
            plt.gcf().set_size_inches(8,5)

        plotfile = os.path.join(utils.plot_dir(), 'sb_mag.png')
        plt.savefig(plotfile, bbox_inches='tight', dpi=300)
        plt.close()

        plt.hist(sources[f'min_line_sb_mag'], ec='k')
        plt.xlabel('Min AB Mag')
        plt.ylabel('Counts')
        plt.xscale('log')
        plt.yscale('log')
        plt.gcf().set_size_inches(8,5)

        plotfile = os.path.join(utils.plot_dir(), 'min_sb_mag.png')
        plt.savefig(plotfile, bbox_inches='tight', dpi=300)
        plt.close()

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

    N = len(sources)

    sources['visible_lines'] = np.zeros(N, dtype=np.int)

    zsource = sources['Z']

    for line, wavelength in lines.items():

        # to handle doublets, etc.
        try:
            wavelength = eval(wavelength)
            wavelength = np.mean(wavelength)
        except:
            pass

        shifted = wavelength * (1. + zsource)

        visible = (shifted > blue_lim) & (shifted < red_lim)

        sources[f'{line}_visible'] = visible

        sources['visible_lines'][visible] += 1

    if plot is True:
        plt.hist(sources['visible_lines'])
        plt.xlabel('Number of visible emission lines')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.title(f'Visible lines between {blue_lim} & {red_lim} nm')
        plt.gcf().set_size_inches(8,5)

        outfile = os.path.join(utils.get_plot_dir(), 'visible_lines.png')
        plt.savefig(outfile, bbox_inches='tight', dpi=300)
        plt.close()

    return sources

emission_line_index = {
    'halpha': 24,
    'OII': (3, 4)
}
