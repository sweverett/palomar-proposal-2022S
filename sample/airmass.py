import numpy as np
import os
from multiprocessing import Pool
import matplotlib.pyplot as plt
from astropy.table import Table
from astroplan import Observer
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from astroplan.plots import plot_airmass
from time import time as stime
from argparse import ArgumentParser

# to make plot_airmass quieter...
import warnings
warnings.filterwarnings('ignore')

parser = ArgumentParser()

parser.add_argument('cluster_file', type=str,
                    help='Filepath for cluster file')
parser.add_argument('start', type=str,
                    help='Start time (year-month-day format, e.g. 2022-08-01')
parser.add_argument('end', type=str,
                    help='End time (year-month-day format, e.g. 2023-01-31')
parser.add_argument('-utc_offset', type=int, default=-7,
                    help='UTC offset for observatory timezone (currently ignores daylight savings)')
parser.add_argument('-outfile', type=str, default=None,
                    help='Filepath for outupt file with added good_nights col')
parser.add_argument('-min_airmass', type=float, default=1.5,
                    help='Minimum airmass value to be considered a good target for the night')
parser.add_argument('--overwrite', action='store_true', default=False,
                    help='Set to overwrite outfile')
parser.add_argument('--plot', action='store_true', default=False,
                    help='Set to make plots')

def get_obs_night_times(start, end, utc_offset):
    '''
    start: str
        Starting date (year-month-day format, e.g. 2022-03-19)
    end: str
        Ending date, same as above
    utc_offset: int
        utc offset for timezone (for now, we don't handle daylight savings)

    This function will return a list of astropy Time objects at midnight on each day
    in between the start and end date (inclusive)
    '''

    assert isinstance(utc_offset, int)

    midnight = f'{-utc_offset}:00:00.000'

    print('start: ', start)
    print('end: ', end)

    start = Time(f'{start} {midnight}')
    end = Time(f'{end} {midnight}')
    Ndays = int((end.jd - start.jd) + 1)

    times = []
    for jd in np.linspace(start.jd, end.jd, Ndays):
        times.append(Time(f'{jd}', format='jd'))

    return times

def compute_airmass(targets, observer, time, Ngrid=25):
    '''
    This is a slightly modified subset of the astroplan.plot_airmass() function:
    https://github.com/astropy/astroplan/blob/aa9bf8607027cf5c3132bd568a353d85197e3cfc/astroplan/plots/time_dependent.py#L51
    '''

    from collections.abc import Sequence

    if hasattr(time, 'utcoffset') and use_local_tz:
            tzoffset = time.utcoffset()
            tzname = time.tzname()
            tzinfo = time.tzinfo
    else:
        tzoffset = 0
        tzname = 'UTC'
        tzinfo = None

    # Populate time window if needed.
    time_ut = Time(time)
    if time_ut.isscalar:
        time_ut = time_ut + np.linspace(-12, 12, Ngrid)*u.hour
    elif len(time_ut) == 1:
        warnings.warn('You used a Time array of length 1.  You probably meant '
                      'to use a scalar. (Or maybe a list with length > 1?).',
                      PlotWarning)
    timetoplot = time_ut + tzoffset

    if not isinstance(targets, Sequence):
        targets = [targets]

    for target in targets:
        # Calculate airmass
        airmass = observer.altaz(time_ut, target).secz
        # Mask out nonsense airmasses
        # masked_airmass = np.ma.array(airmass, mask=airmass < 1)
        airmass[airmass < 1] = np.inf

    return airmass

def add_good_nights_col(catalog_file, start, end, utc_offset,
                        ra_tag='RA', dec_tag='DEC', overwrite=False,
                        min_airmass=1.5, outfile=None, plot=False):
    '''
    This function computes the number of nights

    catalog_file: str
        The catalog of sources to loop over.
    start: str
        Starting date (year-month-day format, e.g. 2022-03-19)
    end: str
        Ending date, same as above
    utc_offset: int
        utc offset for timezone (for now, we don't handle daylight savings)
    '''

    palomar = Observer.at_site('palomar')

    catalog = Table.read(catalog_file)
    Nobjs = len(catalog)

    N = len(catalog)
    good = np.zeros(N)

    midnights = get_obs_night_times(start, end, utc_offset)
    Nnights = len(midnights)

    print(f'Considering {Nnights} nights in the observing season')

    print('Computing sunrise & sunset times...')
    ngrid = 10
    sunsets  = [palomar.sun_set_time(
        midnight, which='nearest', n_grid_points=ngrid
        ) + 30*u.minute for midnight in midnights]
    sunrises  = [palomar.sun_rise_time(
        midnight, which='nearest', n_grid_points=ngrid
        ) + 30*u.minute for midnight in midnights]
    # sunrises = palomar.sun_rise_time(midnights, which='nearest', n_grid_points=ngrid) + 30*u.minute

    night_lengths = [sunrises[i] - sunsets[i] for i in range(Nnights)]

    tstart = stime()
    Tsteps = 12
    for i, obj in enumerate(catalog):
        print(f'starting obj {i} of {N}...')

        coords = SkyCoord(obj[ra_tag], obj[dec_tag], frame='icrs', unit='deg')
        target = FixedTarget(coords)

        # add up observing season dates
        for n in range(Nnights):
            sunset = sunsets[n]
            T = night_lengths[n]
            times = sunset + np.linspace(0, T.value, Tsteps)*u.day

            airmass = compute_airmass(target, palomar, times)
            if (airmass < min_airmass).any():
                good[i] += 1

    T = stime() - tstart
    print(f'total time: {T:.2f}s')
    print(f'per obj: {T/N:.2f}s')

    catalog[f'good_nights_{min_airmass:.1f}'] = good

    if outfile is None:
        outfile = catalog_file.replace('.fits', '_airmass.fits')

    catalog.write(outfile, overwrite=overwrite)

    if plot is True:
        plt.hist(good)
        plt.xlabel(f'Nights with airmass < {min_airmass:.1f}')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.gcf().set_size_inches(7,4)
        plt.show()

    return

def main(args):

    cluster_file = args.cluster_file
    start = args.start
    end = args.end
    utc_offset = args.utc_offset
    min_airmass = args.min_airmass
    outfile = args.outfile
    overwrite = args.overwrite
    plot = args.plot

    add_good_nights_col(
        cluster_file, start, end, utc_offset, min_airmass=min_airmass,
        outfile=outfile, overwrite=overwrite, plot=plot
        )

    return

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
