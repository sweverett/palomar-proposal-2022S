import numpy as np
from astropy.table import Table
from astropy.cosmology import LambdaCDM
from astropy import constants
import os
import galsim
from argparse import ArgumentParser
parser = ArgumentParser()
import pudb

parser.add_argument('config_file', type=str,
                    help='Filepath for config file')

def assign_masses(matched, log_M0=14.338, alpha=1.34, lam0=40.):
    '''
    Assign a mass estimate to each cluster given its richness
    and an assumed mass-richness relation
    Return result is predicted mass of cluster, in solar masses.
    Not sure how best to include the mass-richness parameters here,
    but the form and numbers from Melanie's work:
    https://arxiv.org/abs/1603.06953
    M = M0 * (lam/lam0)^alpha
    with MAP parameter values of:
       -- log10 M0 = 14.338
       -- alpha = 1.34
       -- lam0 = 40 (pivot)
    '''

    richness = matched['LAMBDA_cluster']

    log10_mass = log_M0 + alpha * np.log10(richness / lam0)
    mass = 10**(log10_mass)

    return mass

def sigma_crit(z_lens, z_source, H=0.7, Om=0.27, Ode=0.73):
    c = constants.c
    G = Gonstants.G
    cosmo = LambdaCDM(H, Om, Ode)
    norm = c**2 / (4.*np.pi*g)
    dz_lens = cosmo.angular_diameter_distance(z_lens)
    dz_source = cosmo.angular_diameter_distance(z_soruce)

    # if z1 > z2, returns negative distance
    dz = cosmo.angular_distance_z1z2(z_lens, z_source)

    return norm * (dz_source / (dz_lens * dz) )

def compute_shear(matched, zs_col='Z', zl_col='Z_LAMBDA', z_wedge=0.1,
                  Omega_m=0.27, Omega_lam=0.73, concentration=3):
    '''
    Compute estimate of tangential shear at position of source galaxy
    matched: A matched catalog between clusters & emission line sources


    # We need to assume a halo mass-concentration relation.
    # For this mass range that's in the neighborhood of ~3,
    # so just use 3 for now.
    '''

    z_source = matched[zs_col]
    z_lens = matched[zl_col]

    print('Assigning masses...')
    masses = assign_masses(matched)

    # Use GalSim to predict the tangential shear at the position of each source.
    background = z_source > z_lens + z_wedge
    tan_shears = np.zeros_like(z_source)
    tan_shears[~background] = -1.

    print('Creating NFW halos...')
    halos = []
    i = 0
    for mass, z in zip(masses, z_lens):
        if background[i] == True:
            halos.append(
                galsim.NFWHalo(mass, concentration, z, galsim.PositionD(0.,0.), Omega_m, Omega_lam)
                )
        else:
            halos.append(np.nan)
        i += 1

    pos = matched['separation'] * 3600. # convert to arcsec

    print('Shearing sources...')
    shears = np.zeros_like(z_source)
    for i, p in enumerate(pos):
        if background[i] == True:
            pudb.set_trace()
            shears[i] = halos[i].getShear(
                    galsim.PositionD(p), z_source[i]
                )
        else:
            shears[i] = 0.0

    matched['gtan'] = shears

    return matched

def main(args):
    config_file = args.config

    return 0

if __name__ == '__main__':
    args = parser.pars_args()
    rc = main(args)
    if rc == 0:
        print('shear computed succesfully!')
    else:
        print(f'failed w/ a rc = {rc}')
