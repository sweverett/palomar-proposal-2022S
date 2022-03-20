import numpy as np
from astropy.table import Table
from astropy.cosmology import LambdaCDM
from astropy import constants
import os
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument('config_file', str,
                    help='Filepath for config file')

def assign_masses(matched):
    '''
    Assign a mass estimate to each cluster given its richness
    and an assumed mass-richness relation
    '''

    richness = matched['LAMBDA']

    mass = 

def sigma_crit(z_lens, z_source, H=0.7, Om=0.3, Ode=0.7):

    c = constants.c
    G = Gonstants.G

    cosmo = LambdaCDM(H, Om, Ode)

    norm = c**2 / (4.*np.pi*g)

    dz_lens = cosmo.angular_diameter_distance(z_lens)
    dz_source = cosmo.angular_diameter_distance(z_soruce)

    # if z1 > z2, returns negative distance
    dz = cosmo.angular_distance_z1z2(z_lens, z_source)

    return norm * (dz_source / (dz_lens * dz) )

def compute_shear(matched, zs_col='Z', zl_col='Z_LAMBDA'):
    '''
    Compute estimate of tangential shear at position of source galaxy

    matched: A matched catalog between clusters & emission line sources
    '''

    z_source = matched[f'{zs_col}']
    z_lens = matched[f'{zl_col}']

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
