import numpy as np
from astropy.table import Table
from astropy.cosmology import LambdaCDM
import matplotlib.pyplot as plt

cluster_file = '/Users/sweveret/repos/palomar-proposal-2022S/data/redmapper_dr8_public_v6.3_catalog.fits'
cat = Table.read(cluster_file)

z = cat['Z_LAMBDA']
r = cat['R_LAMBDA']

cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

print('computing angular diameter distances...')
a = cosmo.angular_diameter_distance(z).value
print('done!')
# print('a:',a)

print('computing angular separations...')
theta = (r / a) * (180. / np.pi) * 60. # arcmin
print('done!')
# print('theta:',theta)

print('plotting...')
plt.hist(theta, bins=50, ec='k')
plt.xlabel('Angular Extent (R_Lambda / ang_diam_dist); arcmin)')
plt.ylabel('Counts')
plt.title('SDSS DR8 redMaPPer clusters')
plt.yscale('log')
plt.gcf().set_size_inches(9,6)
plt.show()
