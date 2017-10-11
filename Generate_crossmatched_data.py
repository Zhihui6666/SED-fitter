from astropy.io import votable
from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt

# Photometry points wavelength
bandlist = ['GALEX_FUV', 'Swift_UVW2', 'Swift_UVM2', 'GALEX_NUV', 'Swift_UVW1', 'Swift_U', 'SDSS_u', 'SDSS_g',
			'Swift_V', 'SDSS_r', 'SDSS_i', 'SDSS_z', '2MASS_J', '2MASS_H', '2MASS_Ks', 'WISE_W1', 'IRAC_I1',
			'IRAC_I2', 'WISE_W2', 'IRAC_I3', 'IRAC_I4', 'WISE_W3', 'IRS_PB', 'WISE_W4', 'WISE_W4\'', 'IRS_PR',
			'MIPS_M1']
wavelist = [1.5280e-1, 2.0305e-1, 2.2281e-1, 2.2710e-1, 2.5891e-1, 3.5012e-1, 3.5949e-1, 4.6404e-1,
			5.4021e-1, 6.1223e-1, 7.4395e-1, 8.8971e-1, 1.2350, 1.6620, 2.1590, 3.3526, 3.50751,
			4.43658, 4.6028, 5.6281, 7.58916, 11.5608, 15.80, 22.0883, 22.0883, 22.32,
			23.2096]

# mJy = 1e26 * erg/s/angstrom
c = 2.99792458e18 # angstrom/s

# Read Herschel data
hd = votable.parse_single_table('Herschel_data.vot')
# Read photometry data
pht = pd.read_csv('Brown_sample/photometry_minus_a.csv')
# Read cross-matched table
ct = Table.read('Brown_GOALS.ipac', format='ascii.ipac')

namelist = []
datadict = {}
for name in pht['Name']:
	name = name.lstrip(' ')
	name = name.rstrip(' ')
	namelist.append(name)


for targname in ct['Name_x']:
	index = namelist.index(targname)
	photometry = pht[index: index+1]
	obj = targname
	targname = targname.replace(' ', '_')
	print(targname)
	spc = Table.read('Brown_sample/{}_spec.dat'.format(targname), format='ascii')
	akari = spc[spc['col4'] == 2]
	spitzer = spc[spc['col4'] == 3]
	plt.figure(figsize=(7,7))
	plt.step(spitzer['col3']/1e4, spitzer['col2']*1e26*spitzer['col3']**2/c)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('log $\lambda$ / $\mu m$ (in observed frame)')
	plt.ylabel('log flux / mJy')
	plt.title('{}'.format(obj))
	plt.savefig('spectrum/{}.pdf'.format(obj))
