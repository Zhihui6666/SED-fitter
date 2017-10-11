import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from glob import glob
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
# from sedfit.models.model_cat3d_G import Cat3d_G
# from sedfit.models.model_cat3d_H import Cat3d_H
# from sedfit.models.model_clumpy import CLUMPY_intp


def get_dl(z):
	cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)
	dl = cosmo.luminosity_distance(z).value
	return dl


def read_bestfit(filename):
	f = open(filename)
	content = f.readlines()[1:-1]
	f.close()
	dict = {}
	for i in content:
		i = i.strip('\n ')
		line = i.split(' ')
		while '' in line:
			line.remove('')
		dict[line[0]] = float(line[-3])
	return dict


# Load the best-fit table of CLUMPY model
t = Table.read('other_catalog/pg_fit.ipac', format='ascii.ipac')
wave = 10**np.linspace(0., 3.0, 1000)

for name in glob('configs_PG_G/*'):
	filename = name[13:-3]
	print(filename)
	# # Load the best-fit pictures from different torus models
	# f1 = 'result_pg_H_lna_1/{}_result.png'.format(filename)
	# fitresult1 = read_bestfit('result_pg_H_lna_1/{}_bestfit.txt'.format(filename))
	# f2 = 'result_pg_H_lna_5/{}_result.png'.format(filename)
	# fitresult2 = read_bestfit('result_pg_H_lna_5/{}_bestfit.txt'.format(filename))
	# f3 = 'result_pg_G_lna_1/{}_result.png'.format(filename)
	# fitresult3 = read_bestfit('result_pg_G_lna_1/{}_bestfit.txt'.format(filename))
	# f4 = 'result_pg_G_lna_5/{}_result.png'.format(filename)
	# fitresult4 = read_bestfit('result_pg_G_lna_5/{}_bestfit.txt'.format(filename))
	# f5 = 'pg_clu_results/{}_result.png'.format(filename)
	# # Plot 5 subplots of different torus models
	# plt.figure(figsize=(15, 10))
	# plt.figtext(0.42, 0.92, filename, size=25)
	# plt.figtext(0.05, 0.75, 'lna = 1', size=20, rotation='vertical')
	# plt.figtext(0.05, 0.3, 'lna = 5', size=20, rotation='vertical')
	#
	# plt.subplot(231)
	# p1 = mpimg.imread(f1)
	# plt.imshow(p1)
	# plt.axis('off') # Remove the axis frame
	#
	# plt.subplot(234)
	# p2 = mpimg.imread(f2)
	# plt.imshow(p2)
	# plt.axis('off')
	#
	# plt.subplot(232)
	# p3 = mpimg.imread(f3)
	# plt.imshow(p3)
	# plt.axis('off')
	#
	# plt.subplot(235)
	# p4 = mpimg.imread(f4)
	# plt.imshow(p4)
	# plt.axis('off')
	#
	# plt.subplot(233)
	# p5 = mpimg.imread(f5)
	# plt.imshow(p5)
	# plt.axis('off')
	#
	# # Plot 5 torus models in one plot
	# target = t[t['Name'] == filename]
	# z = target['z']
	# dl = get_dl(z)
	# logLtorus = target['logL_C']
	# i = target['i_C']
	# tv = target['tv_C']
	# q = target['q_C']
	# N0 = target['N0_C']
	# sigma = target['sigma_C']
	# Y = target['Y_C']
	# torus_H_1 = Cat3d_H(fitresult1['a'], fitresult1['h'], fitresult1['N0'],
	# 					fitresult1['i'], fitresult1['logL'],dl, z, wave)
	# torus_H_5 = Cat3d_H(fitresult2['a'], fitresult2['h'], fitresult2['N0'],
	# 					fitresult2['i'], fitresult2['logL'], dl, z, wave)
	# torus_G_1 = Cat3d_G(fitresult3['a'], fitresult3['theta'], fitresult3['N0'],
	# 					fitresult3['i'], fitresult3['logL'], dl, z, wave)
	# torus_G_5 = Cat3d_G(fitresult4['a'], fitresult4['theta'], fitresult4['N0'],
	# 					fitresult4['i'], fitresult4['logL'], dl, z, wave)
	# clumpy = CLUMPY_intp(logLtorus, i, tv, q, N0, sigma, Y, wave, dl, z)
	# ax = plt.subplot(236)
	# plt.plot(wave, torus_H_1, c='r', label='Cat3d_H_lna_1')
	# plt.plot(wave, torus_H_5, c='y', label='Cat3d_H_lna_5')
	# plt.plot(wave, torus_G_1, c='g', label='Cat3d_G_lna_1')
	# plt.plot(wave, torus_G_5, c='b', label='Cat3d_G_lna_5')
	# plt.plot(wave, clumpy, c='k', label='CLUMPY')
	# plt.legend()
	# plt.xscale('log')
	# plt.yscale('log')
	# upper = max(clumpy) * 2
	# plt.ylim([5, upper])
	# plt.xlim([1, 100])
	# plt.xlabel('Rest Wavelength ($\mu$m)', size=12)
	# plt.ylabel('$f_\\nu$ (mJy)', size=12)
	# ax.yaxis.set_label_coords(0.08, 0.5)
	# plt.savefig('torus_comparision/{}.png'.format(filename), dpi=200)
	# # plt.show()
	# plt.close()
	# break


	fitresult1 = read_bestfit('result_pg_H_lna_1/{}_bestfit.txt'.format(filename))
	fitresult2 = read_bestfit('result_pg_H_lna_5/{}_bestfit.txt'.format(filename))
	fitresult3 = read_bestfit('result_pg_G_lna_1/{}_bestfit.txt'.format(filename))
	fitresult4 = read_bestfit('result_pg_G_lna_5/{}_bestfit.txt'.format(filename))
	target = t[t['Name'] == filename]
	z = target['z']
	dl = get_dl(z)
	logLtorus = target['logL_C'][0]
	i = target['i_C'][0]
	tv = target['tv_C'][0]
	q = target['q_C'][0]
	N0 = target['N0_C'][0]
	sigma = target['sigma_C'][0]
	Y = target['Y_C'][0]
	logMd = target['logMd_C'][0]
	f = open('torus_comparision/{}.txt'.format(filename), 'w')
	f.writelines('Honig torus model (lna = 1) parameters:\n')
	f.writelines('a: {}\nh: {}\nN0: {}\ni: {}\nlogL: {}\nlogMd: {}\nlna: {}\n\n'.format(fitresult1['a'],
		fitresult1['h'], fitresult1['N0'], fitresult1['i'], fitresult1['logL'], fitresult1['logMd'], fitresult1['lna']))
	f.writelines('Honig torus model (lna = 5) parameters:\n')
	f.writelines('a: {}\nh: {}\nN0: {}\ni: {}\nlogL: {}\nlogMd: {}\nlna: {}\n\n'.format(fitresult2['a'],
		fitresult2['h'], fitresult2['N0'],fitresult2['i'], fitresult2['logL'], fitresult2['logMd'], fitresult2['lna']))
	f.writelines('Gonzalez torus model (lna = 1) parameters: \n')
	f.writelines('a: {}\ntheta: {}\nN0: {}\ni: {}\nlogL: {}\nlogMd: {}\nlna: {}\n\n'.format(fitresult3['a'],
		fitresult3['theta'], fitresult3['N0'], fitresult3['i'], fitresult3['logL'], fitresult3['logMd'], fitresult3['lna']))
	f.writelines('Gonzalez torus model (lna = 5) parameters: \n')
	f.writelines('a: {}\ntheta: {}\nN0: {}\ni: {}\nlogL: {}\nlogMd: {}\nlna: {}\n\n'.format(fitresult4['a'],
		fitresult4['theta'], fitresult4['N0'], fitresult4['i'], fitresult4['logL'], fitresult4['logMd'], fitresult4['lna']))
	f.writelines('CLUMPY torus model parameters:\n')
	f.writelines('logLtorus: {}\ni: {}\ntv: {}\nq: {}\nN0: {}\nsigma: {}\nY: {}\nlogMd: {}\n'.format(
		logLtorus, i, tv, q, N0, sigma, Y, logMd))
	f.close()

# plt.show()