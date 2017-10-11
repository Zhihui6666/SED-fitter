import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii
import numpy as np
from glob import glob

L_sun = 3.828e33 # erg/s
SFR_ratio = 4.5e-44

def filt_list(list, colname, low=None, up=None):
	'''

	:param list: np.ndarray
	:param low: low limit
	:param up: upper limit
	:return: list after filter
	'''
	if low != None:
		list = list[list[colname]>low]
	if up != None:
		list = list[list[colname]<=up]
	return list


# f = open('gas_mass.txt')
# content = f.readlines()
# content = content[35:100]
# table = []
# for eachline in content:
# 	eachline = eachline.rstrip('\n')
# 	temp = eachline.split("\t")
# 	if temp[-3] != 'cdots':
# 		table.append([temp[0], temp[1].replace(' ', ''), temp[-3], temp[-4]])
#
#
# # Save the gass mass to .ipac
#
# table = np.array(table)
# result_table = Table(rows=table, names=("Name1", "Name2", "Gas_mass", "Stellar_mass"), dtype=("S", "S", "double", 'double'))
# result_table['Gas_mass'].unit = "$M_\odot$"
# result_table['Stellar_mass'].unit = "$M_\odot$"
# ascii.write(result_table, "Gas_mass.ipac", format="ipac")

# Save dust mass to .ipac

# GOALS = []
# for filename in glob('bestfit/*'):
# 	# init contribution of AGN
# 	temp1 = 0.
# 	temp4 = 0.
# 	f = open(filename, 'r') # Open bestfit.txt
# 	content = f.readlines()
# 	f.close()
# 	filename = filename.lstrip('bestfit/').rstrip('.txt')[:-9]
# 	h = open('data/' + filename + '.csv', 'r') # Open config
# 	h_content = h.readlines()
# 	h.close()
# 	data = open('configs/' + filename + '_2.py', 'r') # Open data
# 	d_content = data.readlines()
# 	data.close()
# 	for i in h_content:
# 		i = i.split(',')
# 		if i[-1] == '2MASS_H\n':
# 			L_h = float(i[1])
# 	for i in content:
# 		i = i.strip("\n")
# 		i = i.rstrip(" ")
# 		i = i.split(" ")
# 		if i[0] == 'logMd':
# 			# Dust mass
# 			temp = float(i[-1])
# 		if i[0] == 'fraction':
# 			if i[2] == 'CLUMPY':
# 				temp1 += float(i[-1])
# 			if i[2] == 'hot':
# 				temp1 += float(i[-1])
# 		if i[0] == 'Total':
# 			temp2 = float(i[-1])
# 		if i[0] == 'Luminosity':
# 			if i[2] == 'DL07:':
# 				temp3 = float(i[-1])
# 			if i[2] == 'AGN':
# 				temp4 += float(i[-1])
#
# 	a = d_content[12].split(' ')
# 	d = float(a[2])
# 	M_s = 1.034 * np.log10(d * d * L_h / 1000) + 8.073
# 	SFR = SFR_ratio * temp3 * L_sun
# 	GOALS.append([filename,temp, temp1, temp2, temp3, temp4, M_s, SFR])
# GOALS = np.array(GOALS)
# result_table = Table(rows=GOALS, names=("Name", "M_d", "fAGN", "L_IR", "L_DL07", "L_torus", "M_s", "SFR"),
# 					 dtype=("S", "double", "double", "double", "double", "double", "double", "double"))
# ascii.write(result_table, "GOALS.ipac", format="ipac", overwrite=True)

# GOALS = Table.read('GOALS.ipac', format='ascii.ipac')
# gas_mass = Table.read('Gas_mass.ipac', format='ascii.ipac')
# hrs = Table.read('hrs_all.ipac', format='ascii.ipac')
# kingfish = Table.read('kingfish_all.ipac', format='ascii.ipac')
# pg = Table.read('pg_fit.ipac', format='ascii.ipac')
#
# pg_sub = pg[pg['logMd_C']-pg['logMstar'] > -4]
# print(len(pg), len(pg_sub))

# Plot comparision between M_* and M_d

# z1 = np.polyfit(GOALS['M_s'], GOALS['M_d']-GOALS['M_s'], 1)  #一次多项式拟合，相当于线性拟合
# p1 = np.poly1d(z1)
# z2 = np.polyfit(pg_sub['logMstar'], pg_sub['logMd_C']-pg_sub['logMstar'], 1)
# p2 = np.poly1d(z2)
# print(z1)
# print(z2)
# min1 = min(GOALS['M_s'])
# min2 = min(pg_sub['logMstar'])
# max1 = max(GOALS['M_s'])
# max2 = max(pg_sub['logMstar'])
#
# plt.figure(figsize=(7, 7))
# plt.scatter(hrs['logMstar'], hrs['logMdust']-hrs['logMstar'], label='HRS', c='g', s=8)
# plt.scatter(kingfish['logMstar'], kingfish['logMdust']-kingfish['logMstar'], label='KINGFISH', c='b', s=8)
# plt.scatter(pg['logMstar'], pg['logMd_C']-pg['logMstar'], label='PG Quasar', c='k', s=8)
# plt.scatter(filt_list(GOALS, colname='fAGN', low=0.8)['M_s'], filt_list(GOALS, colname='fAGN', low=0.8)['M_d']-filt_list(GOALS, colname='fAGN', low=0.8)['M_s'],
# 			label="GOALS with fAGN > 80%", c='r', s=15, marker='x')
# plt.scatter(filt_list(GOALS, colname='fAGN', low=0.6, up=0.8)['M_s'], filt_list(GOALS, colname='fAGN', low=0.6, up=0.8)['M_d']-filt_list(GOALS, colname='fAGN', low=0.6, up=0.8)['M_s'],
# 			label="GOALS with 60% < fAGN <= 80%", c='r', s=15, marker='s')
# plt.scatter(filt_list(GOALS, colname='fAGN', low=0.4, up=0.6)['M_s'], filt_list(GOALS, colname='fAGN', low=0.4, up=0.6)['M_d']-filt_list(GOALS, colname='fAGN', low=0.4, up=0.6)['M_s'],
# 			label="GOALS with 40% < fAGN <= 60%", c='r', s=15, marker='*')
# plt.scatter(filt_list(GOALS, colname='fAGN', low=0.2, up=0.4)['M_s'], filt_list(GOALS, colname='fAGN', low=0.2, up=0.4)['M_d']-filt_list(GOALS, colname='fAGN', low=0.2, up=0.4)['M_s'],
# 			label="GOALS with 20% < fAGN <= 40%", c='r', s=15, marker='+')
# plt.scatter(filt_list(GOALS, colname='fAGN', up=0.2)['M_s'], filt_list(GOALS, colname='fAGN', up=0.2)['M_d']-filt_list(GOALS, colname='fAGN', up=0.2)['M_s'],
# 			label="GOALS with fAGN <= 20%", c='r', s=15, marker='D')
# plt.plot([min1, max1], [min1*z1[0]+z1[1], max1*z1[0]+z1[1]], label='GOALS fitting', c='r')
# plt.plot([min2, max2], [min2*z2[0]+z2[1], max2*z2[0]+z2[1]], label='PG fitting', c='k')
# plt.xlabel('log (M$_*$ / M$_\odot$)')
# plt.ylabel('log (M$_d$ / M$_*$)')
# plt.title('M$_d$ versus M$_*$')
# plt.legend()
# plt.show()

# Comparision of gas to dust ratio

# name1 = []
# for i in gas_mass['Name1']:
# 	name1.append(i)
# name2 = []
# for i in gas_mass['Name2']:
# 	name2.append(i)
# dustname = []
# for i in dust_mass['Name']:
# 	dustname.append(i)
# # Collect gas mass
# gas = []
# dust = []
# for i in range(len(dust_mass['Name'])):
# 	name = dustname[i]
# 	index = None
# 	if name[4:] in name1:
# 		index = name1.index(name[4:])
# 	elif name in name2:
# 		index = name2.index(name)
# 	if index:
# 		gas.append(gas_mass[index]['Gas_mass'])
# 		dust.append([dust_mass[i]['M_s'], dust_mass[i]['M_d']])
# gas = np.array(gas)
# dust = np.array(dust)
# hrs_gas = np.log10(10**hrs['logMHI'] + 10**hrs['logMH2'])
# plt.figure(figsize=(7,7))
# plt.scatter(pg['logMstar'], pg['logGDR'], label='PG quasars', c='k')
# plt.scatter(hrs['logMstar'], hrs_gas - hrs['logMdust'], label='HRS', c='g')
# plt.scatter(dust[:, 0], gas - dust[:, 1], label='GOALS', c='r')
# plt.legend()
# plt.xlabel('log (M$_*$ / M$_\odot$)')
# plt.ylabel('log (M$_g$ / M$_d$)')
# plt.show()

# Plot histogram of fAGN

# print(np.mean(GOALS['fAGN']))
# print(max(GOALS['fAGN']))
# print(np.linspace(0, 1, 11))
# plt.hist(GOALS['fAGN'], bins=np.linspace(0,1,11))
# labels = ['0%', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%']
# plt.xticks(np.linspace(0,1,11), labels)
# plt.xlabel('AGN contribution')
# plt.ylabel('Number')
# plt.show()


# Comparision between L_IR and L_torus

# left, width = 0.1, 0.65
# bottom, height = 0.1, 0.65
# bottom_h = left_h = left + width + 0.02
#
# rect_scatter = [left, bottom, width, height]
# rect_histx = [left, bottom_h, width, 0.2]
# rect_histy = [left_h, bottom, 0.2, height]
#
# plt.figure(figsize=(7,7))
#
# axScatter = plt.axes(rect_scatter)
# axHistx = plt.axes(rect_histx)
# axHisty = plt.axes(rect_histy)
#
# axScatter.scatter(np.log10(GOALS['L_IR']), np.log10(GOALS['L_torus']), label='GOALS', c='r', s=10)
# axScatter.scatter(np.log10(pg['Ltotal_C']/L_sun), np.log10(pg['Ltorus_C']/L_sun), label='PG', c='b', s=10)
# axScatter.set_xlabel('log (L$_{IR}$ / L$_\odot$)', fontsize=15)
# axScatter.set_ylabel('log (L$_{torus}$ / L$_\odot$)', fontsize=15)
# # axScatter.set_xscale('log')
# # axScatter.set_yscale('log')
# axScatter.legend()
#
#
# axHistx.hist([np.reshape(np.log10(GOALS['L_IR']), (-1, 1)), np.reshape(np.log10(pg['Ltotal_C']/L_sun), (-1, 1))],
# 			 color=['r', 'b'], histtype='step', stacked=True, fill=False)
# print(np.reshape(np.log10(GOALS['L_torus']), (-1, 1)))
# print(np.reshape(np.log10(pg['Ltorus_C']/L_sun), (-1, 1)))
# axHisty.hist([np.reshape(np.log10(GOALS['L_torus']), (-1, 1)), np.reshape(np.log10(pg['Ltorus_C']/L_sun), (-1, 1))],
# 			 orientation='horizontal', color=['r', 'b'], histtype='step', stacked=True, fill=False)
#
# # Turn off tick labels on marginals
# plt.setp(axHistx.get_xticklabels(), visible=False)
# plt.setp(axHisty.get_yticklabels(), visible=False)
#
# # Set labels on marginals
# axHisty.set_xlabel('Number')
# axHistx.set_ylabel('Number')
# plt.show()


# Comparision between L_IR and L_DL07

# left, width = 0.1, 0.65
# bottom, height = 0.1, 0.65
# bottom_h = left_h = left + width + 0.02
#
# rect_scatter = [left, bottom, width, height]
# rect_histx = [left, bottom_h, width, 0.2]
# rect_histy = [left_h, bottom, 0.2, height]
#
# plt.figure(figsize=(7,7))
#
# axScatter = plt.axes(rect_scatter)
# axHistx = plt.axes(rect_histx)
# axHisty = plt.axes(rect_histy)
#
# axScatter.scatter(np.log10(GOALS['L_IR']), np.log10(GOALS['L_DL07']), label='GOALS', c='r', s=10)
# axScatter.scatter(np.log10(pg['Ltotal_C']/L_sun), np.log10(pg['Ldl07_C']/L_sun), label='PG', c='b', s=10)
# axScatter.set_xlabel('log (L$_{IR}$ / L$_\odot$)', fontsize=15)
# axScatter.set_ylabel('log (L$_{DL07}$ / L$_\odot$)', fontsize=15)
# # axScatter.set_xscale('log')
# # axScatter.set_yscale('log')
# axScatter.legend()
#
#
# axHistx.hist([np.reshape(np.log10(GOALS['L_IR']), (-1, 1)), np.reshape(np.log10(pg['Ltotal_C']/L_sun), (-1, 1))],
# 			 color=['r', 'b'], histtype='step', stacked=True, fill=False)
# axHisty.hist([np.reshape(np.log10(GOALS['L_DL07']), (-1, 1)), np.reshape(np.log10(pg['Ldl07_C']/L_sun), (-1, 1))],
# 			 orientation='horizontal', color=['r', 'b'], histtype='step', stacked=True, fill=False)
#
# # Turn off tick labels on marginals
# plt.setp(axHistx.get_xticklabels(), visible=False)
# plt.setp(axHisty.get_yticklabels(), visible=False)
#
# # Set labels on marginals
# axHisty.set_xlabel('Number')
# axHistx.set_ylabel('Number')
# plt.show()


# Comparision between M_dust and SFR

# plt.figure(figsize=(7,7))
# plt.scatter(hrs['logMdust'], hrs['logSFR'], label='HRS', c='g', s=10)
# plt.scatter(kingfish['logMdust'], kingfish['logSFR'], label='KINGFISH', c='y', s=10)
# plt.scatter(pg['logMd_C'], np.log10(pg['SFR_IR_C']), label='PG', c='b', s=10)
# plt.scatter(GOALS['M_d'], np.log10(GOALS['SFR']), label='GOALS', c='r', s=10)
# plt.ylabel('log (SFR / M$_\odot$ year$^{-1}$)')
# plt.xlabel('log (M$_d$ / M$_\odot$)')
# plt.legend()
# plt.show()
