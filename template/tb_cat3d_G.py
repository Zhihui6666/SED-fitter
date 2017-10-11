import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
from sedfit.fitter.template import Template
from scipy.interpolate import splrep, splev
from sklearn.neighbors import KDTree
from astropy.table import Table


aList = [-1.75, -1.50, -1.25, -1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50]
thetaList = [30, 45, 60]
N0List = [2.5, 5.0, 7.5, 10.0, 12.5]
tauV = 50
Rout = 450
iList = [0, 15, 30, 45, 60, 75, 90]
length = len(aList) * len(thetaList) * len(N0List) * len(iList)
print(length)

# test the difference of random distribution of clouds

# for i in range(10):
# 	f = Table.read('CAT3D_torus_models_garcia17_isotropic/CAT3D.N2.5_a0.00_theta30_Rout450_tauV50_DISTRIB0{0}.sed'.format(i), format='ascii')
# 	flux = f['col3']/f['col1']
# 	wave = f['col2']
# 	plt.plot(np.log10(wave), np.log10(flux), lw=0.5)
# plt.show()

# Conclusion: random distribution of clouds result in no significant difference!

XList = []
tckList = []
counter = 0
for N0 in N0List:
	if N0 == 10.0:
		newaList = [0.5, 0.00, -0.25, -0.50, -0.75, -1.00, -1.25, -1.50, -1.75]
	elif N0 == 12.5:
		newaList = aList[:-3]
	else:
		newaList = aList
	for theta in thetaList:
		for a in newaList:
			f = Table.read('template/CAT3D_torus_models_garcia17_isotropic/CAT3D.N{0:.1f}_a{1:.2f}_theta{2}_Rout450_tauV50_DISTRIB00.sed'.format(N0, a, theta), format='ascii')
			for i in iList:
				index = i/15 + 3
				flux = f['col{0}'.format(index)]/f['col1']
				tck = splrep(f['col2'], flux)
				tckList.append(tck)
				XList.append([a, theta, N0, i])
				print("[{0}%]".format(100. * (counter + 1) / length))
				counter += 1

wave = f['col2']
kdt = KDTree(XList)
print("Interpolation finishes!")
modelInfo = {
	"a": aList,
	"N0": N0List,
	"theta": thetaList,
	"i": iList,
	"wavelength": wave,
}

parFormat = ["a", "theta", "N0", "i"]
readMe = '''
	This template is from: http://www.sungrazer.org/cat3d.html
	The interpolation is tested well!
	'''
templateDict = {
	"tckList": tckList,
	"kdTree": kdt,
	"parList": XList,
	"modelInfo": modelInfo,
	"parFormat": parFormat,
	"readMe": readMe
}
print("haha")
t = Template(tckList=tckList, kdTree=kdt, parList=XList, modelInfo=modelInfo,
			 parFormat=parFormat, readMe=readMe)
print("haha")
t = Template(**templateDict)
print("haha")

fp = open("template/Cat3d_G.tmplt", "w")
# pickle.dump(t, fp)
pickle.dump(templateDict, fp)
fp.close()

# test of template

fp = open("template/Cat3d_G.tmplt", "r")
tpDict = pickle.load(fp)
fp.close()
t = Template(**tpDict)

counter = 0
for N0 in N0List:
	if N0 == 10.0:
		newaList = aList[:-2]
	elif N0 == 12.5:
		newaList = aList[:-3]
	else:
		newaList = aList
	for theta in thetaList:
		for a in newaList:
			f = Table.read('template/CAT3D_torus_models_garcia17_isotropic/CAT3D.N{0:.1f}_a{1:.2f}_theta{2}_Rout450_tauV50_DISTRIB00.sed'.format(N0, a, theta), format='ascii')
			for i in iList:
				index = i/15 + 3
				flux = f['col{0}'.format(index)]/f['col1']
				tck = splrep(f['col2'], flux)
				tckList.append(tck)
				pars = [a, theta, N0, i]
				flux_intp = t(wave, pars)
				print(np.max(abs(flux - flux_intp) / flux_intp))
				print("[{0}%]".format(100. * (counter + 1) / length))
				counter += 1
	if counter > 100:
		break
print(t.get_parFormat())
print(t.readme())
print(t.get_nearestParameters(pars))