from glob import glob
import numpy as np
import sedfit.SED_Toolkit as sedt
from sedfit.fitter import basicclass as bc
from sedfit import sedclass as sedsc
from astropy.table import Table
import matplotlib.pyplot as plt



IRAC4 = []
MIPS24 = []
IRAC4_e = []
MIPS24_e = []
IRS_IRAC4  = []
IRS_MIPS24  = []
Namelist = []
for filename in glob('configs/*_3.py'):
	config = open(filename)
	content = config.readlines()
	config.close()
	targname = filename.lstrip('configs/')[:-5]
	print(targname)
	z = float(content[11].split(' ')[2])
	data = 'data/' + targname + '.csv'
	sedPck = sedt.Load_SED(data)
	sed = sedPck["pht"]
	spc = sedPck["spc"]
	bandList = ["Spitzer_IRAC4", "Spitzer_MIPS_24"]
	sedName = "Phot"
	spcName = "IRS"
	sed = sedt.SED_to_restframe(sed, z)
	sed = sedt.SED_select_band(sed, bandList)
	spc = sedt.SED_to_restframe(spc, z)
	sedwave = sed[0]
	sedflux = sed[1]
	sedsigma = sed[2]
	sedband = sed[3]
	spcwave = np.array(spc[0])
	spcflux = np.array(spc[1])
	spcsigma = np.array(spc[2])

	if spcwave != [] and 7.872 > spcwave[0] and 24 < spcwave[-1] and len(sedflux) == 2:
		sedflag = np.ones_like(sedwave)
		sedDataType = ["name", "wavelength", "flux", "error", "flag"]
		phtData = {sedName: bc.DiscreteSet(sedband, sedwave, sedflux, sedsigma, sedflag, sedDataType)}
		spcflag = np.ones_like(spcwave)
		spcDataType = ["wavelength", "flux", "error", "flag"]
		spcData = {"IRS": bc.ContinueSet(spcwave, spcflux, spcsigma, spcflag, spcDataType)}
		sedData = sedsc.SedClass(targname, z, phtDict=phtData, spcDict=spcData)
		sedData.set_bandpass(sedband, sedwave)

		irac4 = sedflux[0]
		mips24 = sedflux[1]
		irac4_e = sedsigma[0]
		mips24_e = sedsigma[1]

		# Calculate the W3 and W4 from IRS spectra
		sp_irac4 = sedData.filtering("Spitzer_IRAC4", spcwave, spcflux)[1]
		sp_mips24 = sedData.filtering("Spitzer_MIPS_24", spcwave, spcflux)[1]
		# Collect the data
		IRAC4.append(irac4)
		MIPS24.append(mips24)
		IRAC4_e.append(irac4_e)
		MIPS24_e.append(mips24_e)
		IRS_IRAC4.append(sp_irac4)
		IRS_MIPS24.append(sp_mips24)
		Namelist.append(targname)

index = IRAC4.index(min(IRAC4))
print(Namelist[index])
# Compare the flux between observation and spectrum
plt.figure(figsize=(7, 7))
plt.errorbar(x=IRAC4, y=IRS_IRAC4, xerr=IRAC4_e, linestyle='none', marker='o')
plt.plot((min(IRS_IRAC4), max(IRS_IRAC4)), (min(IRS_IRAC4), max(IRS_IRAC4)), linestyle='--', c='b')
plt.xscale('log')
plt.yscale('log')
plt.title('IRAC4')
plt.xlabel('IRAC4 flux (mJy)')
plt.ylabel('IRS IRAC4 (mJy)')
plt.savefig('IRAC4.png', dpi=200)

plt.figure(figsize=(7, 7))
plt.errorbar(MIPS24, IRS_MIPS24, xerr=MIPS24_e, linestyle='none', marker='o')
plt.plot((min(IRS_MIPS24), max(IRS_MIPS24)), (min(IRS_MIPS24), max(IRS_MIPS24)), linestyle='--', c='b')
plt.xscale('log')
plt.yscale('log')
plt.title('MIPS24')
plt.xlabel('MIPS24 flux (mJy)')
plt.ylabel('IRS MIPS24 (mJy)')
plt.savefig('MIPS24.png', dpi=200)
