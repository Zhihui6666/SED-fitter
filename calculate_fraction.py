from __future__ import print_function
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib_version = eval(matplotlib.__version__.split(".")[0])
if matplotlib_version > 1:
    plt.style.use("classic")
import sedfit.SED_Toolkit as sedt
from sedfit import model_functions as sedmf
from gsf_core import configImporter
from sedfit import model_functions
from scipy import interpolate
import glob
import numpy as np
import pandas as pd

Mpc2m = 3.086 * 10**22
L_sun = 3.828 * 10**26 # unit: Watt

def getSedData(targname, redshift, distance, dataDict, sedPck, silent=True):
    pht = sedPck["pht"]
    spc = sedPck["spc"]
    #->Settle into the rest frame
    frame = dataDict.get("frame", "rest") #The coordinate frame of the SED; "rest"
                                          #by default.
    if frame == "obs":
        pht = sedt.SED_to_restframe(pht, redshift)
        spc = sedt.SED_to_restframe(spc, redshift)
        if not silent:
            print("[setSedData]: The input SED is in the observed frame!")
    elif frame == "frame":
        if not silent:
            print("[setSedData]: The input SED is in the rest frame!")
    else:
        if not silent:
            print("[setSedData]: The input SED frame ({0}) is not recognised!".format(frame))
    #->Select bands
    bandList_use = dataDict.get("bandList_use", []) #The list of bands to incorporate;
                                                    #use all the available bands if empty.
    bandList_ignore = dataDict.get("bandList_ignore", []) #The list of bands to be
                                                          #ignored from the bands to use.
    pht = sedt.SED_select_band(pht, bandList_use, bandList_ignore, silent)
    phtwave  = pht[0]
    phtflux  = pht[1]
    phtsigma = pht[2]
    phtband  = pht[3]
    spcwave  = spc[0]
    spcflux  = spc[1]
    spcsigma = spc[2]
    dict = {"phtwave":phtwave, "phtflux": phtflux, "phtsigma": phtsigma,
            "spcwave": spcwave, "spcflux": spcflux, "spcsigma": spcsigma}
    return dict


for filename in glob.glob("configs/*.py"):
    config = configImporter(filename)
    filename = filename.lstrip("configs/").rstrip(".py")
    f = open("result/"+filename+"_bestfit.txt", "a+")
    # f = open(filename + "_bestfit.txt", "a+")
    print(filename)
    sedFile = config.sedFile
    sedPck = sedt.Load_SED(sedFile)
    sedData = getSedData(config.targname, config.redshift, config.distance, config.dataDict, sedPck, False)
    modelDict = config.modelDict
    funcLib = sedmf.funcLib
    waveModel = config.waveModel
    parAddDict_all = {}
    parAddDict_all["DL"] = config.distance
    parAddDict_all["z"] = config.redshift
    parAddDict_all["frame"] = "rest"
    frequency = 2.9979e14 / waveModel
    lib = model_functions.funcLib

    parfit = {}
    for i in f:
        i = i.strip("\n")
        i = i.rstrip(" ")
        # if i.find("total") != -1:
        #     clue = "yes"
        i = i.split(" ")
        if i[0] != "Name" and i[0] != "reduced" and i[0] != "total" and i[0] != "fraction" and i[0] != "":
            parfit[i[0]] = float(i[-1])

    string = ""

    par = {}

    for i in modelDict.keys():
        par[i]={}
        for a in modelDict[i].keys():
            if a != "function":
                if modelDict[i][a]["vary"]:
                    par[i][a] = parfit[a]
                else:
                    par[i][a] = modelDict[i][a]["value"]
    model_a = [0.0 for i in range(1000)]
    model_b = [0.0 for i in range(1000)]
    model_c = [0.0 for i in range(1000)]
    model_d = [0.0 for i in range(1000)]
    for i in par.keys():
        if i == "BC03":
            model_a += lib["BC03"]["function"](wave=waveModel, logMs=par[i]["logMs"], age=par[i]["age"], DL=parAddDict_all["DL"], z=parAddDict_all["z"])
        elif i == "DL07":
            model_b += lib["Extinction"]["function"](wave=waveModel, logtau=par["Extinction"]["logtau"])*lib["DL07"]["function"](wave=waveModel, logumin=par[i]["logumin"], logumax=par[i]["logumax"], qpah=par[i]["qpah"], loggamma=par[i]["loggamma"], logMd=par[i]["logMd"], DL=parAddDict_all["DL"], z=parAddDict_all["z"])
        elif i == "CLUMPY":
            model_c += lib["Extinction"]["function"](wave=waveModel, logtau=par["Extinction"]["logtau"])*lib["CLUMPY_intp"]["function"](wave=waveModel, logL=par[i]["logL"], i=par[i]["i"], tv=par[i]["tv"], q=par[i]["q"], N0=par[i]["N0"], sigma=par[i]["sigma"], Y=par[i]["Y"], DL=parAddDict_all["DL"], z=parAddDict_all["z"])
        elif i == "Hot_Dust":
            model_d += lib["Extinction"]["function"](wave=waveModel, logtau=par["Extinction"]["logtau"])*lib["BlackBody"]["function"](wave=waveModel, logOmega=par[i]["logOmega"], T=par[i]["T"])
        elif i == "Cat3d":
            model_c += lib["Extinction"]["function"](wave=waveModel, logtau=par["Extinction"]["logtau"])*lib["Cat3d"]["function"](wave=waveModel, logL=par[i]["logL"], i=par[i]["i"], theta=par[i]["theta"], N0=par[i]["N0"],a=par[i]["a"], DL=parAddDict_all["DL"], z=parAddDict_all["z"])
    y = model_a + model_b + model_c + model_d
    k = interpolate.interp1d(waveModel, y)

    # Calculate the reduced chi^2 for all the objects

    if filename != 'NGC3690_2' and filename != 'NGC3690_3' and filename != 'NGC3690_4':
    # Not enough data points to calculate reduced chi^2 for this object

        chi_2 = 0.
        if sedData["spcwave"] != []:           # For those who have IRS spectrum
            k_spc = k(sedData["spcwave"])
            for i in range(len(k_spc)):
                chi_2 += (k_spc[i]-sedData["spcflux"][i])**2/sedData["spcsigma"][i]**2
            chi_2_spc = chi_2/(len(k_spc)-len(parfit)+1)
            print ("reduced chi^2 for spectrum: {0:.3}".format(chi_2_spc))
            f.writelines("\nreduced chi^2 for spectrum: {0:.3}\n".format(chi_2_spc))
        # For all the sources
        k_pht = k(sedData["phtwave"])
        for i in range(len(k_pht)):
            chi_2 += (k_pht[i] - sedData["phtflux"][i]) ** 2 / sedData["phtsigma"][i] ** 2
        chi_2_total = float(chi_2) / (len(k_pht) + len(sedData["spcwave"]) - len(parfit) + 1)
        print("reduced chi^2 for whole SED: {0:.3}".format(chi_2_total))
        f.writelines("reduced chi^2 for whole SED: {0:.3}\n".format(chi_2_total))

    # Calculate the contribution of each component

    total = np.trapz(y[waveModel > 8.], frequency[waveModel > 8.])
    f_a = np.trapz(model_a[waveModel > 8.], frequency[waveModel > 8.]) / total
    print("fraction of BC03 = {0:.3e}".format(f_a))
    f.writelines("fraction of BC03 = {0:.3e}\n".format(f_a))
    f_b = np.trapz(model_b[waveModel > 8.], frequency[waveModel > 8.]) / total
    print("fraction of DL07 = {0:.3e}".format(f_b))
    f.writelines("fraction of DL07 = {0:.3e}\n".format(f_b))
    f_c = 0
    if sum(model_c) != 0:
        f_c = np.trapz(model_c[waveModel > 8.], frequency[waveModel > 8.]) / total
        print("fraction of Torus = {0:.3e}".format(f_c))
        f.writelines("fraction of Torus = {0:.3e}\n".format(f_c))
    f_d = 0
    if sum(model_d) != 0:
        f_d = np.trapz(model_d[waveModel > 8.], frequency[waveModel > 8.]) / total
        print("fraction of hot dust = {0:.3e}".format(f_d))
        f.writelines("fraction of hot dust = {0:.3e}\n".format(f_d))

    # Calculate luminosity of AGN and DL07
    L_tot = -total * 10**(-29) * 4 * np.pi * (config.distance * Mpc2m)**2 / L_sun
    print('Total IR Luminosity: {0:.3e}'.format(float(L_tot)))
    f.writelines('Total IR Luminosity: {0:.3e}\n'.format(float(L_tot)))
    L_DL07 = L_tot * f_b
    L_AGN = L_tot * (f_c + f_d)
    print('Luminosity of DL07: {0:.3e}'.format(float(L_DL07)))
    f.writelines('Luminosity of DL07: {0:.3e}\n'.format(float(L_DL07)))
    if f_c != 0:
        print('Luminosity of AGN torus: {0:.3e}'.format(float(L_AGN)))
        f.writelines('Luminosity of AGN torus: {0:.3e}'.format(float(L_AGN)))

    f.close()
