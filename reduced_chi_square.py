from __future__ import print_function
import matplotlib
matplotlib.use("Agg")
import sedfit.SED_Toolkit as sedt
from sedfit import model_functions as sedmf
from gsf_core import configImporter
from sedfit import model_functions
from sedfit import sedclass as sedsc
from sedfit.fitter import basicclass as bc
import glob
from sedfit import fit_functions as fit

dataDict = {
    "phtName": "Phot",
    "spcName": "IRS",
    "bandList_use": ['Spitzer_IRAC1', 'Spitzer_IRAC2'],
#    "bandList_use": ['2MASS_J', '2MASS_H', '2MASS_Ks',
#                      'Spitzer_IRAC1', 'Spitzer_IRAC2', 'IRAS_60',
#                      'Herschel_PACS_70', 'IRAS_100', 'Herschel_PACS_160','Herschel_SPIRE_250','Herschel_SPIRE_350','Herschel_SPIRE_500',
#                      'JCMT_SCUBA1_450', 'JCMT_SCUBA1_850','Spitzer_IRAC4','Spitzer_IRAC3','IRAS_12','Spitzer_MIPS_24', 'IRAS_25'],
    "bandList_ignore":[ '2MASS_J', '2MASS_H', '2MASS_Ks','Spitzer_IRAC4','Spitzer_IRAC3','IRAS_12','Spitzer_MIPS_24', 'IRAS_25', 'IRAS_60',
                      'Herschel_PACS_70', 'IRAS_100', 'Herschel_PACS_160','Herschel_SPIRE_250','Herschel_SPIRE_350','Herschel_SPIRE_500',
                      'JCMT_SCUBA1_450', 'JCMT_SCUBA1_850'],
#    "bandList_ignore":[ ],
    "frame": "obs",
    #'FUV', 'NUV', 'U', 'B', 'V', 'R', 'I', 'JCMT_SCUBA1_450'
}

for filename in glob.glob("configs/*.py"):
    config = configImporter(filename)
    filename = filename.lstrip("configs/").rstrip(".py")
    f = open("result/"+filename+"_bestfit.txt", "a+")
    # f = open(filename + "_bestfit.txt", "a+")
    print (filename)
    sedFile = config.sedFile
    sedPck = sedt.Load_SED(sedFile)
    sedData1 = sedsc.setSedData(config.targname, config.redshift, config.distance, dataDict, sedPck, False)
    sedData2 = sedsc.setSedData(config.targname, config.redshift, config.distance, config.dataDict, sedPck, False)
    modelDict = config.modelDict
    funcLib = sedmf.funcLib
    waveModel = config.waveModel
    parAddDict_all = {}
    parAddDict_all["DL"] = config.distance
    parAddDict_all["z"] = config.redshift
    parAddDict_all["frame"] = "rest"
    frequency = 1/waveModel
    lib = model_functions.funcLib

    sedModel = bc.Model_Generator(modelDict, funcLib, waveModel, parAddDict_all)

    parfit = {}
    parlist = []
    for i in f:
        i = i.strip("\n")
        i = i.rstrip(" ")
        # if i.find("total") != -1:
        #     clue = "yes"
        i = i.split(" ")
        if i[0] != "Name" and i[0] != "reduced" and i[0] != "total" and i[0] != "fraction" and i[0] != "":
            parfit[i[0]] = float(i[-1])
            parlist.append(float(i[-1]))
    if len(sedData2.get_List("y")) > 100:
        reduced_chi_sq_spc = (fit.logLFunc(parlist, sedData1, sedModel) * -2) / (len(sedData1.get_List("y")) - len(parlist) + 1)
        reduced_chi_sq_all = (fit.logLFunc(parlist, sedData2, sedModel) * -2) / (len(sedData2.get_List("y")) - len(parlist) + 1)
        print (reduced_chi_sq_spc, reduced_chi_sq_all)
        f.writelines("\nreduced chi^2 from 3 to 40 micron: "+str(reduced_chi_sq_spc)+"\n")
        f.writelines("reduced chi^2 for whole SED: "+str(reduced_chi_sq_all)+"\n")
    f.close()
