from __future__ import print_function
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib_version = eval(matplotlib.__version__.split(".")[0])
if matplotlib_version > 1:
    plt.style.use("classic")
import sys
import types
import numpy as np
import importlib
from time import time
import cPickle as pickle
import sedfit.SED_Toolkit as sedt
from sedfit.fitter import basicclass as bc
from sedfit.mcmc import mcmc_emcee as mcmc
from sedfit import sedclass as sedsc
from sedfit import model_functions as sedmf
from matplotlib.ticker import FuncFormatter, FormatStrFormatter

__all__ = ["configImporter", "fitter", "gsf_fitter"]

def mjrFormatter(x, pos):
    """
    Define the function to setup the major axis tick label.
    """
    return "$10^{{{0:.0f}}}$".format(np.log10(x))


def ticksFinder(ymin, ymax,
                yTicksTry=np.array([0, 1e0, 1e1, 1e2, 1e3, 1e4])):
    """
    Find the proper ticklabel from ymin to ymax in logscale.
    """
    yTicksLabels = yTicksTry[(yTicksTry>ymin) & (yTicksTry<ymax)]
    if len(yTicksLabels) > 1:
        midTick = (np.log10(ymax)+np.log10(ymin))/2.0
        fltr_label = np.argmin(np.abs(np.log10(yTicksLabels) - midTick))
        yTicksLabel = yTicksLabels[fltr_label]
    else:
        yTicksLabel = yTicksLabels[0]
    return yTicksLabel

def configImporter(configfile):
    """
    This function import the provided configure file.

    Parameters
    ----------
    configfile : string
        The name of the configure file (with the path).

    Returns
    -------
    config : module object
        The imported module.

    Notes
    -----
    None.
    """
    pathList = configfile.split("/")
    configPath = "/".join(pathList[0:-1])
    sys.path.append(configPath)
    configName = pathList[-1].split(".")[0]
    config = importlib.import_module(configName)
    return config

def fitter(sedData, sedModel, unctDict, parTruth, emceeDict, mpi_pool=None):
    """
    This function is run the SED fitting with the MCMC method.

    Parameters
    ----------
    sedData : SEDClass object
        The data set of SED.
    sedModel : ModelCombiner object
        The combined model. The parameters are set to generate the mock SED.
    unctDict : dict
        {
            "lnf" : float, (-inf, lnf_max]
                The ln of f, the imperfectness of the model.
            "lna" : float, (-inf, lnf_max]
                The ln of a, the amplitude of the residual correlation.
            "lntau" : float, (-inf, lnf_max]
                The ln of tau, the scale length of the residual correlation.
        }
    parTruth : bool
        The toggle whether to provide the truth of the model.
    emceeDict : dict
        The dict containing the parameters for emcee to sample the parameter space.
    mpi_pool : (optional) emcee.mpi_pool.MPIPool object
        The pool of MPI to run, if provided.

    Returns
    -------
    em : EmceeModel object
        The object of EmceeModel.

    Notes
    -----
    None.
    """
    #->Prepare to run the iteration
    t0 = time()
    setupKeys = emceeDict["Setup"].keys()
    print( "\n#{:-^50}#".format("emcee Setups") )
    if not mpi_pool is None:
        setupKeys.remove("threads")
        print("**MPI mode")
    for keys in setupKeys:
        print("{0}: {1}".format(keys, emceeDict["Setup"][keys]))
    threads   = emceeDict["Setup"]["threads"]
    printFrac = emceeDict["Setup"]["printfrac"]
    psLow     = emceeDict["Setup"]["pslow"]
    psCenter  = emceeDict["Setup"]["pscenter"]
    psHigh    = emceeDict["Setup"]["pshigh"]
    #->Start the iteration
    runList = emceeDict.keys()
    runList.remove("Setup")
    for loop_run in range(len(runList)):
        runName = runList[loop_run]
        #->Print the fitting stage.
        runDict = emceeDict[runName]
        runKeys = runDict.keys()
        SamplerType = runDict.get("sampler", "EnsembleSampler")
        nwalkers    = runDict.get("nwalkers", 100)
        iteration   = runDict.get("iteration", [500, 500])
        thin        = runDict.get("thin", 1)
        ballR       = runDict.get("ball-r", 0.1)
        print( "\n#{:-^50}#".format( " {0} ".format(runName) ) )
        if (SamplerType == "EnsembleSampler") & ("ntemps" in runKeys):
            runKeys.remove("ntemps")
        for keys in runKeys:
            print("{0}: {1}".format(keys, runDict[keys]))
        #->Setup the sampler
        if unctDict is None:
            modelUnct = False
        else:
            modelUnct = True
        em = mcmc.EmceeModel(sedData, sedModel, modelUnct, unctDict, SamplerType)
        if SamplerType == "EnsembleSampler":
            if mpi_pool is None:
                sampler = em.EnsembleSampler(nwalkers, threads=threads)
            else:
                sampler = em.EnsembleSampler(nwalkers, pool=mpi_pool)
            if loop_run == 0: #If it is the first iteration, the initial position of the walkers are set.
                p0 = [em.from_prior() for i in range(nwalkers)]
            else:
                p0 = em.p_ball(pmax, ratio=ballR)
        elif SamplerType == "PTSampler":
            ntemps = runDict["ntemps"]
            if mpi_pool is None:
                sampler = em.PTSampler(ntemps, nwalkers, threads=threads)
            else:
                sampler = em.PTSampler(ntemps, nwalkers, pool=mpi_pool)
            if loop_run == 0:#If it is the first iteration, the initial position of the walkers are set.
                p0 = []
                for i in range(ntemps):
                    p0.append([em.from_prior() for i in range(nwalkers)])
            else:
                p0 = em.p_ball(pmax, ratio=ballR)
        #->Run the MCMC sampling
        for i in range(len(iteration)):
            em.reset()
            steps = iteration[i]
            print( "\n{:*^35}".format(" {0}th {1} ".format(i, runName)) )
            em.run_mcmc(p0, iterations=steps, printFrac=printFrac, thin=thin)
            em.diagnose()
            pmax = em.p_logl_max()
            em.print_parameters(truths=parTruth, burnin=0)
            em.plot_lnlike(filename="gsf_temp_lnprob.png", histtype="step")
            print( "**Time ellapse: {0:.3f} hour".format( (time() - t0)/3600. ) )
            p0 = em.p_ball(pmax, ratio=ballR)
    return em


def gsf_fitter(configName, targname=None, redshift=None, distance=None, sedFile=None, mpi_pool=None):
    """
    The wrapper of fitter() function. If the targname, redshift and sedFile are
    provided as arguments, they will be used overriding the values in the config
    file saved in configName. If they are not provided, then, the values in the
    config file will be used.

    Parameters
    ----------
    configName : str
        The full path of the config file.
    targname : str or None by default
        The name of the target.
    redshift : float or None by default
        The redshift of the target.
    distance : float or None by default
        The distance of the source from the Sun.
    sedFile : str or None by default
        The full path of the sed data file.
    mpi_pool : (optional) emcee.mpi_pool.MPIPool object
        The pool of MPI to run, if provided.

    Returns
    -------
    None.

    Notes
    -----
    None.
    """
    ############################################################################
    #                                Setup                                     #
    ############################################################################
    config = configImporter(configName)
    if targname is None:
        assert redshift is None
        assert distance is None
        assert sedFile is None
        targname = config.targname
        redshift = config.redshift
        distance = config.distance
        sedFile  = config.sedFile
    else:
        assert not redshift is None
        assert not sedFile is None
    print("#--------------------------------#")
    print("Target:      {0}".format(targname))
    print("Redshift:    {0}".format(redshift))
    print("Distance:    {0}".format(distance))
    print("SED file:    {0}".format(sedFile))
    print("Config file: {0}".format(configName))
    print("#--------------------------------#")

    try:
        silent = config.silent
    except:
        silent = False

    #->Setup the data Data
    dataDict = config.dataDict
    sedPck = sedt.Load_SED(sedFile)
    sedData = sedsc.setSedData(targname, redshift, distance, dataDict, sedPck, silent)

    #->Setup the model
    modelDict = config.modelDict
    print("The model info:")
    parCounter = 0
    for modelName in modelDict.keys():
        print("[{0}]".format(modelName))
        model = modelDict[modelName]
        for parName in model.keys():
            param = model[parName]
            if not isinstance(param, types.DictType):
                continue
            elif param["vary"]:
                print("-- {0}, {1}".format(parName, param["type"]))
                parCounter += 1
            else:
                pass
    print("Varying parameter number: {0}".format(parCounter))
    print("#--------------------------------#")
    funcLib   = sedmf.funcLib
    waveModel = config.waveModel
    try:
        parAddDict_all = config.parAddDict_all
    except:
        parAddDict_all = {}
    parAddDict_all["DL"]    = sedData.dl
    parAddDict_all["z"]     = redshift
    parAddDict_all["frame"] = "rest"
    sedModel  = bc.Model_Generator(modelDict, funcLib, waveModel, parAddDict_all)

    ############################################################################
    #                                   Fit                                    #
    ############################################################################
    modelUnct = config.modelUnct #Whether to consider the model uncertainty in the fitting
    parTruth  = config.parTruth  #Whether to provide the truth of the model
    unctDict = config.unctDict
    emceeDict = config.emceeDict
    em = fitter(sedData, sedModel, unctDict, parTruth, emceeDict, mpi_pool)

    ############################################################################
    #                              Post process                                #
    ############################################################################
    try:
        ppDict = config.ppDict
    except:
        print("[gsf] Warning: cannot find ppDict in the configure file!")
        ppDict = {}
    psLow    = ppDict.get("low", 16)
    psCenter = ppDict.get("center", 50)
    psHigh   = ppDict.get("high", 84)
    nuisance = ppDict.get("nuisance", True)
    fraction = ppDict.get("fraction", 0)
    burnIn   = ppDict.get("burn-in", 50)

    dataPck = {
        "targname": targname,
        "redshift": redshift,
        "distance": sedData.dl,
        "sedPck": sedPck,
        "dataDict": dataDict
    }
    modelPck = {
        "modelDict": modelDict,
        "waveModel": waveModel,
        "parAddDict_all": parAddDict_all,
        "parTruth": parTruth,
        "modelUnct": modelUnct
    }
    fitrs = {
        "dataPck": dataPck,
        "modelPck": modelPck,
        "ppDict": ppDict,
        "posterior_sample": em.posterior_sample(burnin=burnIn, fraction=fraction),
        "chain": em.sampler.chain,
        "lnprobability": em.sampler.lnprobability
    }
    fp = open("result/{0}.fitrs".format(targname), "w")
    pickle.dump(fitrs, fp)
    fp.close()
    #->Save the best-fit parameters
    em.Save_BestFit("result/{0}_bestfit.txt".format(targname), low=psLow, center=psCenter, high=psHigh,
                    burnin=burnIn, fraction=fraction)
    #->Plot the chain of the final run
    em.plot_chain(filename="result/{0}_chain.png".format(targname), truths=parTruth)
    #->Plot the SED fitting result figure
    # Plot the SED data and fit
    ps = fitrs['posterior_sample']
    sedwave = sedData.get_List("x")
    sedflux = sedData.get_List("y")
    spcwave = sedData.get_csList("x")
    spcflux = sedData.get_csList("y")
    if sedData.check_csData():
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)  # The big subplot
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        # ->Plot the upper panel
        xmin = np.nanmin(spcwave) * 0.9
        xmax = np.nanmax(spcwave) * 1.1
        ymin = np.nanmin(spcflux) * 0.8
        ymax = np.nanmax(spcflux) * 1.2
        xlim = [xmin, xmax]
        ylim = [ymin, ymax]
        em.plot_fit(truths=parTruth, FigAx=(fig, ax1), xlim=xlim, ylim=ylim, nSamples=100,
                    burnin=burnIn, fraction=fraction, ps=ps, showLegend=False)
        # plot uncertainty of spectrum
        # spcsig = sedPck['spc'][2]
        # spcflux = np.array(spcflux)
        # spcsig = np.array(spcsig)
        # ax1.fill_between(spcwave, spcflux-spcsig, spcflux+spcsig, facecolors='grey', alpha=0.2)
        # -->Set the labels
        xTickLabels = [10., 20.]
        ax1.set_xticks(xTickLabels)
        ax1.set_xticklabels(xTickLabels)
        ax1.xaxis.set_major_formatter(FormatStrFormatter("%d"))
        yTL = ticksFinder(ymin, ymax, yTicksTry=np.linspace(ymin, ymax, 20))
        yTickLabels = [np.around(yTL, decimals=-1 * int(np.log10(yTL)))]
        ax1.set_yticks(yTickLabels)
        ax1.set_yticklabels(yTickLabels)
        ax1.yaxis.set_major_formatter(FormatStrFormatter("%d"))
        ax1.set_xlabel("")
        ax1.set_ylabel("")
        ax1.tick_params(axis="both", which="major", length=8, labelsize=18)
        ax1.tick_params(axis="both", which="minor", length=5)
        ax1.text(0.05, 0.9, targname,
                 verticalalignment='bottom', horizontalalignment='left',
                 transform=ax.transAxes, fontsize=24,
                 bbox=dict(facecolor='white', alpha=0.5, edgecolor="none"))
        # -->Set the legend
        phtName = dataDict["phtName"]
        spcName = dataDict["spcName"]
        handles, labels = ax1.get_legend_handles_labels()
        handleUse = []
        labelUse = []
        for loop in range(len(labels)):
            lb = labels[loop]
            hd = handles[loop]
            if lb == "Hot_Dust":
                lb = "BB"
            # if lb == "CLUMPY":
            #    lb = "CLU"
            if lb == phtName:
                hd = hd[0]
            if lb != spcName:
                labelUse.append(lb)
                handleUse.append(hd)
            else:
                label_spc = lb
                handle_spc = hd
        labelUse.append(label_spc)
        handleUse.append(handle_spc)
        ax1.legend(handleUse, labelUse, loc="lower right", ncol=2,
                   framealpha=0.9, edgecolor="white",  # frameon=False, #
                   fontsize=16, labelspacing=0.3, columnspacing=0.5,
                   handletextpad=0.3, numpoints=1, handlelength=(4. / 3.))
        # plotName = r"PG {0}${1}${2}".format(targname[2:6], targname[6], targname[7:])
        # ax1.text(0.05, 0.8, "{0}".format(plotName),
        #          verticalalignment='bottom', horizontalalignment='left',
        #          transform=ax1.transAxes, fontsize=24,
        #          bbox=dict(facecolor='white', alpha=0.5, edgecolor="none"))
        # ->Plot the lower panel
        xmin = np.min(sedwave) * 0.9
        xmax = np.max(sedwave) * 1.1
        ymin = np.min(sedflux) * 0.5
        ymax = np.max(sedflux) * 2.0
        xlim = [xmin, xmax]
        ylim = [ymin, ymax]
        em.plot_fit(truths=parTruth, FigAx=(fig, ax2), xlim=xlim, ylim=ylim, nSamples=100,
                    burnin=burnIn, fraction=fraction, ps=ps, showLegend=False)
        # plot uncertainty of spectrum
        # ax2.fill_between(spcwave, spcflux - spcsig, spcflux + spcsig, facecolors='grey', alpha=0.2)
        ax2.set_xlabel("")
        ax2.set_ylabel("")
        ax2.tick_params(axis="both", which="major", length=8, labelsize=18)
        ax2.tick_params(axis="both", which="minor", length=5)
        # -->Set the labels
        yTicksLabels = [ticksFinder(ymin, ymax)]  # [1e1, 1e2, 1e3] #
        ax2.set_yticks(yTicksLabels)
        ax2.set_yticklabels(yTicksLabels)
        ax2.yaxis.set_major_formatter(FuncFormatter(mjrFormatter))
        plt.tight_layout(pad=1.8)
        # ->Setup the shared axis label.
        ax.set_xlabel(r"Rest Wavelength ($\mu$m)", fontsize=24)
        ax.set_ylabel(r"$f_\nu \mathrm{(mJy)}$", fontsize=24)
        ax.xaxis.set_label_coords(0.5, -0.05)
        ax.yaxis.set_label_coords(-0.06, 0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both',  # changes apply to the x-axis
                       which='both',  # both major and minor ticks are affected
                       bottom='off',  # ticks along the bottom edge are off
                       top='off',  # ticks along the top edge are off
                       labelbottom='off',  # labels along the bottom edge are off)
                       labelleft="off")
    else:
        fig = plt.figure(figsize=(7, 7))
        ax = plt.gca()
        xmin = np.min(sedwave) * 0.9
        xmax = np.max(sedwave) * 1.1
        ymin = np.min(sedflux) * 0.001
        ymax = np.max(sedflux) * 2.0
        xlim = [xmin, xmax]
        ylim = [ymin, ymax]
        em.plot_fit(truths=parTruth, FigAx=(fig, ax), xlim=xlim, ylim=ylim, nSamples=100,
                    burnin=burnIn, fraction=fraction, ps=ps)
        ax.text(0.05, 0.9, targname,
                verticalalignment='bottom', horizontalalignment='left',
                transform=ax.transAxes, fontsize=24,
                bbox=dict(facecolor='white', alpha=0.5, edgecolor="none"))
        # plotName = r"PG {0}${1}${2}".format(targname[2:6], targname[6], targname[7:])
        # ax.text(0.58, 0.88, "{0}".format(plotName),
        #         verticalalignment='bottom', horizontalalignment='left',
        #         transform=ax.transAxes, fontsize=24,
        #         bbox=dict(facecolor='white', alpha=0.5, edgecolor="none"))

        # -->Set the legend
        phtName = dataDict["phtName"]
        spcName = dataDict["spcName"]
        handles, labels = ax.get_legend_handles_labels()
        handleUse = []
        labelUse = []
        for loop in range(len(labels)):
            lb = labels[loop]
            hd = handles[loop]
            if lb == "Hot_Dust":
                lb = "BB"
            if lb == "CLUMPY":
                lb = "CLU"
            if lb == "Cat3d":
                lb = "Cat3d"
            if lb == phtName:
                hd = hd[0]
            labelUse.append(lb)
            handleUse.append(hd)
        plt.legend(handleUse, labelUse, loc="lower left", fontsize=18, numpoints=1,
                   handletextpad=0.3, handlelength=(4. / 3.))
    plt.savefig("result/{0}_result.pdf".format(targname), bbox_inches="tight")
    plt.close()
    #->Plot the posterior probability distribution
    em.plot_corner(filename="result/{0}_triangle.png".format(targname), burnin=burnIn,
                   nuisance=nuisance, truths=parTruth, fraction=fraction,
                   quantiles=[psLow/100., psCenter/100., psHigh/100.], show_titles=True,
                   title_kwargs={"fontsize": 20})
    print("Post-processed!")
