import os

clumpy_path = "/Users/zhihui/SEDfitting/Fitter_new/template/"
#->Obtain the current path
pathList = os.path.abspath(__file__).split("/")
#->Create the path to the filters
pathList[-1] = "filters/"
filter_path = "/".join(pathList)
#->Create the path to the templates
pathList[-2] = "template/"
template_path = "/".join(pathList[0:-1])
#template_path = "/Users/jinyi/Work/mcmc/Fitter/template/"
