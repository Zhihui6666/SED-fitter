import glob
import numpy as np
import pandas as pd
from pprint import pprint

# output of chi square and fraction

# namelist = []
# for filename in glob.glob("data/*.csv"):
#     filename = filename.lstrip("data/").rstrip(".csv")
#     namelist.append(filename)
#
#
# data = pd.DataFrame(np.zeros((len(namelist), 3)) ,index=namelist, columns=["2-component", "3-component", "4-component"])
#
# for filename in glob.glob("configs/*.py"):
#     filename = filename.lstrip("configs/").rstrip(".py")
#     index = namelist.index(filename[:-2])
#     f = open("result/"+filename+"_bestfit.txt")
#     num1 = 100; num2 = 100
#     content = f.readlines()
#     for i in range(len(content)):
#         if content[i].find("reduced") != -1:
#             num1 = i
#             break
#     for i in range(len(content)):
#         if content[i].find("total") != -1:
#             num2 = i
#             break
#     # print num1, num2, index
#     string = ""
#     for i in content[min(num2, num1):]:
#         string += i
#     if filename[-1] == "2":
#         data.iloc[index, 0] = string
#     elif filename[-1] == "3":
#         data.iloc[index, 1] = string
#     elif filename[-1] == "4":
#         data.iloc[index, 2] = string
# data = data.sort_index()
# data.to_csv("output.csv")


# output of bestfit parameters

col = []
f = open("result/CGCG011-076_4_bestfit.txt")
content = f.readlines()
for i in content:
    i = i.strip("\n")
    i = i.rstrip(" ")
    i = i.split(" ")
    if i[0] != "Name" and i[0] != "reduced" and i[0] != "total" and i[0] != "fraction" and i[0] != "" and i[0] != "#lnL_max:":
        col.append(i[0])
for i in content:
    if i[:7] == "reduced":
        i = i.split(": ")
        col.append(i[0])
    elif i.find("=") != -1:
        i = i.split(" = ")
        col.append(i[0])
f.close()

objlist = []
for filename in glob.glob("configs/*.py"):
    filename = filename.lstrip("configs/").rstrip(".py")
    objlist.append(filename)

bestfitPar = pd.DataFrame(index=objlist, columns=col)

for filename in objlist:
    f = open("result/"+filename+"_bestfit.txt")
    index = objlist.index(filename)
    content = f.readlines()
    for i in content:
        i = i.strip("\n")
        i = i.rstrip(" ")
        i = i.split(" ")
        if i[0] != "Name" and i[0] != "reduced" and i[0] != "total" and i[0] != "fraction" and i[0] != "" and i[
            0] != "#lnL_max:":
            bestfitPar.iloc[index,col.index(i[0])] = float(i[-1])
    for i in content:
        if i[:7] == "reduced":
            i = i.split(": ")
            bestfitPar.iloc[index, col.index(i[0])] = round(float(i[-1]), 2)
        elif i.find("=") != -1:
            i = i.split(" = ")
            bestfitPar.iloc[index, col.index(i[0])] = str(round(float(i[-1]), 4) * 100) + "%"
    f.close()

bestfitPar = bestfitPar.sort_index()
# bestfitPar.to_csv("bestfitPar.csv")

