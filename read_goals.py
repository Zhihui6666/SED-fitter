import pandas
from pprint import pprint
import glob

f = open("ALLGOALS.csv")
content = f.readlines()
f.close()
list = []
namelist = []

for i in content:
    list.append(i.rstrip("\r\n").split(","))

for i in list:
    if i[36][0] != "F":
        i[36] = "F"+i[36]
    namelist.append(i[36])

# write down the data

# df = pandas.read_csv("data(old)/F00073+2538.csv")
# num = 0
# for name in namelist:
#     # JHK band
#     for i in range(3):
#         df.iloc[i, 1] = float(list[num][2*i])
#         df.iloc[i, 2] = float(list[num][2*i+1])
#
#     # WISE_1 - 4
#     df.iloc[3, 1] = float(list[num][6])
#     df.iloc[3, 2] = float(list[num][7])
#     df.iloc[6, 1] = float(list[num][8])
#     df.iloc[6, 2] = float(list[num][9])
#     df.iloc[9, 1] = float(list[num][10])
#     df.iloc[9, 2] = float(list[num][11])
#     df.iloc[10, 1] = float(list[num][12])
#     df.iloc[10, 2] = float(list[num][13])
#
#     # Herschel 1 - 6
#     for i in range(6):
#         df.iloc[i+12, 1] = float(list[num][2*i+14])
#         df.iloc[i+12, 2] = float(list[num][2*i+15])
#
#
#     # IRAC 1 - 4
#     df.iloc[4, 1] = float(list[num][26])
#     df.iloc[4, 2] = float(list[num][27])
#     df.iloc[5, 1] = float(list[num][28])
#     df.iloc[5, 2] = float(list[num][29])
#     df.iloc[7, 1] = float(list[num][30])
#     df.iloc[7, 2] = float(list[num][31])
#     df.iloc[8, 1] = float(list[num][32])
#     df.iloc[8, 2] = float(list[num][33])
#
#     # MIPS_24
#     df.iloc[11, 1] = float(list[num][34])
#     df.iloc[11, 2] = float(list[num][35])
#
#     df.to_csv("data/"+str(name)+".csv", index=False)
#     df.close()
#     num += 1

# write the config file
model = open("configs(old)/F00073+2538_4.py")
content = model.readlines()
model.close()
num = 0
for name in namelist:
    content[10] = "targname = " + "\"" + str(name) + "\"" + "\n"
    content[11] = "redshift = " + list[num][-2] + "\n"
    content[12] = "distance = " + list[num][-1] + "\n"
    content[14] = "sedFile = "+ "\"data/" + str(name) + ".csv" + "\"" + "\n"
    content[235] = "        (\"iteration\", [1000, 1000, 600])," + "\n"
    num += 1
    f2 = open("configs/"+str(name)+"_2.py", "w")
    f3 = open("configs/"+str(name)+"_3.py", "w")
    f4 = open("configs/"+str(name)+"_4.py", "w")
    f2.writelines(content[:102]+content[175:])
    f3.writelines(content[:156]+content[175:])
    f4.writelines(content)
    f2.close()
    f3.close()
    f4.close()

# f = open("1.sh", "w")
# for filename in glob.glob("configs/*.py"):
#     filename = filename
#     f.writelines("python UniFit.py "+filename+" -m 64"+"\n")
# f.close()