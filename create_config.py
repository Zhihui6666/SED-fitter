import glob
from sedfit.SED_Toolkit import Load_SED

# write config file
# add a condition if no IRS spectrum, use photometry points in MIR

# a = Load_SED("data/CGCG049-057.csv")
# print(a["spc"][0])

# a = open("configs/CGCG011-076_2.py", "r")
# b = open("configs/CGCG011-076_3.py", "r")
# c = open("configs/CGCG011-076_4.py", "r")
# content = a.readlines()
# beginning = content[:10]
# print(beginning)
# datadict = content[15:39]
# body1 = content[39:]
# body2 = b.readlines()[34:]
# body3 = c.readlines()[34:]
# chainlength = content[-53:]
# print(chainlength)
# lna = content[-94]
# a.close()
# b.close()
# c.close()
#
# for filename in glob.glob("configs/*.py"):
# 	filename = filename.lstrip("configs/")
# 	if filename[:-5] != "configs/CGCG011-076":
# 		a = open(filename, "r")
# 		content = a.readlines()
# 		targetinfo1 = content[10:13]
# 		targetinfo2 = content[14]
# 		a.close()
		# temp = open(filename, "w")
# 		temp.writelines(beginning)
# 		temp.writelines(targetinfo1)
# 		temp.writelines(targetinfo2)
# 		temp.writelines(datadict)
# 		if filename[-4] == "2":
# 			temp.writelines(body1)
# 		elif filename[-4] == "3":
# 			temp.writelines(body2)
# 		else:
# 			temp.writelines(body3)
# 		temp.writelines(content[:-53])
# 		temp.writelines(chainlength)
# 		temp.close()
# 	a = open(filename, "r")
# 	content = a.readlines()
# # 	# print(content[14])
# # 	b = content[:-53]
# 	a.close()
# 	# print(content[:-94])
# 	# print(content[-93:])
# 	temp = open(filename, "w")
# 	# temp.writelines(b)
# 	# temp.writelines(chainlength)
# 	temp.writelines(content[:-94])
# 	temp.writelines(lna)
# 	temp.writelines(content[-93:])
# 	temp.close()

# write .sh file

# f = open("1.sh", "w")
# for filename in glob.glob("configs/*.py"):
# 	filename = filename.lstrip("configs/")
# 	f.writelines("python UniFit.py configs/"+filename+" -m 44\n")
# f.close()


# Deal with newCLUMPY model
f = open('configs/CGCG011-076_3_for_pg.py')
content = f.readlines()
f.close()
new = content[107:146] + content[187:]
print(new)

for filename in glob.glob('configs/*_3.py'):
	f = open(filename)
	content = f.readlines()
	f.close()
	info = content[:107]
	f = open(filename, 'w')
	f.writelines(info)
	f.writelines(new)
	f.close()
print('Done!')
