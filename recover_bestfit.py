import glob

for filename in glob.glob("result/*bestfit.txt"):
	f = open(filename, "r")
	content = f.readlines()
	f.close()
	for i in range(len(content)):
		if content[i][:8] == "#lnL_max":
			num = i
	f = open(filename, 'w')
	f.writelines(content[:num+1])
	f.close()