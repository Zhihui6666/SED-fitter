from glob import glob
import pandas as pd

redshift = pd.read_csv('redshift.csv')

model = open('configs_test/F00073+2538.py')
content = model.readlines()
model.close()
head = content[:7]
body = content[11:]
print(head)
print(body)

for i in range(len(redshift)):
	targname = redshift.loc[i, 'IRAS Name']
	f = open('configs_test/{}.py'.format(targname), 'w')
	f.writelines(head)
	f.writelines('targname = \"{}\"\n'.format(targname))
	print('targname = \"{}\"\n'.format(targname))
	f.writelines('redshift = {0:.6f}\n'.format(redshift.loc[i, 'z']))
	print('redshift = {0:.6f}\n'.format(redshift.loc[i, 'z']))
	f.writelines('distance = None\n')
	print('distance = None\n')
	f.writelines('sedFile = \"data_test/{}.csv\"\n'.format(targname))
	print('sedFile = \"data_test/{}.csv\"\n'.format(targname))
	f.writelines(body)
	f.close()