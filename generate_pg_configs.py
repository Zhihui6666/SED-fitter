from astropy.table import Table

# I fix the extinction to a very small value for PG quasars
rl_table = Table.read('data_PG/targlist_clu_rl.ipac', format='ascii.ipac')
rq_table = Table.read('data_PG/targlist_clu_rq.ipac', format='ascii.ipac')

f = open('configs/CGCG011-076_3_for_pg.py')
content = f.readlines()
f.close()
# cat3d_G_and_extinction = content[107: 158]
cat3d_H_and_extinction = content[147: 198]
rl = open('data_PG/config_clu_rl.py')
rl_content = rl.readlines()
rl.close()
rq = open('data_PG/config_clu_rq.py')
rq_content = rq.readlines()
rq.close()
head = rl_content[:8]
mid = rl_content[10:46]
rl_end = rl_content[117:]
rq_end = rq_content[117:]

i = 0
for targname in rl_table['Name']:
	f = open('configs_PG_H/{}.py'.format(targname), 'w')
	f.writelines(head)
	f.writelines('targname = \"{}\"\n'.format(targname))
	f.writelines('redshift = {0:.3f}\n'.format(rl_table['z'][i]))
	f.writelines(mid)
	f.writelines(cat3d_H_and_extinction)
	f.writelines(rl_end)
	f.close()
	i += 1

i = 0
for targname in rq_table['Name']:
	f = open('configs_PG_H/{}.py'.format(targname), 'w')
	f.writelines(head)
	f.writelines('targname = \"{}\"\n'.format(targname))
	f.writelines('redshift = {0:.3f}\n'.format(rq_table['z'][i]))
	f.writelines(mid)
	f.writelines(cat3d_H_and_extinction)
	f.writelines(rq_end)
	f.close()
	i += 1
