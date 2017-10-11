from astropy.table import Table, vstack
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import numpy as np


gamma = 0.247


def Drude_profile(wave, b):
	I = b * gamma ** 2 / ((wave / 18 - 18 / wave) ** 2 + gamma ** 2)
	tau = 0.9 * I + 0.1 * (9.7 / wave) ** 1.7
	return tau


def residual(b, y, x):
	return y - Drude_profile(x, b)

t = Table.read('tau_lambda_kemper.txt', format='ascii')
t_use = t[t['col1'] > 12.7]

x = t_use['col1']
y = t_use['col2']
b0 = 10.
plsq = leastsq(residual, b0, args=(y, x))
b = plsq[0]
print(b)
x_new = np.logspace(1.2, 3, 1000)
y_new = Drude_profile(x_new, b)
plt.plot(x_new, y_new, label='fitting')
# plt.plot(t['col1'], t['col2'], label='model')
plt.plot(x, y, label='model')
plt.xlim([12, 50])
plt.legend()
plt.show()

upper = t['col1'][-1]
x = []
y = []
for i in range(len(t['col1'])):
	x.append(t['col1'][i])
	y.append(t['col2'][i])
wave = x_new[x_new > upper]
tau = y_new[x_new > upper]

for i in range(len(wave)):
	x.append(wave[i])
	y.append(tau[i])
#
# f = open('tau_lambda_kemper_new.txt', 'w')
# for i in range(len(x)):
# 	f.writelines('{0:.2f}\t{1:.4f}\n'.format(x[i], y[i]))
# f.close()
