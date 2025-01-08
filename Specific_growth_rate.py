# -*- coding: utf-8 -*-
"""
@author: Minjeong Kang
"""
# code description: Growth curve fitting and calculation of specific growth rate using the Gompertz model for data from inosine media.

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
def gompertz(t, asymptote, lag, mu):
    y = asymptote * np.exp(-lag * np.exp(-mu * t))
    return y
input_file = open("xdhB_WT_inosine media.csv", 'r')
temp = input_file.readlines()
input_file.close()
time = []
OD_read1 = []
for i in temp:
    time.append(float(i.split(',')[0]))
    OD_read1.append(float(i.split(',')[1][:-1]))

numpy_time = np.array(time)
initial_OD = OD_read1[0]
L = []
y = []
for i in OD_read1:
    y.append(np.log(float(i))-np.log(float(initial_OD)))

popt, pcov = curve_fit(gompertz, numpy_time, y)

bio_mu = popt[0]*popt[2]/np.exp(1)
fitted_time = np.arange(time[0], time[-1], 0.1)
fitted_y = []
for i in fitted_time:
    fitted_y.append(gompertz(i, *popt))
plt.figure(figsize=(10,5))
#plt.plot(time, OD_read)
#plt.show()

plt.plot(fitted_time, fitted_y, color = 'red')
plt.plot(time, y, color = 'black')
plt.show()
plt.savefig("sgr")
print ('Specific growth rate is ' + str(bio_mu))
print (popt)
