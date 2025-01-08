# -*- coding: utf-8 -*-
"""
@author: MinjeongKang
"""
# code description: Growth curve fitting for a 48-multi-well plate using the Gompertz model.

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
def gompertz(t, asymptote, lag, mu):
    y = asymptote * np.exp(-lag * np.exp(-mu * t))
    return y
input_file = open("test.csv", 'r')
temp = input_file.readlines()
input_file.close()
time = []
OD_read = {}
for i in range(1,49):
    OD_read[i] = []
for i in temp:
    time.append(float(i.split(',')[0]))
    for j in range(1,48):
        OD_read[j].append(float(i.split(',')[j]))
    OD_read[48].append(float(i.split(',')[48][:-1]))

file = open("96well.txt", 'r')
temp = file.readlines()
file.close()
carbon = []
for i in temp:
    carbon.append(i.split('\n')[:-1])
  
numpy_time = np.array(time)
import matplotlib.gridspec as gridspec
p = 0
q = 0
fig = plt.figure(1, figsize=(120,80))
gs = gridspec.GridSpec(8, 12)
for i in range(1,49):
    y= []    
    for j in OD_read[i]:
        y.append(np.log(float(j))-np.log(float(OD_read[i][0])))
    q += 1
    if q>12:
        q = 1
        p += 1
    try:    
        popt, pcov = curve_fit(gompertz, numpy_time, y)
        bio_mu = popt[0]*popt[2]/np.exp(1)
        fitted_time = np.arange(time[0], time[-1], 0.1)
        fitted_y = []
        for k in fitted_time:
            fitted_y.append(gompertz(k, *popt))
        plt.ylim(0,5)
        ax1 = plt.subplot(gs[p, q-1])
        ax1.plot(fitted_time, fitted_y, color = 'red')
        fig.add_subplot(ax1)
    except:
        pass
    ax2 = plt.subplot(gs[p, q-1])
    ax2.plot(time, y, color = 'black')
    ax2.set_title(carbon[i-1])   
    fig.add_subplot(ax2)
    print (str(bio_mu))
   #plt.show()
plt.savefig("xdhB_glycerol.pdf")
