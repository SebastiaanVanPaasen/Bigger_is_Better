# -*- coding: utf-8 -*-
"""
Created on Mon May 20 21:05:25 2019

@author: floyd
"""
from matplotlib import pyplot as plt
import numpy as np

cr = np.array([6.15, 5.85, 6.04, 6.36, 5.53, 5.71, 11.11])
ct = np.array([1.9, 1.81 ,1.8 ,1.89, 1.71, 1.56, 2.94])
b = np.array([52.36, 49.75, 50.95, 61.93, 50.66 ,34.32, 60.90])
sweep = np.array([27.69, 27.69,31.07, 31.07, 27.69, 29.39, 35.19])
color = ['C0','C1','C2','C3','C4', 'C5', 'C6']
labels = ["HIGH SDD", "HIGH DD", "LOW SDD", "HIGH DD STRUT", "HIGH SDD STRUT", "737-800", "777-200"]
plt.figure()
for i in range(len(cr)):
    coordinates = [0, b[i]/2*np.tan(np.radians(sweep[i])),b[i]/2*np.tan(np.radians(sweep[i]))+ct[i],cr[i]]
    y = [0,b[i]/2,b[i]/2,0]
    for j in range(len(coordinates)-1):
        if j == 1:
            mylabel = labels[i]
        else:
            mylabel = ""
        line = [coordinates[j],coordinates[j+1]]
        line2 = [y[j],y[j+1]]

        plt.plot(line2,line,color=color[i], label=mylabel)
plt.legend()
plt.xlabel("Half span [m]")
plt.ylabel("Length [m]")
plt.axis('scaled')
plt.grid()
plt.show()