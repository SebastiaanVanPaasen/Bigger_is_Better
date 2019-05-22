# -*- coding: utf-8 -*-
"""
Created on Mon May 20 21:05:25 2019

@author: floyd
"""
from matplotlib import pyplot as plt
import numpy as np

cr = np.array([6.1, 6.62, 7.16, 8.62, 4.87, 7.88, 12.876])
ct = np.array([1.89,2.65,2.22,3.45, 1.95, 1.25, 1.919])
b = np.array([55.87, 41.7,39.87, 48.28, 61.3,34.32, 60.90])
sweep = np.array([27.69, 2.73,28.77, 3.07, 1.36, 29.39, 35.19])
color = ['C0','C1','C2','C3','C4', 'C5', 'C6']
labels = ["AERO", "DD2E", "DD4E", "HBPE", "STRW", "737-800", "777-200"]
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