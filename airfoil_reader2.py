# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 11:42:41 2019

@author: floyd
"""
from matplotlib import pyplot as plt
cl = []
cd = []
file = open("C:/Users/floyd/Desktop/fus_coor.txt",'r')
for line in file:
     cl.append([(float(line.split()[0]))])
     cd.append([(float(line.split()[1]))])
#plt.plot(cl,cd) 
#plt.axis('equal')