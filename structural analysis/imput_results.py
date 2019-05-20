# -*- coding: utf-8 -*-
"""
Created on Fri May 17 10:09:04 2019

@author: Mels
"""

##import code from design_results:
import sys
sys.path.append("C:/Users/Mels/Desktop/3e jaar TUDelft/DSE/code/Bigger_is_Better/design_results")
f = open('Aerodynamic_concept.txt','r')
lines = f.readlines()
for line in lines:
    x = line.split(':')
    x_first = x[0]
    print(x_first)
    x_second = x[1]
    print(x_second)

#    for i in range(len(lines)):
#        line[i][1] = line[i][:-1]
