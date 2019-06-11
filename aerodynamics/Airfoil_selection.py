#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 11:11:06 2019

@author: Max
"""

import matplotlib.pyplot as plt
import numpy as np

f = open('SC(2)-0612.txt', 'r')
lines = f.readlines()

f.close
SC12 = []
for line in lines: 
    x = line.split()
    SC12.append(x)
    for i in range(len(SC12)):
        for j in range(len(SC12[i])):
            SC12[i][j] = float(SC12[i][j])


f = open('SC(2)-0614.txt', 'r')
lines = f.readlines()

f.close
SC14 = []
for line in lines: 
    x = line.split()
    SC14.append(x) 
    for i in range(len(SC14)):
        for j in range(len(SC14[i])):
            SC14[i][j] = float(SC14[i][j])



for i in range(len(SC12)):
    plt.plot(SC12[i][0],SC12[i][1], 'xb')

    
SC16 = []
for i in range(len(SC14)):
    SC16.append([SC14[i][0], SC14[i][1]*(16/14)])


print(SC16)

for i in range(len(SC14)):
    plt.plot(SC14[i][0],SC14[i][1], 'ob')
    
    
    
#for i in range(len(SC12_14)):
#    plt.plot(SC12_14[i][0],SC12_14[i][1], 'xr')    
#    


#with open('SC(2)-0616.txt', 'w') as f:
#    for i in SC16:
#        f.write( "%s\n" % i)

   
np.savetxt('SC(2)-0616.txt', SC16)     
        

    
plt.show()    
print (SC12, SC14)    