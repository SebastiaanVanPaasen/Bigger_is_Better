# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:52:58 2019

@author: Mathilde
"""
import numpy as np
from matplotlib import pyplot as plt

R = 4. #radius of the fuselage
sigma = 502*10**6 #ultimate stress of aluminium 
M =60000000 #combined moment of Mz and My on fuselage 
theta = list(np.array([-15,-10,-5,0,5,10,15])*np.pi/180) #angle between the moment and the z-axis

#define the circumference of the fuselage

s = 2*np.pi*R
N = 200 #number of point evaluated
ds = s/N
z_pos = []
y_pos = []
z = 0.
y = 0.

    

alpha = np.arange(0,361,1) # angle between the z-axis and an point evaluated in circumference of cross-section

#first approximation of the fuselage Izz, assuming circular with constant thickness 
#I_zz = np.pi*R**3*t

for i in range(len(theta)):
    
    t_circ = []
    for j in range(len(alpha)):
        t = M * np.sin(alpha[j]*np.pi/180-theta[i])/(np.pi*R**2*sigma)
        t_circ.append(abs(t*1000))
        
    print(max(t_circ))   
    #stability_trend = np.polyfit(x_cg, Sh_S_stability, 1)
    plt.plot(alpha,t_circ)
    plt.show()
    
    
    
    
    