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

    
#Iteration 1 of the variable thickness
alpha = np.arange(0,361,1) # angle between the z-axis and an point evaluated in circumference of cross-section
#first approximation of the fuselage Izz, assuming circular with constant thickness 
#I_zz = np.pi*R**3*t

for i in range(len(theta)):
    
    t_circ = []
    for j in range(len(alpha)):
        t = M * np.sin(alpha[j]*np.pi/180-theta[i])/(np.pi*R**2*sigma)
        t_circ.append(abs(t*1000))
        
    t_max = max(t_circ)   
    #stability_trend = np.polyfit(x_cg, Sh_S_stability, 1)
    plt.plot(alpha,t_circ)
    plt.show()
    
#Iteration 2   
    
beta = 5
n_segments = 360/beta
segment_length = np.ones((1,int(n_segments)))*(beta/2*np.pi)*2*np.pi*R
alpha_segment = np.arange(beta/2,360,beta)*np.pi/180

print(alpha_segment)

t_iter_1=[]
for i in range(len(alpha_segment)):
    if alpha_segment[i]>0 and alpha_segment[i]<(0.5*np.pi+min(theta)):
        t = abs(t_max*np.sin(alpha_segment[i]+max(theta)))
        t_iter_1.append(t)
    if (0.5*np.pi+min(theta))< alpha_segment[i] and (0.5*np.pi+max(theta))>alpha_segment[i]:
        t = t_max
        t_iter_1.append(t)
    if (0.5*np.pi+max(theta))< alpha_segment[i] and np.pi > alpha_segment[i]:
        t = abs(t_max*np.sin(alpha_segment[i]+max(theta)))
        t_iter_1.append(t)
    if np.pi < alpha_segment[i] and (1.5*np.pi+min(theta))> alpha_segment[i]:
        t = abs(t_max*np.sin(alpha_segment[i]+max(theta)))
        t_iter_1.append(t)
    if (1.5*np.pi+min(theta))< alpha_segment[i] and (1.5*np.pi+max(theta))> alpha_segment[i]:
        t = t_max
        t_iter_1.append(t)
    if (1.5*np.pi+max(theta))< alpha_segment[i] and 2*np.pi> alpha_segment[i]:
        t = abs(t_max*np.sin(alpha_segment[i]+max(theta)))
        t_iter_1.append(t)       
        
print(t_iter_1)      
plt.plot(alpha_segment,t_iter_1)
plt.show()
#Iteration 2 of variable thickness


for i in range(int(n_segments)):
    I_xx += t[i]*segment_length[i]*np.sin(alpha_segment[i]+0.5*np.pi)**2/12 + t[i]*segment_length[i]*R**2*np.sin(alpha_segment[i])**2
    I_yy += t[i]*segment_length[i]*np.cos(alpha_segment[i]+0.5*np.pi)**2/12 + t[i]*segment_length[i]*R**2*np.cos(alpha_segment[i])**2
    
    
    
    