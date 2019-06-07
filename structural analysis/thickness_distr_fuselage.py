# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:52:58 2019

@author: Mathilde
"""
import numpy as np
from matplotlib import pyplot as plt

R = 2. #radius of the fuselage
sigma = 483*10**6 #ultimate stress of aluminium 
M = 60000000 #combined moment of Mz and My on fuselage 
theta = list(np.array([-15,-10,-5,0,5,10,15])*np.pi/180) #angle between the moment and the z-axis

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
    
#Iteration 2   
  
beta = 5
n_segments = 360/beta
segment_length = np.ones((1,int(n_segments)))*(beta/360)*2*np.pi*R
alpha_segment = np.arange(beta/2,360,beta)*np.pi/180


t_iter=[]
for i in range(len(alpha_segment)):
    if alpha_segment[i]>=0 and alpha_segment[i]<=(0.5*np.pi+min(theta)):
        t = abs(t_max*np.sin(alpha_segment[i]+max(theta)))
        t_iter.append(t)
    if (0.5*np.pi+min(theta))< alpha_segment[i] and (0.5*np.pi+max(theta))>=alpha_segment[i]:
        t = t_max
        t_iter.append(t)
    if (0.5*np.pi+max(theta))< alpha_segment[i] and np.pi > alpha_segment[i]:
        t = abs(t_max*np.sin(alpha_segment[i]+min(theta)))
        t_iter.append(t)
    if np.pi < alpha_segment[i] and (1.5*np.pi+min(theta))> alpha_segment[i]:
        t = abs(t_max*np.sin(alpha_segment[i]+max(theta)))
        t_iter.append(t)
    if (1.5*np.pi+min(theta))< alpha_segment[i] and (1.5*np.pi+max(theta))> alpha_segment[i]:
        t = t_max
        t_iter.append(t)
    if (1.5*np.pi+max(theta))< alpha_segment[i] and 2*np.pi>= alpha_segment[i]:
        t = abs(t_max*np.sin(alpha_segment[i]+min(theta)))
        t_iter.append(t)       
        
t_iter = np.array(t_iter)/1000
#print(t_iter) 
plt.plot(alpha_segment*180/np.pi,t_iter*1000)
#plt.show()
#Iteration 2 of variable thickness

I_xx = 0
I_yy = 0

def max_stress(theta,alpha_segment,I_xx,I_yy):
    stress_ratios = np.zeros((len(alpha_segment), len(theta)))
    max_stress_ratios = []
    for i in range(len(theta)): 
        for j in range(len(alpha_segment)):    
            stress_ratio = (M*np.cos(theta[i])*R*np.sin(alpha_segment[j])/I_xx + -M*np.sin(theta[i])*R*np.cos(alpha_segment[j])/I_yy)/sigma
            stress_ratios[j][i] = stress_ratio
    
    for i in range(len(stress_ratios)):
        max_stress_ratios.append(max(abs(stress_ratios[i])))
    max_stress_ratios = np.asarray(max_stress_ratios)
    return(max_stress_ratios)

finalthicknessreached = False
c1ounter = 0
while finalthicknessreached == False:# and counter < 100:
    I_xx = t_iter*segment_length[0]**3*np.sin(alpha_segment+0.5*np.pi)**2/12 + t_iter*segment_length[0]*R**2*np.sin(alpha_segment)**2
    I_yy = t_iter*segment_length[0]**3*np.cos(alpha_segment+0.5*np.pi)**2/12 + t_iter*segment_length[0]*R**2*np.cos(alpha_segment)**2
    I_xx = sum(I_xx)
    I_yy = sum(I_yy) 
    stress_ratio = max_stress(theta,alpha_segment,I_xx,I_yy)
    
    factor = 1/(counter+1)
    t_new = t_iter*stress_ratio*(1+0.2*factor)
    t_iter = t_new
#    counter += 1
#    print(counter)
    print(stress_ratio[0])
    print(t_iter[0])
#    
    for i in stress_ratio:
#        print(i)
        if i<=1.01:
            finalthicknessreached = True
        else:
            finalthicknessreached = False
            break

plt.plot(alpha_segment*180/np.pi,t_iter*1000)
plt.show()



#max_stress_ratios = np.zeros((1,len(t_iter)))

#for i in range(len(alpha_segment)):
#    I_xx += t_iter[i]*segment_length[0][i]**3*np.sin(alpha_segment[i]+0.5*np.pi)**2/12 + t_iter[i]*segment_length[0][i]*R**2*np.sin(alpha_segment[i])**2
#    I_yy += t_iter[i]*segment_length[0][i]**3*np.cos(alpha_segment[i]+0.5*np.pi)**2/12 + t_iter[i]*segment_length[0][i]*R**2*np.cos(alpha_segment[i])**2
#
#    stress_ratios = np.zeros((len(alpha_segment), len(theta)))
#        
#for i in range(len(theta)): 
#    for j in range(len(alpha_segment)):    
#        stress_ratio = (M*np.cos(theta[i])*R*np.sin(alpha_segment[j])/I_xx + -M*np.sin(theta[i])*R*np.cos(alpha_segment[j])/I_yy)/sigma
#        stress_ratios[j][i] = stress_ratio
#    #    print(np.shape(stress_ratios))
#    #    print(stress_ratios)
##max_stress_ratios=[]
##    print(t_old)
#print(t_iter[0:5])
#
#for i in range(len(stress_ratios)):
#    max_stress_ratios.append(max(abs(stress_ratios[i])))
#
#    t_new = t_iter[i]*max(abs(stress_ratios[i]))    
#    t_iter[i] = t_new
#    
#print(max_stress_ratios[0:5])   
#print(t_iter[0:5])
#finalthicknessreached = False


#for i in range(len(max_stress_ratios)):
#    I_xx = 0
#    I_yy = 0
#    
#    while finalthicknessreached == False:
#        
#        for p in range(len(alpha_segment)):
#            I_xx += t_iter[p]*segment_length[0][p]**3*np.sin(alpha_segment[p]+0.5*np.pi)**2/12 + t_iter[p]*segment_length[0][p]*R**2*np.sin(alpha_segment[p])**2
#            I_yy += t_iter[p]*segment_length[0][p]**3*np.cos(alpha_segment[p]+0.5*np.pi)**2/12 + t_iter[p]*segment_length[0][p]*R**2*np.cos(alpha_segment[p])**2
#
#        if 0.92>abs(max_stress_ratios[i]) or abs(max_stress_ratios[i])>0.99:
##            stress_ratios = np.zeros((len(alpha_segment), len(theta)))
##            for k in range(len(theta)): 
#            for j in range(len(alpha_segment)):#len(alpha_segment)):    
#                    stress_ratio = (M*np.cos(theta[k])*R*np.sin(alpha_segment[j])/I_xx + -M*np.sin(theta[k])*R*np.cos(alpha_segment[j])/I_yy)/sigma
#                    stress_ratios[j][k] = stress_ratio
#    
#            max_stress_ratios[i]= max(abs(stress_ratios[i]))
#            t_new = t_iter[i]*max(abs(stress_ratios[i]))    
#            t_iter[i] = t_new
#            print(max_stress_ratios[i])
#            print("thickness",t_iter[i])
#        
#        elif abs(max_stress_ratios[i])>= 0.92 and abs(max_stress_ratios[i])<= 0.99:
#            finalthicknessreached = True