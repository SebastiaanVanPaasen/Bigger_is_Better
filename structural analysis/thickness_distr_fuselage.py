# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:52:58 2019

@author: Mathilde
"""
import numpy as np
from matplotlib import pyplot as plt
import math as m


def t_fus_due_to_moment(sigma, Mx, My):
    
    x_pos = []
    y_pos = []
    file = open("C:/Users/mathi/Documents/DSE/Bigger_is_Better/structural analysis/fus_coor.txt", 'r')
    for line in file:
        x_pos.append([(float(line.split()[0]))])
        y_pos.append([(float(line.split()[1]))])
    x_pos = np.array(x_pos)
    y_pos = np.array(y_pos)    
    R_2 = (x_pos**2 + y_pos**2)**0.5/1000
#    print(R_2)

    theta = list(np.array([-30, -15, 0, 15, 30])*np.pi/180) #angle between the moment and the z-axis
    alpha_segment = []
    #Iteration 1 of the variable thickness
    for i in range(len(y_pos)):
        if y_pos[i]>0 and x_pos[i]>0:
            alpha_segment.append(np.arctan(y_pos[i]/x_pos[i]))
        elif y_pos[i]>0 and x_pos[i]<0:
            alpha_segment.append(np.arctan(-x_pos[i]/y_pos[i]) + np.pi/2)
        elif y_pos[i]<0 and x_pos[i]<0:
            alpha_segment.append(np.arctan(y_pos[i]/x_pos[i]) + np.pi)
        elif y_pos[i]<0 and x_pos[i]>0:
            alpha_segment.append(np.arctan(x_pos[i]/-y_pos[i]) + np.pi*3/2)#*180/np.pi#*np.pi/180
    
    segment_length = np.zeros((1, len(alpha_segment)))

    alpha_segment = np.array(alpha_segment)
  
    t_circ = np.zeros((len(theta), len(alpha_segment)))
    t_max_all = []
    
    for i in range(len(theta)):
        
        for j in range(len(alpha_segment)):
            t = ((Mx * np.cos(theta[i]) + My * np.sin(theta[i]))*np.sin(alpha_segment[j]) + (My * np.cos(theta[i]) - Mx * np.sin(theta[i]))*np.cos(alpha_segment[j]))/(np.pi*R_2[j]**2*sigma)
            t_circ[i][j] = abs(t[0])
     
    
    for i in range(len(alpha_segment)):
        segment_values = []
        for j in range(len(theta)):
            segment_values.append(t_circ[j][i])
        
        t_max = max(segment_values)
        t_max_all.append(t_max)
   
    for i in range(len(R_2)-1):
        
        segment_length[0][i] = alpha_segment[i]*(R_2[i+1][0]+R_2[i][0])/2
#        print(segment_length)
#    print(R_2[0])
    segment_length[0][-1] = alpha_segment[-1] * (R_2[0][0] + R_2[-1][0])/2
#    print(segment_length)
#  print(t_circ[1]*1000) 
    plt.plot(alpha_segment*180/np.pi,t_circ[0]*1000,label = 'theta = -30')
    plt.plot(alpha_segment*180/np.pi,t_circ[1]*1000, label = 'theta = -15')
    plt.plot(alpha_segment*180/np.pi,t_circ[2]*1000, label = 'theta = 0')
    plt.plot(alpha_segment*180/np.pi,t_circ[3]*1000, label = 'theta = 15')
    plt.plot(alpha_segment*180/np.pi,t_circ[4]*1000, label = 'theta = 30')
    plt.plot(alpha_segment*180/np.pi,np.array(t_max_all)*1000, label = 'maximum thickness')

    plt.xlabel("Alpha [degrees]")
    plt.ylabel("Thickness distribution [mm]")

    plt.legend()
    plt.show()
    return(t_max_all, alpha_segment, R_2, segment_length, theta)
    #Iteration 2 of variable thickness
    
I_xx = 0
I_yy = 0
    
def max_stress(theta,alpha_segment,I_xx,I_yy, Mx, My):
    stress_ratios = np.zeros((len(alpha_segment), len(theta)))
    max_stress_ratios = []
    for i in range(len(theta)): 
        for j in range(len(alpha_segment)):    
            stress_ratio = ((Mx * np.cos(theta[i]) + My * np.sin(theta[i])) * R_2[j][0] * np.sin(alpha_segment[j]) / I_xx + ( My * np.cos(theta[i]) - Mx * np.sin(theta[i]))*R_2[j][0]*np.cos(alpha_segment[j])/I_yy)/sigma
            stress_ratios[j][i] = stress_ratio
    
    for i in range(len(stress_ratios)):
        max_stress_ratios.append(max(abs(stress_ratios[i])))
    
    max_stress_ratios = np.asarray(max_stress_ratios)
    
    return(max_stress_ratios)

finalthicknessreached = False
counter = 0

sigma = 105*10**6 #ultimate stress of aluminium 
Mx = -37799412
My = 216.17*10**3*7

t_iter, alpha_segment, R_2, segment_length, theta = t_fus_due_to_moment(sigma, Mx, My)

while finalthicknessreached == False:# and counter < 100:
    I_xx = t_iter*segment_length[0]**3*np.sin(alpha_segment[0]+0.5*np.pi)**2/12 + t_iter*segment_length[0]*R_2[0]**2*np.sin(alpha_segment[0])**2
    I_yy = t_iter*segment_length[0]**3*np.cos(alpha_segment[0]+0.5*np.pi)**2/12 + t_iter*segment_length[0]*R_2[0]**2*np.cos(alpha_segment[0])**2
    I_xx_tot = sum(I_xx)
    I_yy_tot = sum(I_yy) 
    stress_ratio = max_stress(theta,alpha_segment,I_xx_tot,I_yy_tot, Mx, My)
    
    factor = 1/(counter+1)
    t_new = t_iter*stress_ratio#*(1+0.2*factor)
    t_iter = t_new
#    counter += 1
#    print(counter)
    print(stress_ratio)
    print(t_iter[0])
#    
    for i in stress_ratio:
#        print(i)
        if 0.5 < i <= 1.01:
            finalthicknessreached = True
        else:
            finalthicknessreached = False
            break
        
plt.figure(2)
plt.plot(alpha_segment*180/np.pi, t_iter*1000)
plt.show()




sigma_fatigue_hoop = 350 * 10**6 # look up from graph
sigma_fatigue_long = 300 * 10**6
internal_p = 78.2 * 10**3 
external_p = 30.1 * 10**3
R = 5.95/2

def t_fus_due_to_pressure(sigma_fatigue_hoop, sigma_fatigue_long,R,internal_p, external_p):
    delta_p = internal_p - external_p
    t_min_hoop = delta_p*R/sigma_fatigue_hoop
    t_min_long = delta_p*R/(2*sigma_fatigue_long)
    return t_min_hoop, t_min_long

t_min_hoop,t_min_long = t_fus_due_to_pressure(sigma_fatigue_hoop, sigma_fatigue_long,R,internal_p, external_p)
#print(t_min_hoop)
#print(t_min_long)


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