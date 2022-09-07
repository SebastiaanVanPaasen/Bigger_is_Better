# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 16:10:44 2019

@author: mathi
"""

import thickness_distr_fuselage as tf
import numpy as np
from matplotlib import pyplot as plt


x_pos = []
y_pos = []
file = open("C:/Users/mathi/Documents/DSE/Bigger_is_Better/structural analysis/fus_coor.txt", 'r')
for line in file:
    x_pos.append([(float(line.split()[0]))])
    y_pos.append([(float(line.split()[1]))])
x_pos = np.array(x_pos[::-1])/1000
y_pos = -np.array(y_pos[::-1])/1000

#plt.figure()
#plt.plot(x_pos,y_pos)
#plt.show()

def fuselage_centroid(boom_area, y_pos, x_pos):
    
    y_pos = np.reshape(y_pos, (1,len(boom_area)))
    x_pos = np.reshape(x_pos, (1,len(boom_area)))

#    print(np.shape(y_pos))
    y_centroid = sum(boom_area*y_pos[0])/sum(boom_area)
    x_centroid = sum(boom_area*x_pos[0])/sum(boom_area)
    
    return y_centroid, x_centroid
    
#set_boom_area = 0.02

def fuselage_booms(sigma, Mx, My):
    
    t_max, alpha_segment, R, segment_length, theta, t_circ = tf.t_fus_due_to_moment(sigma, Mx, My)
    t_iter = np.array(tf.t_iter)
    req_boom_area = t_iter*segment_length
    actual_boom_area = req_boom_area - (segment_length*0.002)
    
    for i in range(len(actual_boom_area)):
        if actual_boom_area[i] <= 0.001:
            actual_boom_area[i] = segment_length[i]*0.002
        elif actual_boom_area[i] > 0.001: 
            actual_boom_area[i] = 0.005
#        elif actual_boom_area[i] > 0 and actual_boom_area[i] < 0.001 : 
#            actual_boom_area[i] = 0.001
#        elif actual_boom_area[i] > 0.001 and actual_boom_area[i] < 0.002 : 
#            actual_boom_area[i] = 0.002
#        elif actual_boom_area[i] > 0.002 and actual_boom_area[i] < 0.003 : 
#            actual_boom_area[i] = 0.003
#        elif actual_boom_area[i] > 0.003 and actual_boom_area[i] < 0.004 : 
#            actual_boom_area[i] = 0.004
        else:
            actual_boom_area[i] = 0.0008
    
    return actual_boom_area, t_iter




def fuselage_inertia(boom_area, y_centroid, x_centroid, y_pos, x_pos, My, Mx, sigma):

    t_max, alpha_segment, R, segment_length, theta, t_circ = tf.t_fus_due_to_moment(sigma, Mx, My)
    y_pos = np.reshape(y_pos, (1,len(boom_area)))
    x_pos = np.reshape(x_pos, (1,len(boom_area)))
    
    I_xx = sum(boom_area*(y_pos[0] - y_centroid)**2)
    I_yy = sum(boom_area*(x_pos[0] - x_centroid)**2)
#    I_xy = sum(boom_area*(y_pos[0] - y_centroid)*(x_pos[0] - x_centroid))
#    print("ixx", I_xx)
    final_stress = np.zeros((len(theta), len(alpha_segment)))
#    part_mx = (Mx * np.cos(theta[i]) + My * np.sin(theta[i])) * R_2[j] * np.sin(alpha_segment[j]) / I_xx
#    part_my = (My * np.cos(theta[i]) - Mx * np.sin(theta[i])) * R_2[j] * np.cos(alpha_segment[j]) / I_yy

    for i in range(len(theta)):
        
        for j in range(len(alpha_segment)): 
            
            final_stress[i][j] = (Mx * np.cos(theta[i]) + My * np.sin(theta[i])) * (y_pos[0][j] - y_centroid) / I_xx + ( My * np.cos(theta[i]) - Mx * np.sin(theta[i]))*(x_pos[0][j] - x_centroid)/I_yy


    return I_xx, I_yy, final_stress


boom_area_actual, t_iteratrion = fuselage_booms(tf.sigma, tf.Mx, tf.My)
print("req boom area", boom_area_actual*10000)

y_centroid, x_centroid = fuselage_centroid(boom_area_actual, y_pos, x_pos)

alpha_segment = tf.alpha_segment
theta = tf.theta
print(y_centroid, x_centroid)

plt.figure()
plt.plot(alpha_segment * 180 / np.pi, boom_area_actual * 10000)
plt.xlabel("\u03B1 [degrees]")
plt.ylabel("Areas per segment [mm]")
plt.legend()
plt.show()

I_xx, I_yy, sigma_final = fuselage_inertia(boom_area_actual, y_centroid, x_centroid, y_pos, x_pos, tf.My, tf.Mx, tf.sigma)
plt.figure()
plt.plot(alpha_segment * 180 / np.pi, sigma_final[0], label = "theta ="+ str(round(theta[0]*180/np.pi)))
plt.plot(alpha_segment * 180 / np.pi, sigma_final[1], label = "theta ="+str(theta[1]*180/np.pi))
plt.plot(alpha_segment * 180 / np.pi, sigma_final[2], label = "theta ="+str(theta[2]*180/np.pi))
plt.plot(alpha_segment * 180 / np.pi, sigma_final[3], label = "theta ="+str(theta[3]*180/np.pi))
plt.plot(alpha_segment * 180 / np.pi, sigma_final[4], label = "theta ="+str(theta[4]*180/np.pi))
plt.plot(alpha_segment * 180 / np.pi, sigma_final[5], label = "theta ="+str(theta[5]*180/np.pi))
plt.plot(alpha_segment * 180 / np.pi, sigma_final[6], label = "theta ="+str(round(theta[6]*180/np.pi)))

plt.xlabel("\u03B1 [degrees]")
plt.ylabel("Stress distribution [N/m$^{2}$]")
plt.legend(bbox_to_anchor=(1.03,1), loc="upper left")    
plt.show()


print("Ixx", I_xx, "Iyy", I_yy)
print("stress", sigma_final)

#
def fuselage_stiffener_cost(density, stiffener_length, cost, boom_area):
    #option 1: having different boom areas
    stiffener_volume = sum(boom_area*stiffener_length)
    fuselage_mass = stiffener_volume*density
    fuselage_cost = fuselage_mass*cost
    #option 2: having the same boom areas but different pitch
    return fuselage_mass, fuselage_cost

stiffener_length = 49.
cost = 1.96
density = 2780

mass_fuselage, cost_fuselage = fuselage_stiffener_cost(density, stiffener_length, cost, boom_area_actual)

print("fus mass",mass_fuselage, "fus cost", cost_fuselage)
#    return
#    