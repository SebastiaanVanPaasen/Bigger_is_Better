# -*- coding: utf-8 -*-
"""
Created on Thu May 16 09:55:34 2019

@author: Mathilde
"""

import numpy as np  ### Never use * to import stuff, as it makes it difficult to retrace where functions come from
#import matplotlib.pyplot as plt
#from stress_distribution_wing import load_airfoil



def load_airfoil(filename):
    f = open(filename,'r')
    lines = f.readlines()
    data = []
    data_z = []
    data_y = []
    for line in lines:
        x = line.split()
        data.append(x)
    for i in range(len(data)):
        for j in range(len(data[i])):
            data[i][j] = float(data[i][j])
    for i in range(len(data)):
        data_z.append(data[i][0])
        data_y.append(data[i][1])

#    plt.plot(data_z,data_y)
#    plt.axis('equal')
#    plt.show()
    return data,data_z,data_y
#
#filename = 'SC(2)-0616.txt'
#print(load_airfoil(filename))

def airfoil_geometry(N,b, c, X_root):

    
    data_z_all_sec = []
    data_y_lower_all_sec = []
    data_y_upper_all_sec = []
    
    data_z1, data_y1 = load_airfoil('SC(2)-0616.txt')[1], np.array(load_airfoil('SC(2)-0616.txt')[2]) 
    data_z2, data_y2 = load_airfoil('SC(2)-0612.txt')[1], np.array(load_airfoil('SC(2)-0612.txt')[2])

    tc_1 = 0.16
    tc_2 = 0.12
    
    for i in range(len(X_root)):
             
        if X_root[i]<(b/8):
            
            data_z_order =  np.array(data_z1[0:int((len(data_y1)/2))])*c(X_root[i])
            data_z_all_sec.append(data_z_order)
   
            data_y_lower = np.array(data_y1[(int((len(data_y1)/2))):])*c(X_root[i])*(tc_1 - (X_root[i]/(b/8))*(tc_1-tc_2))/tc_1
            data_y_lower_all_sec.append(data_y_lower)
            data_y_upper = np.array(data_y1[0:int((len(data_y1)/2))])*c(X_root[i])*(tc_1 - (X_root[i]/(b/8))*(tc_1-tc_2))/tc_1
            data_y_upper_all_sec.append(data_y_upper)
        else:
            data_z_order =  np.array(data_z2[0:int((len(data_y2)/2))])*c(X_root[i])
            data_z_all_sec.append(data_z_order)
   
            data_y_lower = np.array(data_y2[(int((len(data_y2)/2))):])*c(X_root[i])
            data_y_lower_all_sec.append(data_y_lower)
            data_y_upper = np.array(data_y2[0:int((len(data_y2)/2))])*c(X_root[i])
            data_y_upper_all_sec.append(data_y_upper)
            
            
        
    data_z_all_sec = np.asarray(data_z_all_sec)
    data_z_all_sec = np.reshape(data_z_all_sec, (len(X_root),len(data_z_order)))
    data_y_upper_all_sec = np.asarray(data_y_upper_all_sec)
    data_y_upper_all_sec = np.reshape(data_y_upper_all_sec, (len(X_root),len(data_z_order)))
    data_y_lower_all_sec = np.asarray(data_y_lower_all_sec)
    data_y_lower_all_sec = np.reshape(data_y_lower_all_sec,(len(X_root),len(data_z_order)))
    
#    print(data_y_upper_all_sec[0])
#    print(data_y_lower_all_sec[10])
    
    return data_z_all_sec, data_y_upper_all_sec, data_y_lower_all_sec

#dx = 0.1
#b = 56.3
#X_root = np.arange(0, (b/2)+dx, dx)
#
#cr = 5.53
#ct = 1.97
#L_wing = b/2
#
#def calc_chord(x):
#    return cr - ((cr - ct) / L_wing) * x
#
#data_z_all_sec, data_y_upper_all_sec, data_y_lower_all_sec = airfoil_geometry(100,b,calc_chord,X_root)
#
#
#
##data = load_airfoil('NACA3414.txt')
##for i in range(int((b/8)/dx)+1):
###    print(data_y_upper_all_sec[i])#int((b/8)/dx)+1])
##    if i/10 ==0:
#plt.axis('equal')
#plt.plot(data_z_all_sec[0], data_y_upper_all_sec[0])
#plt.show()

#plt.scatter(data_z_all_sec[0], data_y_upper_all_sec[0])
#plt.scatter(data_z_all_sec[0], data_y_lower_all_sec[0])
#plt.show()
##print(data_z_all_sec[0])
#print(data_y_upper_all_sec[0])
##print(data_y_lower_all_sec[0])
#print(np.shape(data_z_all_sec))