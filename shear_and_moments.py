# -*- coding: utf-8 -*-
"""
Created on Fri May 10 16:50:26 2019

@author: mathi
"""

import numpy as np
import subprocess
import os
from matplotlib import pyplot as plt
from math import *

from lift_distr import *

output_avl = lift_distribution(2.5)
mac = 8.

def CN_CT(output_avl,mac):
    
    x_total,cl_total,cd_total = get_correct_data(output_avl,mac)[0], get_correct_data(output_avl,mac)[1], get_correct_data(output_avl,mac)[2] 
    dx = []
    aoa = []
    CL = []
    CD = []
    CN = []
    CT = []
    angle = (30/180) *np.pi
    
    x_pos = x_total[0:int((len(x_total)/2))]
    cl = cl_total[0:int((len(x_total)/2))]
    cd = cd_total[0:int((len(x_total)/2))]
    
    for i in range(len(x_pos)-1):
        diffx = abs(x_pos[i+1] - x_pos[i])
        C_L = (cl[i+1] + cl[i])/2.
        C_D = (cd[i+1] + cd[i])/2.
        CL.append(C_L)
        CD.append(C_D)
        dx.append(diffx)
        
    for i in range(len(CL)):
        C_N = CL[i]*np.cos(angle) - CD[i]*np.sin(angle)
        C_T = -CL[i]*np.sin(angle) + CD[i]*np.cos(angle)
        CN.append(C_N)
        CT.append(C_T)
    
    return CN,CT,dx,x_pos


CN = CN_CT(output_avl,mac)[0]
CT = CN_CT(output_avl,mac)[1]
dx = CN_CT(output_avl,mac)[2]
x_pos = CN_CT(output_avl,mac)[3]
rho = 0.5
V = 230
c = 8.
x_eng = 4.
F_eng = 7000000.
angle= 0.523
#plt.plot(x_pos[:-1],CT)
#plt.show()

def external_loads(CN,CT,c,rho,dx,V,x_pos,x_eng,F_eng, angle):
    Fy_right = []
    Fz_right = []
    x_wing = []
    for i in range(len(CL)):
        
       x_r_wing = x_pos[i]+dx[i]
       
       if x_pos[i]<x_eng and x_pos[i+1]>x_eng:
           F_y = -F_eng*np.cos(angle)
           F_z = F_eng*np.sin(angle)
           Fy_right.append(F_y)
           Fz_right.append(F_z)
           x_wing.append(x_eng)
       else: 
           F_y = 0.5*rho*dx[i]*c*CN[i]*V**2
           F_z = 0.5*rho*dx[i]*c*CT[i]*V**2
           Fy_right.append(F_y)
           Fz_right.append(F_z)
           x_wing.append(x_r_wing)
    #plot the external loading on single side of the wing

    #calculating the internal loading 
    Vy_root = -sum(Fy_right)
    Vz_root = -sum(Fz_right)
    Vy = np.zeros(len(x_wing))
    for i in range(len(x_wing)):
        
        if i==0:    
            Vy[i]= Vy_root + Fy_right[i]
        else:
            Vy[i] = Vy[i-1] + Fy_right[i]

    x_l_wing = list(np.array(x_wing[::-1])*-1)
    Fy_left = Fy_right[::-1]
    Fz_left = Fz_right[::-1]
    x_total = x_l_wing + x_wing
    Fy_tot = Fy_left + Fy_right
    Fz_tot = Fz_left + Fz_right
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    plt.subplot(2,1,1)
    plt.plot(x_total, Fy_tot)
    plt.subplot(2,1,2)
#    plt.plot(x_wing, Vy)
    plt.plot(x_total, Fz_tot)
    plt.show()
    
    return

print(external_loads(CN,CT,c,rho,dx,V,x_pos,x_eng,F_eng, angle))

