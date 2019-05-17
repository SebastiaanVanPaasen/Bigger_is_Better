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

def input_CL(n,S,V,rho,W):
    input_CL = 2.5*W/(0.5*rho*V**2*S)
    return input_CL

n = 2.5
S = 200
V = 250
rho = 0.5
c = 8.
W = 2000000.
CD0 = 0.02
b=64.
output_avl = lift_distribution(input_CL(n,S,V,rho,W))

def CN_CT(output_avl,c, CD0,b):
    
    x_total,cl_total,cd_total = get_correct_data(output_avl,c)[0], get_correct_data(output_avl,c)[1], get_correct_data(output_avl,c)[2] 
    dx = []
    aoa = []
    CL = []
    CD = []
    CN = []
    CT = []
    angle = (30/180) *np.pi #moet input van AVL worden 
    
    x_pos = x_total[0:int((len(x_total)/2))] 
    cl = cl_total[0:int((len(x_total)/2))]
    cd = cd_total[0:int((len(x_total)/2))]
    
    for i in range(len(x_pos)-1):
        diffx = abs(x_pos[i+1] - x_pos[i])
        C_L = (cl[i+1] + cl[i])/2.
        C_D = ((cd[i+1] + cd[i])/2.) + CD0
        CL.append(C_L)
        CD.append(C_D)
        dx.append(diffx)
        
    for i in range(len(CL)):
        C_N = CL[i]*np.cos(angle) - CD[i]*np.sin(angle)
        C_T = -CL[i]*np.sin(angle) + CD[i]*np.cos(angle)
        CN.append(C_N)
        CT.append(C_T)
    
    return CN,CT,dx,x_pos,CD


CN = CN_CT(output_avl,c, CD0,b)[0]
CT = CN_CT(output_avl,c, CD0,b)[1]
dx = CN_CT(output_avl,c, CD0,b)[2]
x_pos = CN_CT(output_avl,c, CD0,b)[3]
CD = CN_CT(output_avl,c, CD0,b)[4]
print(x_pos)
x_eng = 6.
F_eng = 70000.
angle= 0.523
x_fuelstart = 0.
x_fuelend = 5
#plt.plot(x_pos[:-1],CT)
#plt.show()

def fuel_surface():
    S_begin = 10
    S_end = 6
    S_slope = S_end/S_begin
    fuel_surface = np.zeros(len(dx))
    for i in range(len(dx)):
        if i ==0:
            
            fuel_surface[i] = S_begin 
        else: 
            fuel_surface[i]= fuel_surface[i-1]+ S_slope*dx[i]
            
    return fuel_surface

def external_loads(CN,CT,c,rho,dx,V,x_pos,x_eng,F_eng,angle,x_fuelstart,x_fuelend):
    Fy_right = []
    Fz_right = []
    x_wing = []
<<<<<<< Updated upstream
    fuel_density = 0.8
    g = 9.81
    
    for i in range(len(CN)):
    
       if x_pos[i]>x_fuelstart and x_pos[i]<x_fuelend:
           F_fuel = dx[i]*fuel_surface()[i]*fuel_density*g
           x_fuel = x_pos[i] + 0.5*dx[i]
           F_y = -F_fuel*np.cos(angle) + 0.5*rho*dx[i]*c*CN[i]*V**2
           F_z = F_fuel*np.sin(angle) + 0.5*rho*dx[i]*c*CT[i]*V**2
=======
    for i in range(len(CN)):
        
       x_r_wing = x_pos[i]+dx[i]
       
       if x_pos[i]<x_eng and x_pos[i+1]>x_eng:
           F_y = -F_eng*np.cos(angle)
           F_z = F_eng*np.sin(angle)
>>>>>>> Stashed changes
           Fy_right.append(F_y)
           Fz_right.append(F_z)
           x_wing.append(x_fuel)
           print("fuel=",x_fuel)
           
       if x_pos[i]<x_eng and x_pos[i+1]>x_eng:
           
           if  x_pos[i]>x_fuelstart and x_pos[i]<x_fuelend:
               F_fuel = dx[i]*fuel_surface()[i]*fuel_density*g
               x_fuel = x_pos[i] + 0.5*dx[i]
               F_y = -F_eng*np.cos(angle) -F_fuel*np.cos(angle) + 0.5*rho*dx[i]*c*CN[i]*V**2
               F_z = F_eng*np.sin(angle)+ F_fuel*np.sin(angle) + 0.5*rho*dx[i]*c*CT[i]*V**2
               Fy_right.append(F_y)
               Fz_right.append(F_z)
               x_wing.append(x_eng)
               print("hoi",x_eng)
               
           else:
               F_y = -F_eng*np.cos(angle) + 0.5*rho*dx[i]*c*CN[i]*V**2
               F_z = F_eng*np.sin(angle) + 0.5*rho*dx[i]*c*CT[i]*V**2
               Fy_right.append(F_y)
               Fz_right.append(F_z)
               x_wing.append(x_eng)
               print(x_eng)
           
       elif x_pos[i]<x_fuelstart or x_pos[i]>x_fuelend and x_pos[i]<x_eng or x_pos[i]>x_eng: 
           F_y = 0.5*rho*dx[i]*c*CN[i]*V**2
           F_z = 0.5*rho*dx[i]*c*CT[i]*V**2
           Fy_right.append(F_y)
           Fz_right.append(F_z)
           x_r_wing = x_pos[i]+ 0.5*dx[i]
           x_wing.append(x_r_wing)
           print(x_r_wing)
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
    
    
    print(x_l_wing)
    
    
    
    
    
    
    
    
    
    
    
#    plt.subplot(2,1,1)
#    plt.plot(x_total, Fy_tot)
#    plt.subplot(2,1,2)
#    plt.plot(x_wing, Vy)
#    plt.plot(x_total, Fz_tot)
#    plt.show()
    
    return(Fy_left,Fz_left,x_l_wing)


<<<<<<< Updated upstream
print(external_loads(CN,CT,c,rho,dx,V,x_pos,x_eng,F_eng, angle,x_fuelstart,x_fuelend))
=======
>>>>>>> Stashed changes

