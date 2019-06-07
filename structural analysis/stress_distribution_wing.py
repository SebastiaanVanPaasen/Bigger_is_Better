# -*- coding: utf-8 -*-
"""
Created on Wed May 15 09:51:57 2019

@author: Mathilde

"""
from matplotlib import pyplot as plt
import numpy as np  ### Never use * to import stuff, as it makes it difficult to retrace where functions come from
from airfoil_geometry import airfoil_geometry
import loading_and_moment_diagrams as lm
from loading_and_moment_diagrams import load_diagrams,c  #   return data_z_all_sec, data_y_upper_all_sec, data_y_lower_all_sec
#from loading_and_moment_diagrams import load_diagrams

N= 100
b = lm.b
HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)

My_wing = lm.load_diagrams(N)[1]
Mz_wing = lm.load_diagrams(N)[0]
T_wing = lm.load_diagrams(N)[2]

I_yy_wing =  0.1
I_zz_wing = 0.1
I_zy_wing = 0.01

def wing_stress(My_wing, Mz_wing, T_wing, HalfspanValues, I_yy_wing, I_zz_wing, I_zy_wing):

    z_nodes = airfoil_geometry(N,b)[0]
    y_up_nodes = airfoil_geometry(N,b)[1]
    y_low_nodes = airfoil_geometry(N,b)[2]

    
    local_stress_up = np.zeros((len(HalfspanValues), len(z_nodes[0])))
    local_stress_low = np.zeros((len(HalfspanValues), len(z_nodes[0])))
    dist_qc = np.zeros((len(HalfspanValues), len(z_nodes[0])))

    for i in range(len(HalfspanValues)):
        for j in range(len(z_nodes[0])):
            dist_qc[i][j] = z_nodes[i][j] - 0.25*lm.c(HalfspanValues[i])
            local_stress_up[i][j] = ((-My_wing[i]*I_zz_wing - Mz_wing[i]*I_zy_wing)*dist_qc[i][j] + (-Mz_wing[i]*I_yy_wing - My_wing[i]*I_zy_wing)*y_up_nodes[i][j])/(I_zz_wing*I_yy_wing - I_zy_wing**2)
            local_stress_low[i][j] = ((-My_wing[i]*I_zz_wing - Mz_wing[i]*I_zy_wing)*dist_qc[i][j] + (-Mz_wing[i]*I_yy_wing - My_wing[i]*I_zy_wing)*y_low_nodes[i][j])/(I_zz_wing*I_yy_wing - I_zy_wing**2)

    return dist_qc, local_stress_up[0],local_stress_low[0]

dist_qc, stress_up, stress_low = wing_stress(My_wing, Mz_wing, T_wing, HalfspanValues, I_yy_wing, I_zz_wing, I_zy_wing)
#plt.plot(dist_qc[0], stress_up)
#plt.plot(z_nodes[0], stress_low)
#plt.show()

R = 2.5
fus_sec = list(np.arange(0,31,1))
I_xx_fus = 2
I_yy_fus = 2
I_xy_fus = 1
My_fus = np.arange(30000000,33100000,100000)
Mx_fus = np.arange(10000000,13100000,100000)

def fuselage_stress(R, fus_sec, My_fus, Mx_fus,I_xx_fus, I_yy_fus, I_xy_fus):
    
    alpha_deg = list(np.arange(0,365,5)*np.pi/180)#first approximation of the fuselage Izz, assuming circular with constant thickness 
    x_pos = []
    y_pos =[]
    for k in alpha_deg:
        x_pos.append(R*np.cos(k))
        y_pos.append(R*np.sin(k))

    local_stress = np.zeros((len(fus_sec),len(alpha_deg)))
    for i in range(len(fus_sec)):
        for j in range(len(alpha_deg)):
            local_stress[i][j] = ((-My_fus[i]*I_xx_fus - Mx_fus[i]*I_xy_fus)*x_pos[j] + (-Mx_fus[i]*I_yy_fus - My_fus[i]*I_xy_fus)*y_pos[j])/(I_xx_fus*I_yy_fus - I_xy_fus**2)

    
    return x_pos, y_pos, local_stress

x_pos,y_pos,local_stress = fuselage_stress(R, fus_sec, My_fus, Mx_fus,I_xx_fus, I_yy_fus, I_xy_fus)
#plt.plot(x_pos[0:35],local_stress[0][0:35])
#plt.show()