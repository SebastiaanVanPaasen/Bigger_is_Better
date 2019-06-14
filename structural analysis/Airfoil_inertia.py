# -*- coding: utf-8 -*-
"""
Created on Wed May 15 16:01:47 2019

@author: floyd
"""
import numpy as np
from airfoil_geometry import airfoil_geometry
import centroid_wing as cw

#from stress_distribution_wing import wing_stress

#N = 100 
#b = lm.b#60.#47.83#39.56#41.76

def s_airfoil(N,b,c):
    
    data_z_all_sec = airfoil_geometry(N,b,c)[0]
    data_y_upper_all_sec = airfoil_geometry(N,b,c)[1]
    data_y_lower_all_sec = airfoil_geometry(N,b,c)[2]
    ds_sec_all = []
    s_all_sec = []
    
    for i in range(len(data_z_all_sec)):
        for j in range(len(data_z_all_sec[0])-1):
            ds_sec = np.sqrt((data_z_all_sec[i][j+1]-data_z_all_sec[i][j])**2+(data_y_upper_all_sec[i][j+1]-data_y_upper_all_sec[i][j])**2) + np.sqrt((data_z_all_sec[i][j+1]-data_z_all_sec[i][j])**2+(data_y_lower_all_sec[i][j+1]-data_y_lower_all_sec[i][j])**2)
            ds_sec_all.append(ds_sec)
    
    ds_sec_all = np.asarray(ds_sec_all)
    ds_sec_all = np.reshape(ds_sec_all, (len(data_z_all_sec),len(data_z_all_sec[0])-1))
    
    for i in range(len(ds_sec_all)):
        s_all_sec.append(sum(ds_sec_all[i]))
  
    return(s_all_sec)


    
    
def I_zz_spars(l_spar_h ,t_spar_v, t_spar_h, N, b, c):
    
    HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
    airfoil_area, z_c_airfoil, y_c_airfoil = cw.get_skin_centroid(b, N, c)
    
    z_centroid_all_sec, y_centroid_all_sec, y_loc_spar_up, y_loc_spar_low, y_loc_stiff_up, y_loc_stiff_low, y_vertical_spar, z_loc_stiff_up, spar_loc_sec, z_loc_stiff_low, spar_areas_verti = cw.wing_centroid(cw.boom_area, cw.spar_areas_hori, cw.t_spar_v, z_c_airfoil, y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, N, b, c)


    I_zz_spar = np.zeros((len(HalfspanValues), 1))
    I_yy_spar = np.zeros((len(HalfspanValues), 1))
    I_yz_spar = np.zeros((len(HalfspanValues), 1))
    

    
    for j in range(len(HalfspanValues)):
        
        
        I_zz_up=0
        I_yy_up=0
        I_yz_up=0
        I_zz_low=0
        I_yy_low=0
        I_yz_low=0
        I_zz_vertical=0
        I_yy_vertical=0
        I_yz_vertical =0
    
        
        for i in range(len(cw.spar_areas_hori)):
            I_zz_up += (1/12) * l_spar_h*t_spar_h**3 + l_spar_h*t_spar_h*(-y_loc_spar_up[j][i]-(-1)*y_centroid_all_sec[j])**2  
            I_yy_up += (1/12) * t_spar_h*l_spar_h**3 + l_spar_h*t_spar_h*(spar_loc_sec[j][i]-z_centroid_all_sec[j])**2
            I_yz_up += l_spar_h*t_spar_h*(-y_loc_spar_up[j][i]-(-1)*y_centroid_all_sec[j])*(spar_loc_sec[j][i]-z_centroid_all_sec[j])
            
            I_zz_low += (1/12) * l_spar_h*t_spar_h**3 + l_spar_h*t_spar_h*(-y_loc_spar_low[j][i] - (-1)*y_centroid_all_sec[j])**2
            I_yy_low += (1/12) * t_spar_h*l_spar_h**3 + l_spar_h*t_spar_h*(spar_loc_sec[j][i]-z_centroid_all_sec[j])**2
            I_yz_low += l_spar_h*t_spar_h*(-y_loc_spar_low[j][i] - (-1)*y_centroid_all_sec[j])*(spar_loc_sec[j][i]-z_centroid_all_sec[j])
            
            l_spar_v = y_loc_spar_up[j][i] - y_loc_spar_low[j][i]
#            print("l_spar_v", l_spar_v)
            I_zz_vertical += (1/12) * t_spar_v* l_spar_v**3 + l_spar_v*t_spar_v*(-y_vertical_spar[j][i] - (-1)*y_centroid_all_sec[j])**2
            I_yy_vertical += (1/12) * l_spar_v* t_spar_v**3 + l_spar_v*t_spar_v*(spar_loc_sec[j][i] - z_centroid_all_sec[j])**2
            I_yz_vertical += l_spar_v* t_spar_v*(-y_vertical_spar[j][i] - (-1)*y_centroid_all_sec[j])*(spar_loc_sec[j][i] - z_centroid_all_sec[j])
            
        I_zz = I_zz_up + I_zz_low + I_zz_vertical
        I_yy = I_yy_up + I_yy_low + I_yy_vertical
        I_yz = I_yz_up + I_yz_low + I_yz_vertical
            
        I_zz_spar[j] = I_zz
        I_yy_spar[j] = I_yy
        I_yz_spar[j] = I_yz

    return I_zz_spar, I_yy_spar, I_yz_spar

#print(I_zz_spars(cw.l_spar_h, cw.t_spar_v, cw.t_spar_h)[0][0])
##print(I_zz_spars(cw.l_spar_h, cw.t_spar_v, cw.t_spar_h)[1][0])
##print(I_zz_spars(cw.l_spar_h, cw.t_spar_v, cw.t_spar_h)[2][0])
#
#I_zz_spars = I_zz_spars(cw.l_spar_h, cw.t_spar_v, cw.t_spar_h)
#I_zz_req = required_Izz(HalfspanValues, data_y_all_sec, y_centroid_all_sec, Mz)

#first define the wing lay out that is required to obtain the right moment of inertia
def wing_geometry(I_zz_req, I_zz_spars, N, b, c):
    
    HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
    airfoil_area, z_c_airfoil, y_c_airfoil = cw.get_skin_centroid(N, b, c)
    z_centroid_all_sec, y_centroid_all_sec, y_loc_spar_up, y_loc_spar_low, y_loc_stiff_up, y_loc_stiff_low, y_vertical_spar, z_loc_stiff_up, spar_loc_sec, z_loc_stiff_low, spar_areas_verti = cw.wing_centroid(cw.boom_area, cw.spar_areas_hori, cw.t_spar_v, z_c_airfoil, y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, N, b, c)

    single_boom_area = np.zeros((len(HalfspanValues), 1))
    
    for i in range(len(HalfspanValues)):
        
        y_2 = 0
        
        for k in range(len(y_loc_stiff_up[0])):
            y_2 += (-y_loc_stiff_up[i][k] - (-1)*y_centroid_all_sec[i])**2 
    
        for j in range(len(y_loc_stiff_low[0])):
            y_2 += (- y_loc_stiff_low[i][j] - (-1)*y_centroid_all_sec[i])**2 
        
        single_boom_area[i][0] = (I_zz_req[i][0] - I_zz_spars[i][0])/y_2
        
#        print(y_2*single_boom_area)
#    plt.plot(HalfspanValues, single_boom_area)
#    plt.show()

    return single_boom_area

#print("single boom area", wing_geometry(I_zz_req, I_zz_spars, y_loc_stiff_up, y_loc_stiff_low)[0]) 
#I_zz_spar = I_zz_spars[0]
#I_yy_spar = I_zz_spars[1]
#I_yz_spar = I_zz_spars[2]
#airfoil_area = cw.airfoil_area
#z_c_airfoil = cw.z_c_airfoil
#y_c_airfoil = cw.y_c_airfoil
#boom_area = wing_geometry(I_zz_req, I_zz_spars, y_loc_stiff_up, y_loc_stiff_low)[1][0]
#calculate the wing moment of inertia per section 
    

def inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, boom_area, N, b, c):
    HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
    airfoil_area, z_c_airfoil, y_c_airfoil = cw.get_skin_centroid(N, b, c)

    z_centroid_all_sec, y_centroid_all_sec, y_loc_spar_up, y_loc_spar_low, y_loc_stiff_up, y_loc_stiff_low, y_vertical_spar, z_loc_stiff_up, spar_loc_sec, z_loc_stiff_low, spar_areas_verti = cw.wing_centroid(cw.boom_area, cw.spar_areas_hori, cw.t_spar_v, z_c_airfoil, y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, N, b, c)

    I_zz = np.zeros((len(HalfspanValues), 1))
    I_yy = np.zeros((len(HalfspanValues), 1))
    I_yz = np.zeros((len(HalfspanValues), 1))
    
    
    for i in range(len(HalfspanValues)):
        
        I_zz_booms=0
        I_yy_booms=0
        I_yz_booms=0
        I_zz_airfoil = 0
        I_yy_airfoil = 0
        I_yz_airfoil = 0
        
        for j in range(len(y_loc_stiff_up[0])):
            
            I_zz_booms += boom_area*(-y_loc_stiff_up[i][j] - (-1)*y_centroid_all_sec[i])**2
            I_yy_booms += boom_area*(z_loc_stiff_up[i][j] - z_centroid_all_sec[i])**2
            I_yz_booms += boom_area*(-y_loc_stiff_up[i][j] - (-1)*y_centroid_all_sec[i])*(z_loc_stiff_up[i][j] - z_centroid_all_sec[i])
        
        for k in range(len(y_loc_stiff_low[0])):
            
            I_zz_booms += boom_area*(-y_loc_stiff_low[i][k] - (-1)*y_centroid_all_sec[i])**2
            I_yy_booms += boom_area*(z_loc_stiff_low[i][k] - z_centroid_all_sec[i])**2
            I_yz_booms += boom_area*(-y_loc_stiff_low[i][k] - (-1)*y_centroid_all_sec[i])*(z_loc_stiff_low[i][k] - z_centroid_all_sec[i])
        
        
        I_zz_airfoil += airfoil_area[i]*(-y_c_airfoil[i] - (-1)*y_centroid_all_sec[i])**2
        I_yy_airfoil += airfoil_area[i]*(z_c_airfoil[i] - z_centroid_all_sec[i])**2
        I_yz_airfoil += airfoil_area[i]*(-y_c_airfoil[i] - (-1)*y_centroid_all_sec[i])*(z_c_airfoil[i] - z_centroid_all_sec[i])

        I_zz[i][0] = I_zz_booms + I_zz_spar[i][0] + I_zz_airfoil
        I_yy[i][0] = I_yy_booms + I_yy_spar[i][0] + I_yy_airfoil
        I_yz[i][0] = I_yz_booms + I_yz_spar[i][0] + I_yz_airfoil
    
#        print(I_zz_booms)
#        print(I_yz_booms)
#        
#    print(I_zz_spars[0][0])

    return I_zz, I_yy, I_yz
    
#print(inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, airfoil_area, z_c_airfoil, y_c_airfoil, boom_area)[0][0])
#print(inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, airfoil_area, z_c_airfoil, y_c_airfoil, boom_area)[1][0])
#print(inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, airfoil_area, z_c_airfoil, y_c_airfoil, boom_area)[2][0])


