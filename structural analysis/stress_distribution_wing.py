from matplotlib import pyplot as plt
import numpy as np  
from airfoil_geometry import airfoil_geometry
import parameter_requirements as pr
import practised_ultimate as prac
import centroid_wing as cw
import Airfoil_inertia as ai


def wing_stress(b, Mz, My, X_root):
    
#    l_spar_h, t_spar_v, t_spar_h = cw.l_spar_h, cw.t_spar_v, cw.t_spar_h
    N = prac.L_wing / prac.dx
    
#    I_zz_spar, I_yy_spar, I_yz_spar = ai.I_zz_spars(l_spar_h, t_spar_v, t_spar_h, N, b ,prac.calc_chord,prac.boom_area_all)
#    I_zz_req = pr.required_Izz(N, b, prac.calc_chord, Mz)
#    
    airfoil_area, z_c_airfoil, y_c_airfoil = cw.get_skin_centroid(N,b, prac.calc_chord, prac.dx)
    z_centroid_all_sec, y_centroid_all_sec, y_loc_spar_up, y_loc_spar_low, y_loc_stiff_up, y_loc_stiff_low, y_vertical_spar, z_loc_stiff_up, spar_loc_sec, z_loc_stiff_low, spar_areas_verti = cw.wing_centroid(prac.boom_area_all[0], cw.spar_areas_hori, cw.t_spar_v, z_c_airfoil, y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, N, b, prac.calc_chord, X_root, prac.dx)
#    print(z_centroid_all_sec)
#    boom_area = ai.wing_geometry(I_zz_req, I_zz_spar, N, b, prac.calc_chord)[0][0]
#    I_zz_wing, I_yy_wing, I_yz_wing = ai.inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, boom_area, N, b, prac.calc_chord)

    I_yy_wing = prac.I_yy_wing
    I_zz_wing = prac.I_zz_sections
    I_yz_wing = prac.I_yz_wing
#    I_yy_wing = inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, boom_area, N, b, c)
#    I_zz_wing = inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, boom_area, N, b, c)
#    I_zy_wing = inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, boom_area, N, b, c)
#    
#    z_centroid = wing_centroid(boom_area, spar_areas_hori, t_spar_v, z_c_airfoil, y_c_airfoil, n_stiff_up, n_stiff_low, N, b, c)
#    y_centroid = wing_centroid(boom_area, spar_areas_hori, t_spar_v, z_c_airfoil, y_c_airfoil, n_stiff_up, n_stiff_low, N, b, c)
#    
    Mz_wing = Mz
    My_wing = My
    
    z_nodes = airfoil_geometry(N,b,prac.calc_chord, X_root)[0] #adapt to centroid
    y_up_nodes = airfoil_geometry(N,b,prac.calc_chord, X_root)[1]  #adapt
    y_low_nodes = airfoil_geometry(N,b,prac.calc_chord, X_root)[2]

    
    local_stress_up = np.zeros((len(X_root), len(z_nodes[0])))
    local_stress_low = np.zeros((len(X_root), len(z_nodes[0])))
    z = np.zeros((len(X_root), len(z_nodes[0])))
    y_up = np.zeros((len(X_root), len(y_up_nodes[0])))
    y_low = np.zeros((len(X_root), len(y_low_nodes[0])))
#    print(z_nodes)
#    print(y_up_nodes[0])
    
    for i in range(len(X_root)):
        for j in range(len(z_nodes[0])):
            z[i][j] = z_nodes[i][j] - z_centroid_all_sec[i]
            y_up[i][j] = -(y_up_nodes[i][j]-y_centroid_all_sec[i])
            y_low[i][j] = (y_centroid_all_sec[i] - y_low_nodes[i][j])
            local_stress_up[i][j] = ((My_wing[i]*I_zz_wing[i] + Mz_wing[i]*I_yz_wing[i])*z[i][j] + (-Mz_wing[i]*I_yy_wing[i] - My_wing[i]*I_yz_wing[i])*y_up[i][j])/(I_zz_wing[i]*I_yy_wing[i] - I_yz_wing[i]**2)
            local_stress_low[i][j] = ((My_wing[i]*I_zz_wing[i] + Mz_wing[i]*I_yz_wing[i])*z[i][j] + (-Mz_wing[i]*I_yy_wing[i] - My_wing[i]*I_yz_wing[i])*y_low[i][j])/(I_zz_wing[i]*I_yy_wing[i] - I_yz_wing[i]**2)
            
#    print(I_zz_wing)
    
    return z, local_stress_up,local_stress_low


Mz = prac.Mz_dist
My = prac.My_dist

#print(wing_stress(60, Mz[0]))
max_stress_up = np.zeros((len(prac.A_S_L),len(prac.X_root))) 
max_stress_low = np.zeros((len(prac.A_S_L),len(prac.X_root))) 
min_stress_up = np.zeros((len(prac.A_S_L),len(prac.X_root)))
min_stress_low = np.zeros((len(prac.A_S_L),len(prac.X_root))) 

for i in range(len(prac.A_S_L)):
    
    z_pos, stress_up, stress_low = wing_stress(52, Mz[i], My[i], prac.X_root)
#    print(np.shape(stress_up))
    for j in range(len(prac.X_root)):
#        print(len(stress_up[j]))
        max_stress_up[i][j] = max(stress_up[j])
        max_stress_low[i][j] = max(stress_low[j])
        min_stress_up[i][j] = min(stress_up[j])
        min_stress_low[i][j] = min(stress_low[j])
        
   
#print(min(stress_up[0]))
#print(max(stress_low[0]))
#    plt.rcParams.update({'font.size': 20})        
##    plt.figure()
#    plt.plot(z_pos[i], stress_up[i], 'y', label = 'Initial stress top')
#    plt.plot(z_pos[i], stress_low[i], 'g', label = 'Initial stress bottom')
#    plt.xlabel("Chordwise position with respect to the centroid [m]")
#    plt.ylabel("Stress [N/m$^2$]")
#
#    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left")    
#    plt.show()
    
print("max_stress up",max(max_stress_up[0]))
print("max_stress low",max(max_stress_low[0]))
print("min_stress up",min(min_stress_up[0]))
print("min_stress low",min(min_stress_low[0]))
 
#R = 2.5
#fus_sec = list(np.arange(0,31,1))
#I_xx_fus = 2
#I_yy_fus = 2
#I_xy_fus = 1
#My_fus = np.arange(30000000,33100000,100000)
#Mx_fus = np.arange(10000000,13100000,100000)
#
#def fuselage_stress(R, fus_sec, My_fus, Mx_fus,I_xx_fus, I_yy_fus, I_xy_fus):
#    
#    alpha_deg = list(np.arange(0,365,5)*np.pi/180)#first approximation of the fuselage Izz, assuming circular with constant thickness 
#    x_pos = []
#    y_pos =[]
#    for k in alpha_deg:
#        x_pos.append(R*np.cos(k))
#        y_pos.append(R*np.sin(k))
#
#    local_stress = np.zeros((len(fus_sec),len(alpha_deg)))
#    for i in range(len(fus_sec)):
#        for j in range(len(alpha_deg)):
#            local_stress[i][j] = ((-My_fus[i]*I_xx_fus - Mx_fus[i]*I_xy_fus)*x_pos[j] + (-Mx_fus[i]*I_yy_fus - My_fus[i]*I_xy_fus)*y_pos[j])/(I_xx_fus*I_yy_fus - I_xy_fus**2)
#
#    
#    return x_pos, y_pos, local_stress
#
#x_pos,y_pos,local_stress = fuselage_stress(R, fus_sec, My_fus, Mx_fus,I_xx_fus, I_yy_fus, I_xy_fus)
##plt.plot(x_pos[0:35],local_stress[0][0:35])
#plt.show()