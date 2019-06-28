import numpy as np
import centroid_wing as cw
from airfoil_geometry import airfoil_geometry

#N = 100 
#b = lm.b#60.#47.83#39.56#41.76

def s_airfoil(N,b,c, X_root):
    
    data_z_all_sec = airfoil_geometry(N,b,c, X_root)[0]
    data_y_upper_all_sec = airfoil_geometry(N,b,c, X_root)[1]
    data_y_lower_all_sec = airfoil_geometry(N,b,c, X_root)[2]

    
    ds_sec_all = []
    s_all_sec = []
    
    for i in range(len(data_z_all_sec)):
        ds_sec = []
        
        for j in range(len(data_z_all_sec[0]) - 1):
            ds_sec_up = np.sqrt((data_z_all_sec[i][j + 1] - data_z_all_sec[i][j]) ** 2 + (data_y_upper_all_sec[i][j + 1] - data_y_upper_all_sec[i][j]) ** 2)
            ds_sec.append(ds_sec_up)
            
        ds_sec.append(data_y_upper_all_sec[i][-1] - data_y_lower_all_sec[i][-1])
    
        for j in range(len(data_z_all_sec[0]) - 1):
            ds_sec_low = np.sqrt((data_z_all_sec[i][::-1][j + 1] - data_z_all_sec[i][::-1][j]) ** 2 + (data_y_lower_all_sec[i][::-1][j + 1] - data_y_lower_all_sec[i][::-1][j]) ** 2)
            ds_sec.append(ds_sec_low)
                        
        ds_sec_all.append(ds_sec)
        s_all_sec.append(sum(ds_sec))

  
    return s_all_sec, ds_sec_all


    
    
def I_zz_spars(l_spar_h ,t_spar_v, t_spar_h, N, b, c, boom_area, X_root, dx):
    
    airfoil_area, z_c_airfoil, y_c_airfoil = cw.get_skin_centroid(N, b, c, dx)
    
    z_centroid_all_sec, y_centroid_all_sec, y_loc_spar_up, y_loc_spar_low, y_loc_stiff_up, y_loc_stiff_low, y_vertical_spar, z_loc_stiff_up, spar_loc_sec, z_loc_stiff_low, spar_areas_verti = cw.wing_centroid(boom_area, cw.spar_areas_hori, cw.t_spar_v, z_c_airfoil, y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, N, b, c, X_root, dx)


    I_zz_spar = np.zeros((len(X_root), 1))
    I_yy_spar = np.zeros((len(X_root), 1))
    I_yz_spar = np.zeros((len(X_root), 1))
    

    
    for j in range(len(X_root)):
        
        
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
        
#    print("I_zz in spar",I_zz_spar)
    return I_zz_spar, I_yy_spar, I_yz_spar

#first define the wing lay out that is required to obtain the right moment of inertia
def wing_geometry(I_zz_req, I_zz_spars, N, b, c, boom_area, X_root, dx):
#    dx = 0.1
#    print(X_root)
    airfoil_area, z_c_airfoil, y_c_airfoil = cw.get_skin_centroid(N, b, c, dx)
    z_centroid_all_sec, y_centroid_all_sec, y_loc_spar_up, y_loc_spar_low, y_loc_stiff_up, y_loc_stiff_low, y_vertical_spar, z_loc_stiff_up, spar_loc_sec, z_loc_stiff_low, spar_areas_verti = cw.wing_centroid(boom_area, cw.spar_areas_hori, cw.t_spar_v, z_c_airfoil, y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, N, b, c , X_root, dx)

    single_boom_area = np.zeros((len(X_root), 1))
#    print(z_centroid_all_sec)
#    print(y_centroid_all_sec)
#    print(np.shape(I_zz_spars))

    for i in range(len(X_root)):
        
        y_2 = 0
        I_zz_airfoil = 0
        
        for k in range(len(y_loc_stiff_up[0])):
            y_2 += (-y_loc_stiff_up[i][k] + y_centroid_all_sec[i])**2 
    
        for j in range(len(y_loc_stiff_low[0])):
            y_2 += (- y_loc_stiff_low[i][j] + y_centroid_all_sec[i])**2 
        
            
        I_zz_airfoil += airfoil_area[i]*(-y_c_airfoil[i] - (-1)*y_centroid_all_sec[i])**2
#        I_yy_airfoil += airfoil_area[i]*(z_c_airfoil[i] - z_centroid_all_sec[i])**2
#        I_yz_airfoil += airfoil_area[i]*(-y_c_airfoil[i] - (-1)*y_centroid_all_sec[i])*(z_c_airfoil[i] - z_centroid_all_sec[i])

#        single_boom_area[i] = (I_zz_req[i] - I_zz_spars[i][0] - I_zz_airfoil)/y_2
        
        
        single_boom_area[i] = 0.0050
#    plt.plot(X_root, single_boom_area)
#    plt.show()
#    print("single_boom_area",single_boom_area*10000, len(single_boom_area))
    return max(single_boom_area)



def inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, boom_area, N, b, c, X_root, dx):
    airfoil_area, z_c_airfoil, y_c_airfoil = cw.get_skin_centroid(N, b, c,dx)

    z_centroid_all_sec, y_centroid_all_sec, y_loc_spar_up, y_loc_spar_low, y_loc_stiff_up, y_loc_stiff_low, y_vertical_spar, z_loc_stiff_up, spar_loc_sec, z_loc_stiff_low, spar_areas_verti = cw.wing_centroid(boom_area, cw.spar_areas_hori, cw.t_spar_v, z_c_airfoil, y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, N, b, c, X_root, dx)
   
#    print("boom area in izz", boom_area)
    I_zz = np.zeros(len(X_root))
    I_yy = np.zeros(len(X_root))
    I_yz = np.zeros(len(X_root))
    
    
    for i in range(len(X_root)):
        
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

#        print("boom_area", boom_area)
#        print("I_zz_booms",I_zz_booms)
        I_zz[i] = I_zz_booms + I_zz_spar[i][0] + I_zz_airfoil
        I_yy[i] = I_yy_booms + I_yy_spar[i][0] + I_yy_airfoil
        I_yz[i] = I_yz_booms + I_yz_spar[i][0] + I_yz_airfoil
    

#        print("I_zz used",I_zz_spar[i][0])
#        
#    print("Izz wing",I_zz_airfoil)
    print("I_zz", I_zz[0])
    
    return I_zz, I_yy, I_yz
    
#print(inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, airfoil_area, z_c_airfoil, y_c_airfoil, boom_area)[0][0])
#print(inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, airfoil_area, z_c_airfoil, y_c_airfoil, boom_area)[1][0])
#print(inertia_wing(I_zz_spar, I_yy_spar, I_yz_spar, airfoil_area, z_c_airfoil, y_c_airfoil, boom_area)[2][0])


