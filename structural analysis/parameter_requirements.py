
#import loading_and_moment_diagrams as lm

import centroid_wing as cw
import numpy as np
from airfoil_geometry import airfoil_geometry

#Mz = load_diagrams(100)[0]

def required_Izz(N, b, c, Mz, boom_area, X_root, dx):
    
    data_y_all_sec_up = airfoil_geometry(N, b, c, X_root)[1]
    data_y_all_sec_low = airfoil_geometry(N, b, c, X_root)[2]

    airfoil_area, z_c_airfoil, y_c_airfoil = cw.get_skin_centroid(N, b, c, dx)
    z_centroid_all_sec, y_centroid_all_sec, y_loc_spar_up, y_loc_spar_low, y_loc_stiff_up, y_loc_stiff_low, y_vertical_spar, z_loc_stiff_up, spar_loc_sec,z_loc_stiff_low, spar_areas_verti = cw.wing_centroid(boom_area, cw.spar_areas_hori, cw.t_spar_v, z_c_airfoil, y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, N, b, c, X_root,dx)

    I_zz_req_all_sec = np.zeros(len(X_root))
    sigma_ult = 114 * 10 ** 6
    
    for i in range(len(X_root)):
        y_up_max = -max(data_y_all_sec_up[i]) - (-1)*y_centroid_all_sec[i]
        y_low_max = -min(data_y_all_sec_low[i]) - (-1)*y_centroid_all_sec[i]
#        print(y_low_max)
        
#        print(Mz[i])
        
        y_max = max(abs(y_up_max),abs(y_low_max))
        I_zz_req_all_sec[i] = abs(Mz[i]) * y_max / sigma_ult
#        print("Mz",Mz[i])
#        print("y_max",y_max)
#    print("I_zz_req",I_zz_req_all_sec)

#        print(y_max)
        
    return I_zz_req_all_sec

