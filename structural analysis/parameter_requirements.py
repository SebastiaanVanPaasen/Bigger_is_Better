
#import loading_and_moment_diagrams as lm

import centroid_wing as cw
import numpy as np
from airfoil_geometry import airfoil_geometry

#Mz = load_diagrams(100)[0]

def required_Izz(N, b, c, Mz):
    
    HalfspanValues = np.linspace(0, b / 2 - 0.00001, N)
    data_y_all_sec_up = airfoil_geometry(N, b, c)[1]
    data_y_all_sec_low = airfoil_geometry(N, b, c)[2]

    airfoil_area, z_c_airfoil, y_c_airfoil = cw.get_skin_centroid(N, b, c)
    z_centroid_all_sec, y_centroid_all_sec, y_loc_spar_up, y_loc_spar_low, y_loc_stiff_up, y_loc_stiff_low, y_vertical_spar, z_loc_stiff_up, spar_loc_sec,z_loc_stiff_low = cw.wing_centroid(cw.boom_area, cw.spar_areas_hori, cw.t_spar_v, z_c_airfoil, y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, N, b, c)

    I_zz_req_all_sec = np.zeros((len(HalfspanValues), 1))
    sigma_ult = 110 * 10 ** 6
    
    for i in range(len(HalfspanValues)):
        y_up_max = -max(data_y_all_sec_up[i]) - (-1)*y_centroid_all_sec[i]
        y_low_max = -min(data_y_all_sec_low[i]) - (-1)*y_centroid_all_sec[i]
#        print(y_low_max)
        y_max = max(abs(y_up_max),abs(y_low_max))
        I_zz_req_all_sec[i] = abs(Mz[i]) * y_max / sigma_ult
#        print(y_max)
    return I_zz_req_all_sec

#<<<<<<< HEAD


#print("req Izz",required_Izz(5.8451))

#strut_length = 20.95

#def required_strut_area(strut_length):
#    angle = (20.78/180)*np.pi
#    sigma_carbon = 552 * 10 ** 6#1500*10**6 
#    density_carbon = 2801#1600
#    cost_carbon = 6.6 #100
#    P_strut = lm.strutforce/np.sin(angle)
#    A_req = P_strut/sigma_carbon
#    strut_volume = A_req*strut_length
#    strut_mass = strut_volume*density_carbon
#    strut_cost = strut_mass*cost_carbon
#    print(lm.strutforce)
#    print(P_strut)
#    return 'A_req=',A_req, 'strut_mass=',strut_mass, 'strut_cost=',strut_cost
#
#print(required_strut_area(strut_length))