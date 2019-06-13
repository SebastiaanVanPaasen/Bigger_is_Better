# -*- coding: utf-8 -*-
"""
Created on Tue May 28 09:19:12 2019

@author: mathi
"""

#import loading_and_moment_diagrams as lm
from loading_and_moment_diagrams import load_diagrams
from centroid_wing import wing_centroid
import centroid_wing as cw
import numpy as np

z_centroid_all_sec, y_centroid_all_sec, y_loc_spar_up, y_loc_spar_low, y_loc_stiff_up, y_loc_stiff_low, y_vertical_spar, z_loc_stiff_up, spar_loc_sec,z_loc_stiff_low = wing_centroid(cw.boom_area, cw.spar_areas_hori, cw.t_spar_v, cw.z_c_airfoil, cw.y_c_airfoil, cw.n_stiff_up, cw.n_stiff_low, cw.HalfspanValues)
data_y_all_sec_up = cw.data_y_upper_all_sec
data_y_all_sec_low = cw.data_y_lower_all_sec
Mz = load_diagrams(100)[0]

def required_Izz(HalfspanValues, data_y_all_sec_up, data_y_all_sec_low, y_centroid_all_sec, Mz):
    
    I_zz_req_all_sec = np.zeros((len(HalfspanValues), 1))
    sigma_ult = 110 * 10 ** 6
    
    for i in range(len(HalfspanValues)):
        y_up_max = -max(data_y_all_sec_up[i]) - (-1)*y_centroid_all_sec[i]
        y_low_max = -min(data_y_all_sec_low[i]) - (-1)*y_centroid_all_sec[i]
        print(y_low_max)
        y_max = max(abs(y_up_max),abs(y_low_max))
        I_zz_req_all_sec[i] = abs(Mz[i]) * y_max / sigma_ult
#        print(y_max)
    return I_zz_req_all_sec

#8.54#7.11#6.63#6.06
print("req Izz",required_Izz(cw.HalfspanValues, data_y_all_sec_up, data_y_all_sec_low, y_centroid_all_sec, Mz))
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