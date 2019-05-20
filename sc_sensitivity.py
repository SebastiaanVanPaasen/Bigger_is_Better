# -*- coding: utf-8 -*-
"""
Created on Thu May 16 08:51:34 2019

@author: Hidde
"""
import matplotlib.pyplot as plt
import numpy as np
from iterator_function import *
from class_I.loading_diagram import potato
from class_I.seats import *
from stability_and_control.stability import C_L_alpha_Ah, Sh_S_stability
from stability_and_control.control_curve import Sh_S_control
from stability_and_control.control_stability import control_stability_plot
#from input_files.aerodynamic_concept import *
#------------------------------------------------------------------------------

# ------------------------------------------------------------------------------

# Inputs

opt_Sh_S_list = []
opt_X_LEMAC_list = []
var_list = range(242, 454, 4)
varlist = []
for N_pas_below in var_list:
    print(N_pas_below)
    y_MAC_h, y_MAC_v, MAC_h, MAC_v, taper, c_root, X_LE_root, b_v, QC_sweep, mass_fractions, lg_weight, cg_locations, x_lemac, MAC, fuselage_design, W_TO, V_h_norm, HC_sweep, HC_sweep_h, b, S, l_h, x_le_h, sweep_LE_v, sweep_LE_h, x_le_v, b_v = main_iterator(CL_cruise_input, CD_cruise_input, W_e_frac_input, fuel_fractions_input, mass_fractions_input, N_pas_below)
    
    N_cargo = 2
    cargo_fwdfrac = 0.6  # Fraction of the amount the front compartment holds, if N_cargo = 1 then the value is 1
    
    Safety_margin = 0.10
    xcg_fuel      = 0.55 # %MAC
    l_nosecone = np.sum(fuselage_design[6]) / 2.
    xcg_nlg = l_nosecone  # assumption, use landing gear function !!!
    xcg_mlg = 0.75 * MAC   # assumption, use landing gear function !!!
    W_mlg   = 0.8*lg_weight
    W_nlg   = lg_weight - W_mlg
    
    xcg_seats = cg_seats(fuselage_design[16], fuselage_design[11], fuselage_design[12],l_nosecone)
    
    W_window, W_aisle, W_middle = W_seats(fuselage_design[8], fuselage_design[11], fuselage_design[12], fuselage_design[15])
    #print(W_window, W_aisle, W_middle)
    
    W_wing = mass_fractions[0]*W_TO
    W_ht = W_vt = mass_fractions[1]*W_TO
    W_fuse = mass_fractions[2]*W_TO
    W_nac = mass_fractions[3]*W_TO
    W_prop = mass_fractions[4]*W_TO
    W_fix = mass_fractions[5]*W_TO
    W_payload = mass_fractions[7]*W_TO
    W_fuel = mass_fractions[8]*W_TO
    
    
    W_personkg = W_person*lbs_to_kg*g_0
    
    xcg_fusegroup, xcg_winggroup, xcg_OEW_MAC, min_cg, max_cg, X_LEMAC_range, W_cargo, min_cg_range = potato(x_le_h, sweep_LE_h, y_MAC_h, MAC_h, x_le_v, sweep_LE_v, y_MAC_v, MAC_v, Safety_margin, xcg_fuel, W_fuel, W_fuse, W_nlg, W_vt, W_ht, W_fix, W_wing, W_nac, W_prop, W_mlg,
               cg_locations, xcg_nlg, xcg_mlg, x_lemac, W_payload, N_cargo, fuselage_design[7], MAC,
               cargo_fwdfrac, xcg_seats, W_window, W_aisle, W_middle, N_pas, W_personkg)
    
    M_cruise = M_cruise_input
    #------------------------------------------------------------------------------
    """STABILITY AND CONTROL"""
    
    M_h_cruise = M_cruise*V_h_norm
    #qcsweep    =  0.5515 #quarter chord sweep [rad]
    #hcsweepw    =  0.4869 #half chord sweep [rad]
    #A_w          =  8.67 #aspect ratio
    #S          =  427.8 #Wing surface
    #c          =  8.75 #mean aerodynamic chord
    #b          =  61.9 #wing span
    #c_r        =   12. #root chord
    c_f        =   0.25*MAC #flap chord
    #taper      =   0.149 #taper ratio
    b_f0       =   0.4*b #end of flap span
    b_fi       =   0.1*b #begin of flap span
    #l_f        =   56.4 #fuselage length
    #b_f        =   6.2 #fuselage diameter
    #h_f        =   6.2 #fuselage height
    #l_h        =   24.71 #tailarm
    #A_h        =   4.5 #aspect ratio horizontal tail
    #S_net      =   75. #area in fuselage
    
    """Coefficients"""
    Cm_0       =   -0.1 #airfoil dependent moment coefficient 
    CL_0       =   0.5 #CL of flapped wing at angle of attack of 0 deg 
    #CL_land    =   2.45 #CL during landing full flap configuration
    
    """Others"""
    #T0         =   293. #Temperature during landing  [K]
    #V_app      =   75.   #approach speed [m/s]
    delta_flap =  0.5235 #flap deflection [rad]
    alpha_land =   0.1745 #angle of attach during approach [rad]
    #M_app      =   0.22 #approach speed in Mach
    eta        =  0.95# airfoil efficiency = 0.95
    
    """Retrieve from stability cruve program"""
    x_ac       =  0.3   # aerodynamic centre location over the mac [n/m]
    #V_h_V      =  0.8    #(V_h/V) ratio  [m/s / m/s]
    #CL_alpha_AH =   4.788 # lif#t gradient of tailless aircraft [1/rad]
    #CL_AH       =  2. #lift coefficient A-H
    CL_H        = -0.8 #liftcoefficient adjustable tail
    
    
    
    #------------------------------VALUES GRAPHS--------------------------
    mu_2       = 1. #2D to 3D correction factor from graph  
    mu_3       = 0.025#  correction factor for sweep from graph
    dc_c_f     =  0.5 # flap geometry ratio, see torenbeek book
    
    #M_h_cruise = 0.7395
    #hcsweeph= 0.525
    #A_h=4.5
    #M_w_cruise= 0.87
    #phi=0.3
    #h_wh = 3
    #s_wTEh= 24.71
    b_n_1= 3.2
    b_n_2= 3.2
    b_n_3= 0
    b_n_4= 0
    l_n_1= -2
    l_n_2= -2
    l_n_3=0
    l_n_4=0
    x_ac_wing = 0.3
    #l_fn= 20
    c_g= MAC
    SM=0.10
    
    #----------------------------------------------------------------------
    
    x_cg = np.linspace(min(min_cg), max(max_cg))
    
    if wing_option == 0: 
        if tail_type == 0:
            m_tv = fuselage_design[1]/2.
        else:
            m_tv = fuselage_design[1]/2. + b_v
            l_h = l_h + (x_le_h-x_le_v) + np.tan(sweep_LE_v)*b_v
            
        
    else:
        if tail_type == 0:
            m_tv = -fuselage_design[1]/2. 
        else:
            m_tv = b_v - fuselage_design[1]/2.
            l_h = l_h + (x_le_h-x_le_v) + np.tan(sweep_LE_v)*b_v
            
           
    
    ## CHECK DOWNWASH AND CL ALPHA'S!!!!
    stability_lessmargin_list, stability_list = Sh_S_stability(x_cg, M_h_cruise, eta, HC_sweep_h, HC_sweep, A, A_h, M_cruise, fuselage_design[1], b, S,
                    l_h, QC_sweep, m_tv, V_h_norm, b_n_1, b_n_2, b_n_3, b_n_4, l_n_1, l_n_2, l_n_3, l_n_4, x_ac_wing, 
                    fuselage_design[1], X_LE_root, MAC, c_root, c_g, taper, SM)
    
      
    M_w_landing = V_stall_Landing/np.sqrt(gamma * R_gas * Temp_cruise)
    
    
    C_L_alpha_Ah_landing = C_L_alpha_Ah(M_w_landing, eta, HC_sweep, A, fuselage_design[1], b, S, c_root)  
    
    C_L_Ah_landing = C_L_alpha_Ah_landing*alpha_land
    S_netfus = S - (fuselage_design[1]*c_root)
    control_list = Sh_S_control(CL_H,C_L_Ah_landing,l_h,V_h_norm,x_cg,Cm_0,A,QC_sweep,
                 delta_flap,b,b_f0,b_fi,taper, MAC, c_f,dc_c_f,mu_2,mu_3,x_ac,CL_Landing_max,
                 fuselage_design[1], fuselage_design[1], fuselage_design[7],S,S_netfus,HC_sweep,M_w_landing,eta,CL_0)
    
    
    opt_X_LEMAC, opt_Sh_S = control_stability_plot(x_cg, min_cg, max_cg, X_LEMAC_range, control_list, stability_list)
    
    opt_X_LEMAC_list.append(opt_X_LEMAC)
    opt_Sh_S_list.append(opt_Sh_S)
    varlist.append(fuselage_design[7])

print(opt_X_LEMAC_list)
print(opt_Sh_S_list)


plt.subplot(2,1,1)
plt.plot(var_list, opt_Sh_S_list)
plt.subplot(2,1,2)
plt.plot(varlist, opt_Sh_S_list, 'bx')
plt.show()



W_personkg = W_person * lbs_to_kg * g_0








# ------------------------------VALUES GRAPHS--------------------------
mu_2 = 1.  # 2D to 3D correction factor from graph
mu_3 = 0.025  # correction factor for sweep from graph
dc_c_f = 0.5  # flap geometry ratio, see torenbeek book


