# -*- coding: utf-8 -*-
"""
Created on Fri May 10 09:20:29 2019

@author: mathi
"""
import numpy as np
import math
import scipy as sp
import matplotlib as plt

from control_curve import Sh_S_control
from stability import Sh_S_stability

def control_stability_plot(CL_h,CL_ah,l_h,Cm_ac,x_ac,x_cg, M_h_cruise, eta, hcsweeph, hcsweepw, A_w, A_h, M_w_cruise, b_f, b, S, l_h, qcsweep, phi, h_wh, s_wTEh, V_h_V,b_n_1, b_n_2, b_n_3, b_n_4, l_n_1, l_n_2, l_n_3, l_n_4, x_ac_wing, h_f, l_fn, c, c_r, c_g, taper, SM):
    
    S_control = Sh_S_control(CL_h,CL_ah,l_h,V_h_V,c,Cm_ac,x_ac,x_cg)
    
    S_stability = Sh_S_stability(x_cg, M_h_cruise, eta, hcsweeph, hcsweepw, A_w, A_h, M_w_cruise, b_f, b, S, l_h, qcsweep, phi, h_wh, s_wTEh, V_h_V,b_n_1, b_n_2, b_n_3, b_n_4, l_n_1, l_n_2, l_n_3, l_n_4, x_ac_wing, h_f, l_fn, c, c_r, c_g, taper, SM)
    
    # Calculating the trendline
    contr_trend = np.polyfit(x_cg, S_control, 1)
    stability_trend = np.polyfit(x_cg, S_stability, 1)
    a_control = contr_trend[0]
    b_control = contr_trend[1]
    a_stability = stability_trend[0]
    b_stability = stability_trend[1]
    
    margin = []
    Sh_S = np.linspace(0,1,100)
    for i in range(len(Sh_S)):
        x1 = (Sh_S[i] - b_control)/a_control
        x2 = (Sh_S[i] - b_stability)/a_stability
        margin_sc = x1-x2
        margin.append(margin)
    
    plt.plot(x_cg, S_control, 'r')
    plt.plot(x_cg, S_stability)
    plt.show()
            

        
    