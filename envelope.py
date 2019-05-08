# -*- coding: utf-8 -*-
"""
Created on Mon May  6 13:32:23 2019

@author: Mathilde

"""

import numpy as np
import matplotlib.pyplot as plt

def envelope(W,rho,C_L_max_pos, C_L_max_neg,S,V_C): 

#construct manouevre plot
    lb_to_kg = 0.454
    n_top = [0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,2.5,0]
    n_bottom = [0,-0.2,-0.4,-0.6,-0.8,-1,-1,0]
    V_top = []
    V_bottom = []
    
    V_C = 250
    n_max = 2.5
    n_min = -1
    V_S1 = np.sqrt((2*W)/(rho*C_L_max_pos*S))
    V_H = np.sqrt((2*W)/(rho*C_L_max_neg*S))
    for i in range(len(n_top)-2):
        
        V_Sn = V_S1*np.sqrt(n_top[i])
        V_top.append(V_Sn)
        
        
    for i in range(len(n_bottom)-2):  
        V_Hn = V_H*np.sqrt(abs(n_bottom[i]))
        V_bottom.append(V_Hn)
    
    V_bottom.append(V_C)
    V_D = V_C/0.8
    V_bottom.append(V_D)
    V_top.append(V_D)
    V_top.append(V_D)
    
    plt.plot(V_top,n_top)
    plt.plot(V_bottom, n_bottom)
    plt.show()
    
#construct gust loading plot
    #ρ0= air density at sea level (kg/m3)
    #a=CLα lift slope coefficient (1/rad)
    #UEAS= equivalent gust speed (m/s)
    #VEAS= equivalent aircraft speed (m/s)
    #Kg= gust alleviation coefficient (as function of GH) with 
    #w= aircraft wing loading [W/S] (N/m2)
    #c =  mean geometric chord (m) 
    #ρ = density at altitude (kg/m3)
    #g =acceleration due to gravity (m/s2) 
    
    
    V_gust= [V_S,V_B,V_C,V_D]
    h = 10000
    #interpolation
    alt = 15250-6100
    slope_B = (11.5-20.12)/alt
    slope_C = (7.62-15.25)/alt
    slope_D = (3.18-7.62)/alt
    
    gust_B =20.12- slope_B*(h-6100) 
    gust_C =15.25- slope_C*(h-6100) 
    gust_D =7.62- slope_D*(h-6100)
    
    mu_g = 2*w/(rho*c*a*g)
    K_g = 0.88*mu_g/(5.3+mu_g)
    delta_n_g = 0.5*(rho_0*a)/(m*g/S)*U_EAS*V*K_g*mu_g
    


    
    


