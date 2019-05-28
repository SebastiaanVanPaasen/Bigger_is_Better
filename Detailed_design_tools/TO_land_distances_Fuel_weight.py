# -*- coding: utf-8 -*-
"""
Created on Mon May 27 17:24:10 2019

@author: nikki
"""

#----------------------------IMPORT MODULES-----------------------------------
import numpy as np


#---------------------------DEFINITIONS----------------------------------------



"""Fuel weight: cruise and reserve fuel"""
#Inputs 
#R      =       Design range in [m]
#Ct     =       Thrust specific fuel consumption [kg/N/s]
#g      =       9.81 [m/s^2]
#T      =       Static air temperature at atitude [K]
#T0     =       Static air temperature at sea level [K]
#a0     =       Speed of sound at sea level
#M      =       Mach number 
#CD     =       Drag coefficient
#CL     =       Lift coefficient 
#A      =       Aspect ratio


def Wf_Wto():       #To compute the fuel weight for the trip
    
    #Assumptions:
    #Cruise fuel is the dominant fuel weight
    #Wetted area of the fuselage is the main contributer for fuselage drag   
    
    theta = T/T0
    Wfcr_Wto = 1. - np.exp**(((-R*Ct*g/np.sqrt(theta))/(a0*M))*(CD/CL))
    Wfres_Wto = 0.18 * ((Ct*g/np.sqrt(theta))/np.sqrt(A))
    
    Wf_Wto = Wfcr_Wto + Wfres_Wto
    return Wf_Wto
    
    
"""Take-off distance: all engines functioning"""
#Inputs
#W_to       =       Take-off weight (MTOW)
#S          =       Wing surface area [m^2]
#rho0       =       Sea level density [kg/m^3]
#CL_maxto   =       CL max during take-off
#bypass     =       Bypass ratio of engine
#T_to       =       Take-off thrust [N]
#g          =       9.81 [m/s^2]
#A          =       Aspect ratio

def TO_distance():       #Required distance to pass screen height (30 ft) at a speed of 1.3*VV_stall
    Vs = np.sqrt((W_to*2.)/(S*rho0*CL_maxto))    
    V_LOF = 1.2*Vs
    T_mean = 0.75*((5. + bypass)/(4.+ bypass))*T_to    
    mu_dash = 0.02 + 0.01*CL_maxto
    
    S_run = (V_LOF**2/(2.*g))/((T_mean/W_to)-mu_dash)
    
    gamma_LOF = 0.9*(T_mean/W_to) - (0.3/np.sqrt(A))
    h_to = 10.7 #[m] screen height of 35 ft
    S_air = (V_LOF**2)/(g*np.sqrt(2)) + (h_to/gamma_LOF)
    
    S_to = S_run + S_air
    
    return S_to
    
"""Take-off W/S limit: for given take-off distance"""
#Inputs
#W_to       =       Take-off weight (MTOW)
#S          =       Wing surface area [m^2]
#rho0       =       Sea level density [kg/m^3]
#CL_maxto   =       CL max during take-off
#g          =       9.81 [m/s^2]
#A          =       Aspect ratio


def WS_TO():
    Vs = np.sqrt((W_to*2.)/(S*rho0*CL_maxto))    
    V_LOF = 1.2*Vs
    gamma_LOF = 0.9*(T_mean/W_to) - (0.3/np.sqrt(A))
    V3 = V_LOF * np.sqrt(1. + gamma_LOF*np.sqrt(2.))
    mu_dash = 0.02 + 0.01*CL_maxto
    
    h_to = 10.7 #[m] screen height of 35 ft
    f_to = 1. #or 1.15
    
    a = (S_to/f_to) - (h_to/gamma_LOF)
    b = rho0*g*CL_maxto * (1. + gamma_LOF*np.sqrt(2.))
    c = (V3/Vs)**2 * ((T_mean/W_to-mu_dash)**(-1.) + np.sqrt(2.))

    WS = a * (b/c)
    return WS
 
    
"""TO with engine failure"""
#Inputs
#W_to       =       Take-off weight (MTOW)
#T_to       =       Take-off thrust [N]
#S          =       Wing surface area [m^2]
#rho0       =       Sea level density [kg/m^3]
#CL_maxto   =       CL max during take-off
#g          =       9.81 [m/s^2]
#A          =       Aspect ratio
#CD_to      =       Average CD during take-off
#bypass     =       Bypass ratio of engine
#a_stop     =       Average decelearation usually 0.37g
#a_mean     =       Mean acceleration assume 0.25g??
#h_to       =       10.7 [m]

#gamma2_min =       Correction factor
#           =       0.024; 0.027; 0.030
#           =       2; 3; 4 number of engines

#dt         =       4.5 [s] reaction time
#mu         =       Ground friction coefficient
#           =       0.02 for concrete

def TO_eng_fail():
    Vs = np.sqrt((W_to*2.)/(S*rho0*CL_maxto))    
    V2 = 1.2*Vs
    CL_to = mu*np.pi*A*e
    
    dgamma2 = ((T_to/W_to)-(CL_to/CD_to)**(-1.)) - gamma2_min
    gamma_mean = 0.06 + dgamma2
    
    Vx = V2*( ((1. + 2.*g*h_to/V2**2)/(1. + gamma_mean / (a_mean/g)))**0.5 - ((gamma_mean*g*(dt - 1.))/V2) )
    
    S01 = Vx**2 / (2.*a_mean)                               #Distance covered before engine failure at Vx
    S12 = (1./gamma_mean)*(((V2**2-Vx**2)/(2.*g)) + h_to)   #Distance from engine failure up to save screen height at V2
    Sstop = (Vx**2/(2*a_stop)) + Vx*dt                      #If TO aborted (stop distance needed)
    
    S_continue = S01 + S12
    S_abord = S01 + Sstop
    
    return S_continue, S_abord

def BFL():      #Balanced field length: is when taking off is better than stopping in case of engine failure
    CL_to = mu*np.pi*A*e    
    dgamma2 = ((T_to/W_to)-(CL_to/CD_to)**(-1.)) - gamma2_min    
    CL2 = 0.694*CL_maxto
    T_mean = 0.75*((5. + bypass)/(4.+ bypass))*T_to    
    mu_dash = 0.02 + 0.01*CL_maxto
    
    a = 0.863/(1. + 2.3*dgamma2)    
    b = ((W_to/S)/(rho0*g*CL2)) + h_to    
    c = (1./(T_mean/W_to - mu_dash)) + 2.7
    
    BFL = a*b*c
    return BFL
    
    
"""Landing distance"""
#Inputs
#gamma_mean =       During landing is approx. 0.1
#h_land     =       15.3 [m] screen height for landing
#dn         =       0.1 (change in load factor during landing)
#a_stop/g   =       0.35-0.45 without thrust reversers
#           =       0.40 - 0.50 thrustreversers and spoilers
#           =       0.50 - 0.60 + nose wheel braking

def S_land():
    Va = 1.3*Vs
    Vtd = np.sqrt(Va**2*(1. - (gamma_mean**2/dn)))    

    S_air = (1./gamma_mean)*(((Va**2 - Vtd**2)/(2.*g)) + h_land)
    S_run = Vtd**2 / (2.*a_stop)
    
    SL = S_air + S_run
    return SL
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




