# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 15:50:09 2019

@author: nikki
"""
import numpy as np
import matplotlib.pyplot as plt
#-----------------------------------INPUTS NEEDED-----------------------------
#A_w = wing area [m^2]
#b_w = wingspan [m]
#A_h = horizontal tail area [m^2]
#b_h = horizontal tail span [m]
#A_f = flap area [m^2]
#b_f = flap span [m]
#d = wheel diameter [m]
#n_gear = number of wheels in gear
#f = frequency [Hz], consider range 50-10000 Hz
#theta_flap = 30*np.pi/180 #flap deflection in rad
#r = distance to oberserver as prescribed by the noise measurement points

robs = 1.8 #height of the observer

"""Input values of the B737-800"""
A_w = 124.6
b_w = 34.3
A_h = 32.40
b_h = 13.4
A_f = 33.87
b_f = 0.599*b_w
d = 1.016
n_gear = 4
theta_flap = 30*np.pi/180.
r = 120. 


#-----------------------------------DEFINITIONS------------------------------
"""Effective noise power """
#c = speed of sound
#D = directivity factor
#F = emperical spectral function
#P = power function [W]
def pe_2(rho,M,c,r,P,D,F,theta):
    theta = theta*(np.pi/180.)
    pe_2 = (rho*c*P*D*F)/(4*np.pi*r**2*(1. - M*np.cos(theta))**4)
    return pe_2
    
"""Power function"""
#K&a = dimensionless constants based on emperical data
#G = geometry function
#b_w = wing span [m]
def Power(k,a,M,G):
    P = k*(M**a)*G  
    return P
    
"""Strouhal number"""
#L = length scale characteristics of component considered
#f = frequency considered
def Strouhal(f,L,M,theta,c):
    theta = theta*(np.pi/180.)
    S = (f*L*(1-M*np.cos(theta)))/(M*c)
    return S
    
    
#-------------------------GENERAL CONSTANT INPUTS------------------------------
mu = 1.84e-05       #dynamic viscosity kg/ms 
rho = 1.225
c = 343.2           #speed of sound [m/s]
M = 0.3             #Mach number during approach   
theta = 90.         #emission angle (overhead so 90 degrees)
phi = 0.             #azimuth angle (overhead so 0)
pe0 = 2e-05         #Reference effective pressure [N/m^2] [Pa]


#-------------------------PER COMPONENT NOISE----------------------------------
freq = np.arange(50,10500,500)
SPL = []


for f in freq:
    """TE clean wing"""
    k_w = 4.464e-05
    a_w = 5.
    
    G_wing = 0.37*(A_w/b_w)*((rho*M*c*A_w)/(mu*b_w))**-0.2
    L_wing = G_wing*b_w
    
    #Strouhal number
    S_wing = Strouhal(f,L_wing,M,theta,c)
    #Power function
    P_wing = Power(k_w,a_w,M,G_wing)
    #Spectral function
    F_wing = 0.613*(10*S_wing)**4 * ((10*S_wing)**1.5 + 0.5)**-4
    #Directivity function
    D_wing = 4*np.cos(phi)**2*np.cos(theta/2.)**2
    
    #Effective power
    pe_wing = pe_2(rho,M,c,r,P_wing,D_wing,F_wing,theta)
    
    
    """Horizontal tailplane"""
    k_h = 4.464e-05
    a_h = 5.
    
    G_tail = 0.37*(A_h/b_h)*((rho*M*c*A_h)/(mu*b_h))**-0.2
    L_tail = G_tail*b_h
    
    #Strouhal number
    S_tail = Strouhal(f,L_tail,M,theta,c)
    #Power function
    P_tail = Power(k_h,a_h,M,G_tail)
    #Spectral function
    F_tail = 0.613*(10*S_tail)**4 * ((10*S_tail)**1.5 + 0.5)**-4
    #Directivity function
    D_tail = 4*np.cos(phi)**2*np.cos(theta/2.)**2
    
    #Effective power
    pe_tail = pe_2(rho,M,c,r,P_tail,D_tail,F_tail,theta)
    
    
    """LE slats"""
    k_slats = 4.464e-05
    a_slats = 5.
    
    G_slats = 0.37*(A_w/b_w)*((rho*M*c*A_w)/(mu*b_w))**-0.2
    L_slats = G_slats*b_w
    
    #Strouhal number
    S_slats = Strouhal(f,L_slats,M,theta,c)
    #Power function
    P_slats = Power(k_slats,a_slats,M,G_slats)
    #Spectral function
    F_slats = 0.613*(10*S_slats)**4 * ((10*S_slats)**1.5 + 0.5)**(-4) + 0.613*(2.19*S_slats)**4 * ((2.19*S_slats)**1.5+ 0.5)**(-4)
    #Directivity function
    D_slats = 4*np.cos(phi)**2*np.cos(theta/2.)**2
    
    #Effective power
    pe_slats = pe_2(rho,M,c,r,P_slats,D_slats,F_slats,theta)
    
    
    
    """TE flaps"""
    k_f = 2.787e-04
    a_f = 6.
    
    G_flaps = (A_f/b_w**2)*np.sin(theta_flap)**2
    L_flaps = A_f/b_f
    
    #Strouhal number
    S_flaps = Strouhal(f,L_flaps,M,theta,c)
    #Power function
    P_flaps = Power(k_f,a_f,M,G_flaps)
    #Spectral function
    if S_flaps < 2.:
        F_flaps = 0.0480*S_flaps
    elif 2. <= S_flaps <= 20.:
        F_flaps = 0.1406*(S_flaps)**(-0.55)
    elif S_flaps > 20:
        F_flaps = 216.49*S_flaps**(-3)
    #Directivity function
    D_flaps = 3.*(np.sin(theta_flap)*np.cos(theta) + np.cos(theta_flap)*np.sin(theta)*np.cos(phi))**2
    
    #Effective power
    pe_flaps = pe_2(rho,M,c,r,P_flaps,D_flaps,F_flaps,theta)
    
    """Landing gear"""
    if n_gear <= 2:
        k_gear = 4.349e-04
    else:
        k_gear = 3.414e-04
    
    a_gear = 6
        
    G_gear = n_gear*(d/b_w)**2    
    L_gear = d
    
    #Strouhal number
    S_gear = Strouhal(f,L_gear,M,theta,c)
    #Power function
    P_gear = Power(k_gear,a_gear,M,G_gear)
    #Spectral function
    if n_gear <= 2.:
        F_gear = 13.595*S_gear**2 * (S_gear**2 + 12.5)**(-2.25)
    else:
        F_gear = 0.0577*S_gear**2 * (0.25*S_gear**2 + 1)**(-1.5)
    #Directivity function
    D_gear = (3./2.)*np.sin(theta)**2
    
    #Effective power
    pe_gear = pe_2(rho,M,c,r,P_gear,D_gear,F_gear,theta)
    
    
    """Total sound pressure level SPL in dB"""
    #Use superposition to sum contributions
    SPL_airframe = 10.*np.log10((pe_wing**2 + pe_tail**2 + pe_slats**2 + pe_flaps**2 + pe_gear**2)/pe0**2)

    #Atmospheric absorption from reader of master course Table 3, part II, chapter 3
    alpha = 1e-09*f**2 + 1e-06*f + 0.0008
    SPL_absorp = alpha*r
    
    #Ground effect noise
    dSPL = 10.*np.log10(abs(2 + 2*np.cos((4*np.pi*f)/c)*robs*np.sin(theta)))
       
    SPL.append(SPL_airframe + dSPL - SPL_absorp)
    
     
plt.plot(freq,SPL)
plt.title("Airframe noise including ground effect, doppler, directivity and absorption")
plt.grid(True)
plt.xlabel("Frequency [Hz]")
plt.ylabel("SPL [dB]")
plt.show()


#TRANSFORM INTO PNL AND TEST IT USING THE FAA NOISE DATA SHEETS WITH THE AIRCRAFT IN THERE AND COMPARE VALUES








    
    
    
    