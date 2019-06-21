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

robs = 1.75 #height of the observer

A_w = 35.9
b_w = 35.9


A_h = 32.78
b_h = 14.35

b_f = 0.599*b_w
A_f = 0.3*A_w

d_main = 44.5*2.54/100.
d_nose = 27*2.54/100.

n_gear = 4
n_nose = 2

theta_flap = 5.*np.pi/180.
V =74*1.5

r =300.

theta = 90.*np.pi/180.         #emission angle (overhead so 90 degrees)
phi = 0. #90.*np.pi/180. 

#"""Input values of the fokker 70"""
#A_w = 81.
#b_w = 28.
#
#A_h = 21.72
#b_h = 10.04
#
#A_f = 6.5
#b_f = 16.
#d = 1.1
#
#n_gear = 4.
#n_nose = 2.
#
#theta_flap = 30*np.pi/180.
#V = 68.
#r = 63.5 


#-------------------------GENERAL CONSTANT INPUTS------------------------------
mu = 1.84e-05       #dynamic viscosity kg/ms 
rho = 1.225
c = 343.2           #speed of sound [m/s]
M = V/c             #Mach number during approach   
           #azimuth angle (overhead so 0)
pe0 = 2e-05         #Reference effective pressure [N/m^2] [Pa]



#-----------------------------------DEFINITIONS------------------------------
"""Effective noise power """
#c = speed of sound
#D = directivity factor
#F = emperical spectral function
#P = power function [W]
def pe_2(rho,M,c,r,P,D,F,theta):
    pe_2 = (rho*c*P*D*F)/(4*np.pi*(r**2)*(1. - M*np.cos(theta))**4)
    return pe_2
    
"""Power function"""
#K&a = dimensionless constants based on emperical data
#G = geometry function
#b_w = wing span [m]
def Power(k,a,M,G):
    P = k*(M**a)*G*(rho*(c**3)*b_w**2)  
    return P
    
"""Strouhal number"""
#L = length scale characteristics of component considered
#f = frequency considered
def Strouhal(f,L,M,theta,c):
    S = (f*L*(1. - M*np.cos(theta)))/(M*c)
    return S
    
#-------------------------PER COMPONENT NOISE----------------------------------
freq = np.arange(50,10500,50)

SPL = []
SPL_wing = []
SPL_tail = []
SPL_slats = []
SPL_flaps = []
SPL_gear = []
SPL_nose = []

SPL_cor = []
SPL_wing_cor = []
SPL_tail_cor = []
SPL_slats_cor = []
SPL_flaps_cor = []
SPL_gear_cor = []
SPL_nose_cor = []

F_list = []
S_list = []

for f in freq:
    """TE clean wing"""
    k_w = 4.464e-05
    a_w = 5.
    
    G_wing = 0.37*(A_w/b_w**2)*((rho*M*c*A_w)/(mu*b_w))**(-0.2)
    L_wing = G_wing*b_w
    
    #Strouhal number
    S_wing = Strouhal(f,L_wing,M,theta,c)
    #Power function
    P_wing = Power(k_w,a_w,M,G_wing)
    #Spectral function
    F_wing = 0.613*(10*S_wing)**4 * ((10*S_wing)**1.5 + 0.5)**(-4)
    #Directivity function
    D_wing = 4*np.cos(phi)**2*np.cos(theta/2.)**2
    
    #Effective power
    pe_wing = pe_2(rho,M,c,r,P_wing,D_wing,F_wing,theta)
    SPL_wingi = 10*np.log10(pe_wing/pe0**2)
    PSL_wing = 10*np.log10(P_wing/pe0**2)

    """Horizontal tailplane"""
    k_h = 4.464e-05
    a_h = 5.
    
    G_tail = 0.37*(A_h/b_h**2)*((rho*M*c*A_h)/(mu*b_h))**(-0.2)
    L_tail = G_tail*b_h
    
    #Strouhal number
    S_tail = Strouhal(f,L_tail,M,theta,c)
    #Power function
    P_tail = Power(k_h,a_h,M,G_tail)
    #Spectral function
    F_tail = 0.613*(10*S_tail)**4 * ((10*S_tail)**1.5 + 0.5)**(-4)
    #Directivity function
    D_tail = 4*np.cos(phi)**2*np.cos(theta/2.)**2
    
    #Effective power
    pe_tail = pe_2(rho,M,c,r,P_tail,D_tail,F_tail,theta)
    SPL_taili = 10*np.log10(pe_tail/pe0**2)
    PSL_tail = 10*np.log10(P_tail/pe0**2)
    
    """LE slats"""
    k_slats = 4.464e-05
    a_slats = 5.
    
    G_slats = 0.37*(A_w/b_w**2)*((rho*M*c*A_w)/(mu*b_w))**(-0.2)
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
    SPL_slatsi = 10*np.log10(pe_slats/pe0**2)
    PSL_slats = 10*np.log10(P_slats/pe0**2)
    
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
    SPL_flapsi = 10*np.log10(pe_flaps/pe0**2)
    PSL_flaps = 10*np.log10(P_flaps/pe0**2)
    
    
    """Main landing gear"""
    if n_gear == 2:
        k_gear = 4.349e-04
    else:
        k_gear = 3.414e-04
    
    a_gear = 6
        
    G_gear = n_gear*(d_main/b_w)**2    
    L_gear = d_main
    
    #Strouhal number
    S_gear = Strouhal(f,L_gear,M,theta,c)
    #Power function
    P_gear = Power(k_gear,a_gear,M,G_gear)
    #Spectral function
    if n_gear == 2:
        F_gear = 13.595*S_gear**2 * (S_gear**2 + 12.5)**(-2.25)
    else:
        F_gear = 0.0577*S_gear**2 * (0.25*S_gear**2 + 1)**(-1.5)
    #Directivity function
    D_gear = (3./2.)*np.sin(theta)**2
    
    #Effective power
    pe_gear = pe_2(rho,M,c,r,P_gear,D_gear,F_gear,theta)
    SPL_geari = 10*np.log10(pe_gear/pe0**2)
    PSL_gear = 10*np.log10(P_gear/pe0**2)
   
    #To verify   
    F_list.append(np.log10(F_gear))
    S_list.append(np.log10(S_gear))
    
    
    """Nose gear"""
    if n_nose == 2:
        k_nose = 4.349e-04
    else:
        k_nose = 3.414e-04
    
    a_nose = 6
        
    G_nose = n_nose*(d_nose/b_w)**2    
    L_nose = d_nose
    
    #Strouhal number
    S_nose = Strouhal(f,L_nose,M,theta,c)
    #Power function
    P_nose = Power(k_nose,a_nose,M,G_nose)
    #Spectral function
    if n_nose == 2:
        F_nose = 13.595*S_nose**2 * (S_nose**2 + 12.5)**(-2.25)
    else:
        F_nose = 0.0577*S_nose**2 * (0.25*S_nose**2 + 1)**(-1.5)
    #Directivity function
    D_nose = (3./2.)*np.sin(theta)**2
    
    #Effective power
    pe_nose = pe_2(rho,M,c,r,P_nose,D_nose,F_nose,theta)
    SPL_nosei = 10*np.log10(pe_nose/pe0**2)
    PSL_nose = 10*np.log10(P_nose/pe0**2)
   


    
    """Total sound pressure level SPL in dB"""  
    #Use superposition to sum contributions
    SPL_airframe = 10.*np.log10((pe_wing + pe_tail + pe_slats + pe_flaps + pe_gear)/pe0**2)

    #Atmospheric absorption from reader of master course Table 3, part II, chapter 3
    alpha = 1e-09*f**2 + 1e-06*f + 0.0008
    SPL_absorp = alpha*r
    
    #Ground effect noise
    dSPL = 10.*np.log10(abs(2 + 2*np.cos((4*np.pi*f)/c)*robs*np.sin(theta)))

    #Total and component SPL values    
    SPL.append(SPL_airframe ) 
    SPL_wing.append(SPL_wingi )
    SPL_tail.append(SPL_taili )
    SPL_slats.append(SPL_slatsi )
    SPL_flaps.append(SPL_flapsi)
    SPL_gear.append(SPL_geari )
    SPL_nose.append(SPL_nosei )

    #Total and component SPL values including correction values
    SPL_cor.append(SPL_airframe + dSPL - SPL_absorp)
    SPL_wing_cor.append(SPL_wingi + dSPL - SPL_absorp)
    SPL_tail_cor.append(SPL_taili + dSPL - SPL_absorp)
    SPL_slats_cor.append(SPL_slatsi + dSPL - SPL_absorp)
    SPL_flaps_cor.append(SPL_flapsi + dSPL - SPL_absorp)
    SPL_gear_cor.append(SPL_geari + dSPL - SPL_absorp)
    SPL_nose_cor.append(SPL_nosei + dSPL - SPL_absorp)
   
plt.figure(1)  
plt.plot(freq,SPL,"k",label = "Total")
plt.plot(freq,SPL_wing,"b",label= "Wing")
plt.plot(freq,SPL_tail,"g--",label = "Horizontal tail")
plt.plot(freq,SPL_slats,"r", label = "Slats")
plt.plot(freq,SPL_flaps,"c--", label = "Flaps")
plt.plot(freq, SPL_gear,"orange",label = "Main gear")
plt.plot(freq, SPL_nose,"m--",label = "Nose gear")
#plt.plot(freq,SPL_cor,"y", label = "Tot corr.")

#plt.title("Airframe noise") 
plt.grid(True)
plt.xlabel("1/3 Octave Band Centre Frequency [Hz]",fontsize ='x-large')
plt.ylabel("1/3 Octave Band SPL [dB]",fontsize ='x-large')
plt.xlim(50,10000)
plt.xscale('log')
plt.legend(fontsize ='x-large')

plt.figure(2)
plt.plot(freq,SPL_cor,"k",label = "Total")
plt.plot(freq,SPL_wing_cor,"b",label= "Wing")
plt.plot(freq,SPL_tail_cor, "g--",label = "Horizontal tail")
plt.plot(freq,SPL_slats_cor, "r", label = "Slats")
plt.plot(freq,SPL_flaps_cor, "c--", label = "Flaps")
plt.plot(freq, SPL_gear_cor,"orange",label = "Main gear")
plt.plot(freq, SPL_nose_cor,"m--",label = "Nose gear")

#plt.title("Airframe noise including ground effect and atmospheric absorption") 
plt.grid(True)
plt.xlabel("1/3 Octave Band Centre Frequency [Hz]",fontsize ='x-large')
plt.ylabel("1/3 Octave Band SPL [dB]",fontsize ='x-large')
plt.xlim(50,10000)
plt.xscale('log')
plt.legend(fontsize ='x-large')

plt.show()


#plt.plot(S_list,F_list)
##plt.xscale('log')
##plt.yscale('log')
#plt.show()







    
    
    
    