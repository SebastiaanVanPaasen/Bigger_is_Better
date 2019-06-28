# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 10:33:11 2019

@author: nikki
"""
import numpy as np
import matplotlib.pyplot as plt

#----------------------------------INPUTS--------------------------------------
#m_dot = mass flow rate through fan [kg/s]
#dT = temp. rise across the fan [K]
#N = #number of fan blades
#d_fan = [m] diameter of the fan 

#m_dot = 510.
#dT = 80.926
#N = 30.
#d_fan = 2.84124
#RSS = 100.
##distances
#r = 120.
#robs = 1.76   #height of observer 

m_dot = 535*0.85
dT = 100
N = 32
d_fan = 2.87
RSS = 100.
r = 300.
robs = 1.75
theta = 90.*np.pi/180.
#---------------------------------CONSTANT INPUTS------------------------------
#RSS = roto stator spacing/ fan blade chord, assume > 100, as there are not a lot of blades and a large fan diameter
m0_dot = 0.453 #[kg/s] reference mass flow rate
dT0 = 0.555 #[K] reference fan total temperature
RPM_fan = 2000.  
c = 343.2           #speed of sound [m/s]
M = 0.3

pe0 = 2e-05
W_ref = 10**(-12)  #reference acoustic power

#Octave band frequencies
freq = np.arange(50,10500,50)


#-------------------------------DEFINITIONS-----------------------------------
#F1 = function of the fanblade tip's Mach numbers
#F2 = function of the rotor spacing 
#F3 = function of directivity
#F4 = function of blade passing frequency
#F5 = constant associated with broadband noise
def SPL_fan(m_dot,m0_dot,dT,dT0,F1,F2,F3,F4):
    SPL = 10*np.log10(m_dot/m0_dot) + 20*np.log10(dT/dT0) + F1 + F2 + F3 + F4
    return SPL
    
#--------------------------------NOISE------------------------------------
#To compute the Mach tip numbers
s_tip = np.pi*d_fan
V_tip = RPM_fan*s_tip/60.
M_tip = V_tip/c

#To compute fan blade passing frequency BPF
fb = N*RPM_fan*(1./60.)  #in Hertz

SPL_fan_list = []
SPL_fan_effects = []


for f in freq:
    """Inlet broadband noise"""
    if M_tip < 0.72:
        F1 = 34.
    elif M_tip >= 0.72:
        F1 = 34. - 43.*(M_tip - 0.72)
    
    F2 = -5*np.log10(RSS/300.)
    F3 = -18.5 
    F5 = -0.5* (np.log(f/(4.*fb))/np.log(2.2))**2
    F4 = 10*np.log10(np.exp(-F5))

    SPL_inlet_broad = SPL_fan(m_dot,m0_dot,dT,dT0,F1,F2,F3,F4)

    """Inlet tonal noise"""
    F1 = 42. - 20.*M_tip
    F2 = -10*np.log10(RSS/300.)
    F3 = -15.5
    F4 = 0.
    
    SPL_inlet_tonal = SPL_fan(m_dot,m0_dot,dT,dT0,F1,F2,F3,F4)

    """Combination tones 1/2 BPF"""
    if M_tip < 1.146:    
        F1 = -18. + 46.5*(M_tip - 1.)/0.146
    elif M_tip >= 1.146:
        F1 = 28.5 - 12.*(M_tip - 1.146)/0.854
    
    F2 = 0.
    F3 = -3.9
    
    if f < 0.5*fb:
        F4 = 20.*np.log10(f/(0.5*fb))
    elif f >= 0.5*fb:
        F4 = -20.*np.log10(f/(0.5*fb))
        
    SPL_comb_12 = SPL_fan(m_dot,m0_dot,dT,dT0,F1,F2,F3,F4)
    
    """Combination tones 1/4 BPF"""
    if M_tip < 1.322:
        F1 = -15. + 47.5*(M_tip - 1.)/0.322
    elif M_tip >= 1.322:
        F1 = 32.5 - 9.*(M_tip - 1.322)/0.678
    
    F2 = 0.
    F3 = -3.9
    
    if f < 0.25*fb:
        F4 = 30.*np.log10(f/(0.25*fb))
    elif f >= 0.25*fb:
        F4 = -30.*np.log10(f/(0.25*fb))
        
    SPL_comb_14 = SPL_fan(m_dot,m0_dot,dT,dT0,F1,F2,F3,F4)
    
    """Combination tones 1/8 BPF"""
    if M_tip < 1.61:
        F1 = -12. + 4.7*(M_tip - 1.)/0.61
    elif M_tip >= 1.61:
        F1 = 29.2 - 4.7*(M_tip - 1.61)/0.39
        
    F2 = 0.
    F3 = -3.9
    
    if f < 0.225*fb:
        F4 = 30.*np.log10(f/(0.125*fb))
    elif f >= 0.125*fb:
        F4 = -20.*np.log10(f/(0.125*fb))
        
    SPL_comb_18 = SPL_fan(m_dot,m0_dot,dT,dT0,F1,F2,F3,F4)
    
    """Aft broadband noise"""
    F1 = 34. - 17.*(M_tip - 0.65)
    F2 = -5.*np.log10(RSS/300)
    F3 = -6.9
    F5 = -0.5*(np.log(f/(2.5*fb))/np.log(2.2))**2
    F4 = 10.*np.log10(np.exp(-F5))
    
    SPL_aft_broad = SPL_fan(m_dot,m0_dot,dT,dT0,F1,F2,F3,F4)    
    
    """Aft tonal noise"""
    F1 = 46. - 20.*M_tip
    F2 = -10.*np.log10(100./300.)
    F3 = -7.1
    F4 = 0.
    
    SPL_aft_tonal = SPL_fan(m_dot,m0_dot,dT,dT0,F1,F2,F3,F4)
    
    """Superposition to sum fan contributions"""
    SPL_fantot = 10.*np.log10(10.**(SPL_inlet_broad/10.)+10.**(SPL_inlet_tonal/10.)+10.**(SPL_comb_12/10.)+10.**(SPL_comb_14/10.)+10.**(SPL_comb_18/10.) + 10.**(SPL_aft_broad/10.) + 10.**(SPL_aft_tonal/10.))
    SPL_corrected = SPL_fantot + 20.*np.log10(1./r)  #correct for observer distance 

    """Account for ground and atmospheric effects"""
    #Atmospheric absorption from reader of master course Table 3, part II, chapter 3
    alpha = 1e-09*f**2 + 1e-06*f + 0.0008
    SPL_absorp = alpha*r
    
    #Ground effect noise
    dSPL = 10.*np.log10(abs(2 + 2*np.cos((4*np.pi*f)/c)*robs*np.sin(theta)))
   
    """Final values"""
    SPL_fan_list.append(SPL_corrected - 15.)
    SPL_fan_effects.append(SPL_corrected - SPL_absorp + dSPL - 15.)     #reduction of 15 dB due to geared fan and acoustic lining in fan inlet!!
    
    
plt.figure(7) 
plt.plot(freq,SPL_fan_list)
#plt.title("Fan noise") 
plt.grid(True)
plt.xlabel("1/3 Octave Band Centre Frequency [Hz]",fontsize ='x-large')
plt.ylabel("1/3 Octave Band SPL [dB]",fontsize ='x-large')
plt.xlim(50,10000)
plt.ylim(20,100)
plt.xscale('log')


plt.figure(8)
plt.plot(freq,SPL_fan_effects)
#plt.title("Fan noise including ground effect and atmos. absorption") 
plt.grid(True)
plt.xlabel("1/3 Octave Band Centre Frequency [Hz]",fontsize ='x-large')
plt.ylabel("1/3 Octave Band SPL [dB]",fontsize ='x-large')
plt.xlim(50,10000)
plt.xscale('log')
plt.ylim(0,100)


plt.show()

    

























