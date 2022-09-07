# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 15:15:36 2019

@author: nikki
"""

import numpy as np
import matplotlib.pyplot as plt

#------------------------------INPUTS--------------------------------------------
MTOW       =       242670*9.81
MFW        =       71170.*9.81
Wcr        =       MTOW-0.4*MFW 
T0         =       342.5*1000*2 #127.62*1000  #@ SEA LEVEL!!!
T_climb    =       342.5*1000*2

Vcr        =       256.7

A          =       9.44       
e          =       0.85 
S          =       127. 
CD0        =       0.020



"""NOTE THIS VALUE DEPENDS ON THE ENGINE TYPE"""
m          =       1.3 #factor for the variation of thrust with alitude

#----------------------------DEFINITIONS--------------------------------------      
""" ISA definitions""" 
def ISA_density(h):      # enter height in m
    M = 0.0289644       #kg/mol molar mass of Earth's air
    R = 8.3144590       #universal gas constant Nm/MolK
    g = 9.81
    if h < 11000:
        rho0 = 1.225   #kg/m^3
        T = 288.15     #K
        h0 = 0.         #m
        a = -0.0065    #K/m
        rho = rho0 * (T/(T + a*(h-h0)))**(1. + ((g*M)/(R*a)))
        
    if h >= 11000:
        rho0 = 0.36391 #kg/m^3
        T = 216.65     #K
        h0 = 11000.     #m
        rho = rho0*np.e**((-g*M*(h-h0))/(R*T))
        
    return rho
    
def ISA_temp(h):
    if h < 11000:
        T = 288.15 - 0.0065*h   #in Kelvin
        return T
    if h >= 11000:
        return 216.65    #in Kelvin
    
def Mach(V,h):                  #enter V in m/s
    gamma = 1.4                 #enter h in m
    R = 287 #J/kg/K
    a = np.sqrt(gamma*R*ISA_temp(h))
    M = V/a
    return M 

def Vel(M,h):
    gamma = 1.4                 #enter h in m
    R = 287 #J/kg/K
    a = np.sqrt(gamma*R*ISA_temp(h))
    V = a*M
    return V

#Varying thrust with velocity
def T_alt(T0,h):
    T = T0*(ISA_density(h)/ISA_density(0))**m
    return T

""" Max and min velocities""" 
def Vmax(Wcr,T0,S,A,e,CD0,h):
    #In order to minimise mistakes, split up eq. in more variables
    k = 1. / (np.pi*A*e)
    Tmax = T_alt(T0,h)
    Vmax = (((Tmax/Wcr)*(Wcr/S) + (Wcr/S)*np.sqrt((Tmax/Wcr)**2 - (4*CD0*k)))/(ISA_density(h)*CD0))**0.5
    return Vmax

def Vmin(Wcr,T0,S,A,e,CD0,h,CL_max,Vcr):
    Tmax = T_alt(T0,h)
    #In order to minimise mistakes, split up eq. in more variables
    k = 1. / (np.pi*A*e)
    Vmin = (((Tmax/Wcr)*(Wcr/S) - (Wcr/S)*np.sqrt((Tmax/Wcr)**2 - (4*CD0*k)))/(ISA_density(h)*CD0))**0.5
    Vs = np.sqrt((2*Wcr)/(ISA_density(h)*S*CL_max))
    
    if Vs > Vmin:
        return Vs
    if Vmin > Vs:
        return Vmin


""" Thrust required""" 
def Treq(W,V,S,A,e,h,CD0):
    CL = W / (0.5*ISA_density(h)*(V)**2*S)
    CD = CD0 + (CL**2 / (np.pi*A*e))
    
    Pr = np.sqrt((2*W**3*CD**2)/(ISA_density(h)*S*CL**3))
    
    Tr = Pr/V       #Convert power to thrust by dividing by airspeed
    return Tr   #returns the Power required 



"""Climb performance"""
#Normal rate of climb at a given altitude and airspeed (steady, no acceleration)
def RC(W,T0,V,S,A,e,CD0,h):
    k = 1. / (np.pi*A*e)
    T = T_alt(T0,h)
    RC = V*((T/W) - 0.5*ISA_density(h)*V**2 * (S/W)*CD0 - (W/S)*((2*k)/(ISA_density(h)*V**2)))
    return RC


def RC_unsteady(W,T0,V,S,h,CD):  #at const. EAS thus accelerating during flight
    T = T_alt(T0,h)
    M = Mach(V,h)
    vg_dvdh = 0.5668*M**2            #Constant EAS in tropospere
    D = 0.5*ISA_density(h)*V**2*S*CD
    C = ((T-D)*V)/(W*(1. + vg_dvdh))
    return C


#Maximum climb angle (or steepest climb) and corresponding airspeed and RC 
def steep_climb(T0,W,S,CD0,A,e,h,V):
    T = T_alt(T0,h)
    k = 1. / (np.pi*A*e)
    theta_max = np.arcsin((T/W) - np.sqrt(4*CD0*(k))) 
    V_theta_max = np.sqrt((2/ISA_density(h))*(k/CD0)**0.5 * (W/S)*np.cos(theta_max))
    RC_max_theta = V_theta_max*np.sin(theta_max)
    
    return theta_max*(180/np.pi), V_theta_max,RC_max_theta   #return in degrees


"""Gliding flight: rate of descend"""
def RD(W,S,A,e,CD0,h):
    #Smalles rate of descend is obtained at max. (CL**(3/2))/CD
    k = 1. / (np.pi*A*e)
    CL_CD = 0.25*(3./(k*CD0**(1./3.)))**0.75
    Vv_min = - np.sqrt((2*W)/(ISA_density(h)*S))*(1./CL_CD)

    return Vv_min

def glide_range(L_D,dh):  
    R = (L_D)*dh
    return R


#------------------------------------MAIN PROGRAM---------------------------
H = [0,3048,3048,6096,6096,9144,9144,12192,]
V = [154.33,159.48,174.91,221.21,216.07,221.21,280.37,308.68]

RC_list = []

for i in range(len(H)):
    RCi = RC(Wcr,T0,V[i],S,A,e,CD0,H[i])
    RC_list.append(RCi)
    
plt.plot(H,RC_list)
plt.show()

print RC_list

#Altitude and velocity ranges for the plots
#V = np.arange(50,300,5)              #Criose velocity in m/s
#H = np.arange(3000,,1000)       #Cruise altitude in m
#
#
#plt.subplot(221)
#
#RC_max = [] 
#M_RC_max = []
#V_RC_max = []
#H_RC_max = []
#
#for h in H:
#    RC_list = []
#    M_list = []
#    V_list = []
#    
#    for v in V:
#        RC_list.append(RC(Wcr,T0,v,S,A,e,CD0,h))
#        M_list.append(Mach(v,h))
#        V_list.append(v)
#        
#    plt.plot(M_list,RC_list, label = "%s m" %h)
#    
#    k = RC_list.index(max(RC_list))
#    RC_max.append(max(RC_list))
#    V_RC_max.append(V_list[k])
#    M_RC_max.append(M_list[k])
#    H_RC_max.append(h)
#
#plt.plot(M_RC_max,RC_max," ko",label = " Max. RC" )
#plt.title("Steady rate of climb" )
#plt.xlabel("Mach number")
#plt.ylabel("Rate of climb [m/s]")
#plt.legend()  
#plt.grid(True)
#
#
#plt.subplot(222)
#plt.plot(H,RC_max)
#plt.title("Max. steady rate of climb" )
#plt.xlabel("Altitude [m]" )
#plt.ylabel("Rate of climb [m/s]" )
#plt.grid(True)
#
#plt.subplot(223)
#plt.plot(H,M_RC_max)
#plt.title("Mach number at Max. steady RC" )
#plt.xlabel("Altitude [m]")
#plt.ylabel("Mach number")
#plt.grid(True)
#
#
#plt.show()














