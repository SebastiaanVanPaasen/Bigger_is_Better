# -*- coding: utf-8 -*-
"""
Created on Mon May 13 16:00:12 2019

@author: nikki
"""
#-------------------------------MODULES---------------------------------------
import numpy as np
import matplotlib.pyplot as plt

#-------------------------------INPUTS-----------------------------------------
"""Geometry"""
# S         =       Surface area 


"""Aicraft configuration"""
#CD0        =       Zero lift drag coefficient
#A          =       Aspect ratio
#e          =       Oswald efficiency factor 
#R_des      =       Design range [m]
#Tcr        =       Thrust setting during cruise
#T_W        =       Thrust to weight ratio during cruise


#------------------------------VERIFICATION DATA--------------------------------

"""Inputs unit test based on B737 MAX-8"""

MTOW = 82190.*9.81
OEW  = 45065*9.81
MLW  = 69308.*9.81
MZFW = 65952.* 9.81
MFW  = 20826.*9.81          # Maximum fuel weight (including reserve fuel)
W_fr = MFW/105. * 5.        #reserve fuel

Tcr = 117.3e03  #N
A = 8.45
e = 0.85
CD0 = 0.020
V = 236. #m/s
g = 9.81
S = 124.5 
T_W = 0.3027

h = 10000
#assume weight in cruise is at least MLW to allow for a save landing
Wcr = MLW

#------------------------------DEFINITIONS-----------------------------------

def ISA_density(h):      # enter height in m
    # Temperature
    if h<=11000:
        T0=288.15
        a=-0.0065
        T=T0+a*h
    if 11000<h<=20000:
        T1=216.65

    #Density
    if h<=11000:
        d0=1.225
        g=9.80665
        R=287.00
        D=d0*(T/T0)**((-g/(a*R))-1)
        return D
        
    if 11000<h<=20000:
        d0=1.225
        g=9.80665
        R=287.
        T1=216.65
        D1=d0*np.e**(((-g/(R*T1))*(h)))
        return D1

def drag_plot(h,S,A,e,Wcr):         #Drag VS velocity graph
    V = np.linspace(200,1200,800)   #V range in km/h
    CL_list = []
    CD_list = []
    for i in V:
        CL = Wcr / (0.5*ISA_density(h)*(i/3.6)**2*S)
        CD = CD0 + (CL**2 / (np.pi*A*e))
        CL_list.append(CL)
        CD_list.append(CD)
    
    D = []
    for j in range(len(CL_list)):
        D.append((Wcr * (CD_list[j]/CL_list[j]))/1000.)
        
    plt.plot(V, D)
    plt.xlabel("Speed [km/h]")
    plt.ylabel("Drag [kN]")
    plt.show()
    
    k = D.index(min(D))            
    return V[k]                 #speed at which the minimum drag is experienced


def power_plot(h,Wcr,S,Tcr):
    V = np.linspace(200,1200,800)   #V range in km/h
    Pr_list = []                    #power required list in kW
    Pa_list = []                    #Power available list in kW
    for i in V:
        CL = Wcr / (0.5*ISA_density(h)*(i/3.6)**2*S)
        CD = CD0 + (CL**2 / (np.pi*A*e))
        Pr = Wcr*np.sqrt((2*Wcr*CD**2)/(S*ISA_density(h)*CL**3))
        Pr_list.append(Pr/1000.)
        
        Pa_list.append((i/3.6)*Tcr/1000.)
        
    plt.plot(V,Pr_list,label = "P_req")
    plt.plot(V,Pa_list, label = "P_av")
    plt.legend(loc = "upper right")
    plt.xlabel("Airspeed [km/h]")
    plt.ylabel("Power [kW]")
    plt.show()
    
    Pc_list = []                    # Excess power used to climb in kW
    for i in range(len(Pa_list)):
        Pc = Pa_list[i] - Pr_list[i]
        Pc_list.append(Pc)
    
    RC_list = []                    # rate of climb is proportional to excess power
    for j in Pc_list:
        RC_list.append((j*1000)/Wcr)
    
    plt.plot(V,RC_list)
    plt.xlabel("Airspeed [km/h]")   
    plt.ylabel("Rate of climb [m/s]")
    plt.show()
    
    k = RC_list.index(max(RC_list))
    return V[k], max(RC_list)        #speed at which the maxmimum RC is obtained 


#----------------------------MAIN PROGRAM-----------------------------------
"""Find RC max and corresponding airspeed + fuel consumption at different h"""
h = np.linspace(7000,12000,5000)  #look at altitudes between 7 and 12 km

V_RCmax = []
RC_max = []
for i in h:
    V_RCmaxi = power_plot(i,Wcr,S,Tcr)[0]
    RC_maxi = power_plot(i,Wcr,S,Tcr)[1]
    
    V_RCmax.append(V_RCmaxi)
    RC_max.append(RC_maxi)
    


    
    
    

        
        
    



