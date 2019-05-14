# -*- coding: utf-8 -*-
"""
Created on Tue May 14 08:55:23 2019

@author: nikki
"""

#-------------------------------MODULES---------------------------------------
import numpy as np
import matplotlib.pyplot as plt

#-------------------------------INPUTS-----------------------------------------
"""Weights"""
#Wcr        =       Aircraft weight in N during Cruise

"""Aicraft configuration"""
#CD0        =       Zero lift drag coefficient
#A          =       Aspect ratio
#e          =       Oswald efficiency factor 
#R_des      =       Design range [m]
#S          =       Surface area m^2

#------------------------STATISTICAL INPUTS----------------------------------

Ct           = 12e-06      #Specific fuel conspumtion [kg/N/s] from B737MAX 


#------------------------------VERIFICATION DATA--------------------------------

"""Inputs unit test based on B737 MAX-8"""

MTOW = 82190.*9.81
OEW  = 45065*9.81
MLW  = 69308.*9.81
MZFW = 65952.* 9.81
MFW  = 20826.*9.81          # Maximum fuel weight (including reserve fuel)
W_fr = MFW/105. * 5.        #reserve fuel

A = 8.45
e = 0.85
CD0 = 0.020
V = 236. #m/s
g = 9.81
S = 124.5 

R_range = 11000.  #range of x-axis
R_des = 7000 #[km]
Wcr = MLW #assumption for now
pax_max = 450
n = 0.8 #load factor of numbe rof passengers

#-----------------------------DEFINITIONS-------------------------------------
#Standard air range (SAR) = distances travelled per mass fuel burned
#Fuel burn is related to speed, altitude, thrust (or drag in steady flight)
#unit of SAR (m/s)/(kg/s)

def ISA_density(h):      # enter height in m
    M = 0.0289644       #kg/mol molar mass of Earth's air
    R = 8.3144590       #universal gas constant Nm/MolK
    
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


def SAR(h,A,S,e,CD0,Ct,Wcr):                 #enter V in m/s
    V = np.linspace(400,1100,100)
    SAR = []
    
    for v in V:
        k = 1./(np.pi*A*e)
        q = 0.5*ISA_density(h)*(v/3.6)**2
        SARi = 1./((v/3.6) / ( (CD0 + k *(Wcr/(q*S))**2) *q*S*Ct )) #in kg/m
        SAR.append(SARi*1000.)   #in kg/km
        
    return SAR,V
    

    
#------------------------------MAIN PROGRAM------------------------------------

dh = 500                #step size in altitude
H = range(7000,12000,dh)#altitude range


min_SAR = []
V_minSAR = []


#For a given altitude (in def) run it for different speeds
for h in H:   
    SAR_list = SAR(h,A,S,e,CD0,Ct,Wcr)[0]
    V = SAR(h,A,S,e,CD0,Ct,Wcr)[1]

    min_SAR.append(min(SAR_list))
    i = SAR_list.index(min(SAR_list))
    V_minSAR.append(V[i])    
    
    
    plt.subplot(211)
    plt.plot(V,SAR_list,label='%s altitude [m]' % h)
    plt.title('Fuel consumption w.r.t. airspeed')
    plt.xlabel("Airspeed [km/s]")
    plt.ylabel("Fuel consumption [kg/km]")

plt.legend()


for j in range(len(min_SAR)):
    plt.subplot(212)
    plt.xlabel("Airspeed at minimum SAR [km/s]")
    plt.ylabel("Minimum Fuel consumption [kg/km]")
    plt.plot(V_minSAR[j],min_SAR[j],'o', label = '%s altitude [m]' % H[j])
    plt.title('Minimum fuel consumption with corresponding airspeed and altitude')
  
plt.legend()
plt.show()


#Fuel consumed per km per passenger/seat
for h in H:   
    pax = pax_max*n
    SAR_list = (SAR(h,A,S,e,CD0,Ct,Wcr)[0])
    V = SAR(h,A,S,e,CD0,Ct,Wcr)[1]

    for i in range(len(SAR_list)):
        SAR_list[i] = SAR_list[i]/pax 
        
    plt.plot(V,SAR_list,label='%s altitude [m]' % h)
    plt.title('Fuel consumption per passenger w.r.t. airspeed')
    plt.xlabel("Airspeed [km/s]")
    plt.ylabel("Fuel consumption [kg/km/passenger]")

plt.legend()
plt.show()







