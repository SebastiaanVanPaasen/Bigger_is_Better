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

R_range = 11000.  #range of x-axis
R_des = 7000 #[km]
Wcr = MLW #assumption for now

#-----------------------------DEFINITIONS-------------------------------------
#Standard air range (SAR) = distances travelled per mass fuel burned
#Fuel burn is related to speed, altitude, thrust (or drag in steady flight)

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


def SAR(V,h,A,S,e,CD0,Ct):                 #enter V in m/s
    k = 1./(np.pi*A*e)
    q = 0.5*ISA_density(h)*V**2
    SAR = V / ( (CD0 + k *(Wcr/(q*S))**2) *q*S*Ct )
    return SAR
    
#------------------------------MAIN PROGRAM------------------------------------
    
    

















