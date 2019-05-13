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
    V = np.linspace(200,1000,800)   #V range in km/h
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
    


        
        
    



