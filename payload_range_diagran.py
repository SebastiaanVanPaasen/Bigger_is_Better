# -*- coding: utf-8 -*-
"""
Created on Mon May 13 11:45:46 2019

@author: nikki
"""

#-------------------------------MODULES---------------------------------------
import numpy as np
import matplotlib.pyplot as plt

#-------------------------------INPUTS-----------------------------------------
"""Weights"""
#MTOW       =       Maximum take-off weight [N] (optional as fraction of W_TO)
#MLW        =       Maximum landing weight [N] (optional as fraction of W_TO)
#MZFW       =       Maximum zero fuel weight [N] (optional as fraction of W_TO)
#OEW        =       Operational empty weight [N] (optional as fraction of W_TO)
#MWP        =       Maximum payload weight [N] (optional as fraction of W_TO)
#MFW        =       Maximum fuel weight [N] (optional as fraction of W_TO)

"""Aicraft configuration"""
#CD0        =       Zero lift drag coefficient
#A          =       Aspect ratio
#e          =       Oswald efficiency factor 
#R_des      =       Design range [m]


#------------------------STATISTICAL INPUTS----------------------------------
#W_fr        = 0.05*W_f      #Reserve fuel weight as percentage of fuel weight
Ct           = 12e-06      #Specific fuel conspumtion [kg/N/s] from B737MAX + requirement 


#------------------------------DEFINITIONS-----------------------------------

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

R_range = 20000.  #range of x-axis


def CL_CD(A,e,CD0):  #compute CL/CD at cruise
    CL_CD = (3./4.)*np.sqrt((np.pi*A*e)/(3.*CD0))
    return CL_CD
    
def plot():
    #Finding harmonic range (max. payload, with fuel up to MTOW)
    Wf1 = MTOW-MZFW-W_fr    #fuel weight at max. payload
    R_harmonic = ((V/(g*Ct))*CL_CD(A,e,CD0)*np.log(MTOW/(MTOW-Wf1)))/1000.

    #Max payload line 
    plt.hlines(MZFW,0.,R_harmonic,"k")
    plt.hlines(MZFW+W_fr,0.,R_harmonic,"k")    
    
    #fuel line up to MTOW and R_harmonic (y = ax + b)
    a1 = (MTOW - (MZFW+W_fr))/R_harmonic        
    b1 = MZFW+W_fr
    fuel = [MZFW+W_fr,a1*R_harmonic+b1]     #additional fuel up to MTOW
    range_fuel = [0,R_harmonic]
    plt.plot(range_fuel,fuel)               #plot additional fuel line
    
    #Max range line (increase fuel, decrease payload)
    R_max = ((V/(g*Ct))*CL_CD(A,e,CD0)*np.log(MTOW/(MTOW-MFW)))/1000.
    
    a2 = ((MTOW-MFW) - (MZFW))/(R_max - R_harmonic)
    b2 = MZFW
    weight = [MZFW, a2*R_max+b2]
    range_weight = [R_harmonic, R_max]
    
    #Ferry range (no payload all fuel, no MTOW anymore)
    W_TO = OEW + MFW + W_fr
    Wf2 =  MFW
    R_ferry = ((V/(g*Ct))*CL_CD(A,e,CD0)*np.log(W_TO/(W_TO - Wf2)))/1000.

    a3 = (OEW -(weight[-1]))/(R_ferry-R_max)
    b3 = weight[-1]
    weight.append(a3*R_ferry+b3)
    range_weight.append(R_ferry)
    
    
    plt.plot(range_weight,weight)
    
    
    #Standard weight lines
    plt.hlines(MTOW,0.,R_range,"g","--")
    plt.hlines(OEW,0.,R_range,"r","--")
    plt.hlines(MZFW,0.,R_range,"b","--")
    plt.hlines(MLW,0.,R_range,"y","--")
    
    #plt.ylim(0,11000)
    plt.xlabel("Range [m]")
    plt.ylabel("Weight [N]")
    plt.show()

    
    
    
    
print plot()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
