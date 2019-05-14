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
R_design = 1400e03 #[m]



#------------------------------DEFINITIONS-----------------------------------

def CL_CD(A,e,CD0):  #compute CL/CD at cruise
    CL_CD = (3./4.)*np.sqrt((np.pi*A*e)/(3.*CD0))
    return CL_CD
    
    
def Wf_Wto():              #find V for which Wf/W_TO can be minimised
    Wf_Wto_opt = []
    V = []
    for i in range(200,600):
        V.append(i)
        Wf_Wto = 1. - e**((-R_design)/((i/Ct)*(CL_CD(A,e,CD0))))
        Wf_Wto_opt.append(Wf_Wto)
        
    plt.plot(V,Wf_Wto_opt)
    #plt.show()
 
    
def payload_range():
    """Standard weight lines"""
    plt.hlines(MTOW,0.,R_range,"g","--")
    plt.hlines(OEW,0.,R_range,"r","--")
    plt.hlines(MZFW,0.,R_range,"b","--")
    plt.hlines(MLW,0.,R_range,"y","--")    
    
    
    """Compute different ranges"""
    #Finding harmonic range (max. payload, with fuel up to MTOW)
    Wf1 = MTOW-MZFW-W_fr             #fuel weight at max. payload
    R_harmonic = ((V/(g*Ct))*CL_CD(A,e,CD0)*np.log(MTOW/(MTOW-Wf1)))/1000.

    #Max range line (increase fuel, decrease payload)
    R_max = ((V/(g*Ct))*CL_CD(A,e,CD0)*np.log(MTOW/(MTOW-(MFW-W_fr))))/1000.
    
    #Ferry range (no payload all fuel, no MTOW anymore)
    W_TO = OEW + MFW 
    Wf2 =  MFW - W_fr 
    R_ferry = ((V/(g*Ct))*CL_CD(A,e,CD0)*np.log(W_TO/(W_TO - Wf2)))/1000.
    
    
    """Take-off weight line"""    
    #TO line up to MTOW and R_harmonic 
    TO = [MZFW+W_fr,MTOW]                     #additional fuel up to MTOW
    range_TO = [0,R_harmonic]
    
    #TO line until max Range
    TO.append(MTOW)
    range_TO.append(R_max)
    
    #TO line until ferry range
    TO.append(OEW+MFW)
    range_TO.append(R_ferry)
    
    #Closing the diagram
    TO.append(OEW)
    range_TO.append(R_ferry)
    
    plt.plot(range_TO, TO,"purple")
    

    """Range line without reserve fuel"""
    #Max payload line 
    plt.hlines(MZFW,0.,R_harmonic,"k")   
    
    #Line up to Max range, increasing fuel, decreasing payload
    weight = [MZFW, MTOW-MFW]              #weight line at max range
    range_weight = [R_harmonic,R_max]
    
    #Line up to ferry range, all fuel, no payload
    weight.append(OEW)            #weight line up to ferry range
    range_weight.append(R_ferry)

    plt.plot(range_weight,weight, "k")
    
    """Range line including reserve fuel"""
    #Max payload line 
    plt.hlines(MZFW+W_fr,0.,R_harmonic,"k")   
    
    for i in range(len(weight)):
        weight[i] = weight[i] + W_fr

    plt.plot(range_weight, weight, "k")    
    

    plt.xlabel("Range [m]")
    plt.ylabel("Weight [N]")
    plt.show()



payload_range()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    