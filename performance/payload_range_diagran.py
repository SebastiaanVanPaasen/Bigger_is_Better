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
#W_fr        = 0.05*W_f      #Reserve fuel weight as percentage of fuel weight

"""Aicraft configuration"""
#CD0        =       Zero lift drag coefficient
#A          =       Aspect ratio
#e          =       Oswald efficiency factor 
#R_des      =       Design range [m]


#------------------------STATISTICAL INPUTS----------------------------------

    #Specific fuel conspumtion [kg/N/s] from B737MAX 
  
#------------------------------VERIFICATION DATA--------------------------------

"""Inputs unit test based on B737 MAX-8"""
"""STRW"""
#MTOW = 1712063.74
#OEW  =  900623.75
#
##MLW  = MTOW - 0.3*OEW
#MZFW = OEW + 508587.17
#MFW  =  297908.01*1.4     # Maximum fuel weight (including reserve fuel)
#W_fr = MFW/105. * 5.        #reserve fuel
#A = 18.
#V = 206.53 #m/s
#CD0 = 0.0276
#g = 9.81
#e = 0.9
#R_range = 2500.  #range of x-axis
#R_des = 1400. #[km]
#Ct0 = 55e-06



"""HPBE"""
#MTOW = 1835973.64
#OEW  =  959038.84
#
##MLW  = MTOW - 0.3*OEW
#MZFW = OEW + 508587.17
#MFW  =  363032.71*1.3     # Maximum fuel weight (including reserve fuel)
#W_fr = MFW/105. * 5.        #reserve fuel
#A = 8.
#V = 209.6 #m/s
#CD0 = 0.0219
#g = 9.81
#e = 0.9
#R_range = 3200.  #range of x-axis
#R_des = 1600. #[km]
#Ct0 = 55e-06


"""DD4E"""
#MTOW = 1533186.85
#OEW  =  741915.71
#
##MLW  = MTOW - 0.3*OEW
#MZFW = OEW + 508587.17
#MFW  =  277919.32*1.4     # Maximum fuel weight (including reserve fuel)
#W_fr = MFW/105. * 5.        #reserve fuel
#A = 8.5
#V = 224.58 #m/s
#CD0 = 0.0205
#g = 9.81
#e = 0.9
#R_range = 3200.  #range of x-axis
#R_des = 1600. #[km]
#Ct0 = 55e-06



"""DD2E"""
#MTOW = 1545876.62
#OEW  =  740611.66
#
##MLW  = MTOW - 0.3*OEW
#MZFW = OEW + 508587.17
#MFW  =  292191.26*1.4     # Maximum fuel weight (including reserve fuel)
#W_fr = MFW/105. * 5.        #reserve fuel
#A = 9.
#V = 215.62 #m/s
#CD0 = 0.0202
#g = 9.81
#e = 0.9
#R_range = 3200.  #range of x-axis
#R_des = 1600. #[km]
#Ct0 = 55e-06


"""Aerodynamic design"""
#MTOW = 1828542.22
#OEW  =  1000102.71
##MLW  = MTOW - 0.3*OEW
#MZFW = OEW + 508587.17
#MFW  =  313967.31*1.4     # Maximum fuel weight (including reserve fuel)
#W_fr = MFW/105. * 5.        #reserve fuel
#A = 14.
#V = 221.281 #m/s
#CD0 = 0.0262
#g = 9.81
#e = 0.9
#R_range = 3000.  #range of x-axis
#R_des = 1600. #[km]
#Ct0           = 55e-06  

"""Ref aircraft B737-Max 8 """
MTOW = 82190*9.81
OEW = 45065*9.81
MLW = 69308*9.81
MZFW = 65952*9.81
MFW = 20826*9.81
W_fr = MFW/105 * 5.
A = 8.45
e = 0.85
CD0 = 0.020
V = 236.
g = 9.81

R_range = 11000
R_des = 6590.
Ct0 = 10.7e-06


#------------------------------DEFINITIONS-----------------------------------

def CL_CD(A,e,CD0):  #compute CL/CD at cruise
    CL_CD = (3./4.)*np.sqrt((np.pi*A*e)/(3.*CD0))
    return CL_CD
    
    
def Wf_Wto():              #find V for which Wf/W_TO can be minimised
    Wf_Wto_opt = []
    V = []
    for i in range(200,600):
        V.append(i)
        Wf_Wto = 1. - e**((-R_des*1000)/((i/Ct0)*(CL_CD(A,e,CD0))))
        Wf_Wto_opt.append(Wf_Wto)
        
    plt.plot(V,Wf_Wto_opt)
    #plt.show()
 
    
def payload_range():
    
    """Compute different ranges"""
    #Finding harmonic range (max. payload, with fuel up to MTOW)
    Wf1 =  (MTOW-MZFW-W_fr)             #fuel weight at max. payload
    Ct = (Ct0 / 233.083) * (V) 
    R_harmonic = ((V/(g*Ct))*CL_CD(A,e,CD0)*np.log(MTOW/(MTOW-Wf1)))/1000.

    #Max range line (increase fuel, decrease payload)
    R_max = ((V/(g*Ct))*CL_CD(A,e,CD0)*np.log(MTOW/(MTOW-((MFW-W_fr)*1))))/1000.

    #Ferry range (no payload all fuel, no MTOW anymore)
    W_TO = OEW + MFW 
    Wf2 =  (MFW - W_fr )#*0.4
    R_ferry = ((V/(g*Ct))*CL_CD(A,e,CD0)*np.log(W_TO/(W_TO - Wf2)))/1000.
    
    plt.vlines(R_harmonic,OEW,MTOW,"m","--",label="Harmonic Range")
    plt.vlines(R_max, OEW,MTOW,"crimson","--",label="Maximum Range")
    plt.vlines(R_ferry,OEW,MTOW,"gray","--",label="Ferry Range")
    #plt.vlines(R_des,OEW,MTOW,"orange","--",label="Design Range")
   
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
    
    plt.plot(range_TO, TO,"navy",label = "TO weight")
    

    """Range line without reserve fuel"""
    #Max payload line 
    plt.hlines(MZFW,0.,R_harmonic,"k")   
    
    #Line up to Max range, increasing fuel, decreasing payload
    weight = [MZFW, MTOW-MFW]              #weight line at max range
    range_weight = [R_harmonic,R_max]
    
    #Line up to ferry range, all fuel, no payload
    weight.append(OEW)            #weight line up to ferry range
    range_weight.append(R_ferry)

    plt.plot(range_weight,weight, "k", label = "Payload weight")
    
    """Range line including reserve fuel"""
    #Max payload line 
    plt.hlines(MZFW+W_fr,0.,R_harmonic,"gray")   
    
    for i in range(len(weight)):
        weight[i] = weight[i] + W_fr
    plt.figure(1)
    plt.plot(range_weight, weight, "gray",label = "With reserve fuel")    
    
    """Standard weight lines"""
    plt.hlines(MTOW,0.,R_ferry,"g","--",label = "MTOW")
    plt.hlines(OEW,0.,R_ferry,"r","--",label = "OEW")
    plt.hlines(MZFW,0.,R_ferry,"b","--", label = "MZFW")
    #plt.hlines(MLW,0.,R_range,"y","--", label = "MLW")    
    
    #plt.title('Payload - range diagram')
    plt.xlabel("Range [km]")
    plt.ylabel("Weight [N]")
    plt.ylim(0.95*OEW,1.05*MTOW)
    plt.legend()
    plt.show()
    return R_harmonic,R_max, R_ferry

#a = (payload_range()[0][1]-p6ayload())
print (payload_range())
    
    
    
    
    
    
    
    
    
    
    
    
