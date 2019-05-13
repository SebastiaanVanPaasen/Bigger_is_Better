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

"""Aicraft configuration"""
#CD0        =       Zero lift drag coefficient
#A          =       Aspect ratio
#e          =       Oswald efficiency factor 


#------------------------STATISTICAL INPUTS----------------------------------
#W_fr        = 0.05*W_f      #Reserve fuel weight as percentage of fuel weight
C_t          = 1.2e-05*0.9   #Specific fuel conspumtion [kg/Ns] from B737MAX + requirement 


#------------------------------DEFINITIONS-----------------------------------

def CL_CD(A,e,CD0):  #compute CL/CD at cruise
    CL_CD = (3./4.)*np.sqrt((np.pi*A*e)/(3.*CD0))
    return CL_CD
    
def plot():
    