# -*- coding: utf-8 -*-
"""
Created on Mon May 13 15:34:38 2019

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