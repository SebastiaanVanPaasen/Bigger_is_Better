# -*- coding: utf-8 -*-
"""
Created on Thu May 09 14:36:48 2019

@author: nikki
"""
#-------------------------MODULES---------------------------
import numpy as np
import math
import scipy as sp
import matplotlib as plt

#--------------------------INPUTS---------------------------- 
"""Geometry"""
#qcsweep    =   quarter chord sweep
#hcsweep    =   half chord sweep
#A          =   aspect ratio
#S          =   Wing surface
#c          =   mean aerodynamic chord
#b          =   wing span
#c_r        =   root chord
#c_flap     =   flap chord
#taper      =   taper ratio
#b_f0       =   end of flap span
#b_fi       =   begin of flap span
#l_f        =   fuselage length
#b_f        =   fuselage diameter
#h_f        =   fuselage height
#l_h        =   tailarm
#A_h        =   aspect ratio horizontal tail

"""Coefficients"""
#Cm_0       =   airfoil dependent moment coefficient
#CL_0       =   CL of flapped wing at angle of attack of 0 deg
#CL_land    =   Cl during landing full flap configuration

"""Others"""
#T0         =   Temperature during landing
#V_app      =   approach speed
#mu_2       =   2D to 3D correction factor from graph
#mu_3       =   correction factor for sweep from graph
#delta_flap =   flap deflection

"""Retrieve from stability cruve program"""
#x_ac       =   aerodynamic centre location over the mac
#Vh_V      =   (V_h/V) ratio
#CL_alpha_AH =  lift gradient of tailless aircraft


#------------------------------VALUES---------------------------
eta         =   0.95                #airfoil efficiency


#-----------------------------DEFINITIONS------------------------
def Sh_S_control(CL_H,CL_AH,l_h,Vh_V,c,Cm_ac,x_ac,x_cg):
    is_list = isinstance(x_cg,list)
    den = (CL_H/CL_AH)*(l_h/c)*Vh_V**2.
    
    if is_list == True:
        Sh_S = []
        for i in range(len(x_cg)):
            Sh_S_i = (1./den)*x_cg + (((Cm_ac/CL_AH)-x_ac)/den)
            Sh_S.append(Sh_S_i)
                    
    if is_list == False:
            Sh_S = (1./den)*x_cg + (((Cm_ac/CL_AH)-x_ac)/den)
        
    return Sh_S
    
def Cm_ac(Cm_ac_w,dflap_Cm_ac,dfus_Cm_ac):
    Cm_ac = Cm_ac_w + dflap_Cm_ac + dfus_Cm_ac 
    return Cm_ac
    
def Cm_ac_w():
    
        
    






