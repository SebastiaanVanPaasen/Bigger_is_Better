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
#qcsweep    =   quarter chord sweep [rad]
#hcsweep    =   half chord sweep [rad]
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
#S_net      =   area in fuselage

"""Coefficients"""
#Cm_0       =   airfoil dependent moment coefficient [1/rad]
#CL_0       =   CL of flapped wing at angle of attack of 0 deg [1/rad]
#CL_land    =   Cl during landing full flap configuration [1/rad]

"""Others"""
#T0         =   Temperature during landing  [K]
#V_app      =   approach speed [m/s]
#delta_flap =   flap deflection [rad]
#alpha_land =   angle of attach during approach [rad]
#M_app      =   approach speed in Mach
#eta        =   airfoil efficiency = 0.95

"""Retrieve from stability cruve program"""
#x_ac       =   aerodynamic centre location over the mac [n/m]
#Vh_V      =   (V_h/V) ratio  [m/s / m/s]
#CL_alpha_AH =  lift gradient of tailless aircraft [1/rad]


#------------------------------VALUES GRAPHS--------------------------
#mu_2       =   2D to 3D correction factor from graph  
#mu_3       =   correction factor for sweep from graph
#dc_c_f     =   flap geometry ratio, see torenbeek book


#-----------------------------DEFINITIONS------------------------
def Sh_S_control(CL_H,CL_AH,l_h,Vh_V,c,Cm_ac,x_ac,x_cg):                        #final definition for control curve
    is_list = isinstance(x_cg,list)                                             #as a function of the cg position over mac
    den = (CL_H/CL_AH)*(l_h/c)*Vh_V**2.     #denominator                        #should be ebaluated at landing conditions
    
    if is_list == True:
        Sh_S = []
        for i in range(len(x_cg)):
            Sh_S_i = (1./den)*x_cg + (((Cm_ac/CL_AH)-x_ac)/den)
            Sh_S.append(Sh_S_i)
                    
    if is_list == False:
            Sh_S = (1./den)*x_cg + (((Cm_ac/CL_AH)-x_ac)/den)
        
    return Sh_S
    
def Cm_ac(Cm_ac_w,dflap_Cm_ac,dfus_Cm_ac):                                      #C_mac consists of three parts
    Cm_ac = Cm_ac_w + dflap_Cm_ac + dfus_Cm_ac                                  #due to wing, due to fuselage and due to flaps
    return Cm_ac
    
def Cm_ac_w(Cm_0,A,qcsweep):                                                    #Cm_ac change due to wing
    Cm_ac_w = Cm_0*((A*np.cos(qcsweep)**2.)/(A + 2.*np.cos(qcsweep)))
    return Cm_ac_w
    
def dfus_Cm_ac(b_f,h_f,l_f,S,S_net,A,hcsweep,M_app,eta,CL_0):                   #Cm_ac change due to fuselage
    beta = np.sqrt(1. - M_app**2)
    CL_alpha_w = (2.*np.pi*A)/(2 + np.sqrt(4. + ((A*beta)/eta)**2 + \\
    (1. + (np.tan(hcsweep)**2. / beta**2.))))
    CL_alpha_AH = CL_alpha_w * (1. + 2.15*(b_f/b)) * (S_net/S) + \\
    (np.pi/2.)*(b_f**2./S)
    
    dfus_Cm_ac = -1.8*(1. - (2.5*b_f/b))*((np.pi*b_f*h_f*l_f)/4.*S*c)* \\
    (CL_0/CL_alpha_AH)
    return dfus_Cm_ac

def dflap_Cm_ac():                                                              #Cm_ac due to flaps
    dflap_Cl_max = 2.*delta_flap * np.sin(delta_flap)                           #change in airfoil Cl_max due to flaps
                                                                                #calculated using thin airfoil theory as initial guess
    Swf_S = ((b_f0-b_fi)/b)*(1. + ((1.-taper)/(1.+taper))*(1.-((b_f0+b_fi)/b))) #ratio flap affected area over total are
    theta_f = np.cos(2*(c_f/c) - 1.)**(-1.)                                     #angle on flap geometry depending on chord
    mu_1 = 0.5*(1.-(c_f/c))*((np.sin(theta_f))/(np.pi - \\                      #correction factor for camber effect 3D wing
    (theta_f-np.sin(theta_f)))) 
    
    c_dash_c = 1. + dc_c_f*(c_f/c)                                              #flap chord over wing chord
    
    dflap_Cm_qc = mu_2*(-mu_1*dflap_Cl_max*c_dash_c - \\                        #Cm_quarter chord change due to flaps
    (CL_land + flap_Cl_max*(1-Swf_S))*(1./8.)*c_dash_c*(c_dash_c-1.)) + \\
    0.7*(A/(1.+2./A))*mu_3*dflap_Cl_max*np.tan(qcsweep)
    
    dflap_Cm_ac = dflap_Cm_qc - CL_land*(0.25 - x_ac)                           #Transform to change wrt to aerodyanmic centre
    return dflap_Cm_ac
    
    
    
        
    






