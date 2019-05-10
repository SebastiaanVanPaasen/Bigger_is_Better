# -*- coding: utf-8 -*-
"""
Created on Thu May 09 14:36:48 2019

@author: nikki
"""
#-------------------------MODULES---------------------------
import numpy as np
import math
import scipy as sp
import matplotlib.pyplot as plt

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
#x_ac        =   aerodynamic centre location over the mac [n/m]
#Vh_V        =   (V_h/V) ratio  [m/s / m/s]
#CL_AH       =   lift coefficient A-H during landing
#CL_H        =   liftcoefficient tail during landing


#------------------------------VALUES FROM GRAPHS--------------------------
#mu_2       =   2D to 3D correction factor from graph (Torenbeek)
#mu_3       =   correction factor for sweep from graph (Torenbeek)
#dc_c_f     =   flap geometry ratio, see torenbeek book


#------------------------------------DEFINITIONS-----------------------------
def Sh_S_control(CL_H,CL_AH,l_h,Vh_V,x_cg,Cm_0,A,qcsweep,delta_flap,b,b_f0,b_fi,taper,c,c_f,dc_c_f,mu_2,mu_3,x_ac,CL_land,b_f,h_f,l_f,S,S_net,hcsweep,M_app,eta,CL_0):                        #final definition for control curve
                                                                                #as a function of the cg position over mac
    den = (CL_H/CL_AH)*(l_h/c)*Vh_V     #denominator                            #should be ebaluated at landing conditions
    
    Sh_S = []
    for i in range(len(x_cg)):
        Sh_S_i = (1./den)*x_cg[i]  + (((Cm_ac(Cm_0,A,qcsweep,delta_flap,b,b_f0,b_fi,taper,c,c_f,dc_c_f,mu_2,mu_3,x_ac,CL_land,b_f,h_f,l_f,S,S_net,hcsweep,M_app,eta,CL_0)/CL_AH)-x_ac)/den)
        Sh_S.append(Sh_S_i)
        
    return Sh_S
    
    
    
def Cm_ac(Cm_0,A,qcsweep,delta_flap,b,b_f0,b_fi,taper,c,c_f,dc_c_f,mu_2,mu_3,x_ac,CL_land,b_f,h_f,l_f,S,S_net,hcsweep,M_app,eta,CL_0):                                      #C_mac consists of three parts
    Cm_ac = Cm_ac_w(Cm_0,A,qcsweep) + dflap_Cm_ac(delta_flap,b,b_f0,b_fi,taper,c,c_f,dc_c_f,mu_2,mu_3,A,qcsweep,x_ac,CL_land) + dfus_Cm_ac(b_f,h_f,l_f,b,c,S,S_net,A,hcsweep,M_app,eta,CL_0)                #due to wing, due to fuselage and due to flaps
    return Cm_ac
    
    
    
def Cm_ac_w(Cm_0,A,qcsweep):                                                    #Cm_ac change due to wing 
    Cm_ac_w = Cm_0*((A*np.cos(qcsweep)**2.)/(A + 2.*np.cos(qcsweep)))
    return Cm_ac_w
    
    
    
def dfus_Cm_ac(b_f,h_f,l_f,b,c,S,S_net,A,hcsweep,M_app,eta,CL_0):               #Cm_ac change due to fuselage
    beta = np.sqrt(1. - M_app**2)
    CL_alpha_w = (2.*np.pi*A)/(2. + np.sqrt(4. + ((A*beta)/eta)**2 +  (1. + (np.tan(hcsweep)**2. / beta**2.))))
    CL_alpha_AH = CL_alpha_w * (1. + 2.15*(b_f/b)) * (S_net/S) + (np.pi/2.)*(b_f**2./S)
    
    dfus_Cm_ac = -1.8*(1. - (2.5*b_f/b))*((np.pi*b_f*h_f*l_f)/(4.*S*c))*(CL_0/CL_alpha_AH)
    return dfus_Cm_ac



def dflap_Cm_ac(delta_flap,b,b_f0,b_fi,taper,c,c_f,dc_c_f,mu_2,mu_3,A,qcsweep,x_ac,CL_land):
                                                                                #Cm_ac due to flaps
    dflap_Cl_max = 2.*delta_flap * np.sin(delta_flap)                           #change in airfoil Cl_max due to flaps
                                                                                #calculated using thin airfoil theory as initial guess
    Swf_S = ((b_f0-b_fi)/b)*(1. + ((1.-taper)/(1.+taper))*(1.-((b_f0+b_fi)/b))) #ratio flap affected area over total are
    theta_f = np.cos(2*(c_f/c) - 1.)**(-1.)                                     #angle on flap geometry depending on chord
    mu_1 = 0.5*(1.-(c_f/c))*((np.sin(theta_f))/(np.pi -(theta_f-np.sin(theta_f))))#correction factor for camber effect 3D wing
                                                                                
    c_dash_c = 1. + dc_c_f*(c_f/c)                                              #flap chord over wing chord
                                                                                
    dflap_Cm_qc = mu_2*(-mu_1*dflap_Cl_max*c_dash_c -(CL_land + dflap_Cl_max*(1-Swf_S))*(1./8.)*c_dash_c*(c_dash_c-1.)) +  0.7*(A/(1.+ (2./A)))*mu_3*dflap_Cl_max*np.tan(qcsweep)
                                                                                #Cm_quarter chord change due to flaps
    dflap_Cm_ac = dflap_Cm_qc + CL_land*(0.25 - x_ac)                           #Transform to change wrt to aerodyanmic centre
    return dflap_Cm_ac
    
  

#-----------------------------Unit tests----------------------------------
"""Geometry"""
qcsweep    =  0.5515 #quarter chord sweep [rad]
hcsweep    =  0.4869 #half chord sweep [rad]
A          =  8.67 #aspect ratio
S          =  427.8 #Wing surface
c          =  8.75 #mean aerodynamic chord
b          =  60.9 #wing span
c_r        =   12. #root chord
c_f        =   2.625 #flap chord
taper      =   0.149 #taper ratio
b_f0       =   15.6 #end of flap span
b_fi       =   4.1 #begin of flap span
l_f        =   56.4 #fuselage length
b_f        =   6.2 #fuselage diameter
h_f        =   6.2 #fuselage height
l_h        =   24.71 #tailarm
A_h        =   4.5 #aspect ratio horizontal tail
S_net      =   75. #area in fuselage

"""Coefficients"""
Cm_0       =   -0.1 #airfoil dependent moment coefficient 
CL_0       =   0.5 #CL of flapped wing at angle of attack of 0 deg 
CL_land    =   2.45 #CL during landing full flap configuration

"""Others"""
T0         =   293. #Temperature during landing  [K]
V_app      =   75.   #approach speed [m/s]
delta_flap =  0.5235 #flap deflection [rad]
alpha_land =   0.1745 #angle of attach during approach [rad]
M_app      =   0.22 #approach speed in Mach
eta        =  0.95# airfoil efficiency = 0.95

"""Retrieve from stability cruve program"""
x_ac       =  0.3   # aerodynamic centre location over the mac [n/m]
Vh_V      =  0.8    #(V_h/V) ratio  [m/s / m/s]
CL_alpha_AH =   4.788 # lif#t gradient of tailless aircraft [1/rad]
CL_AH       =  3.25 #lift coefficient A-H
CL_H        = -0.8 #liftcoefficient tail



#-----------------------------------VALUES GRAPHS------------------------------
mu_2       = 1. #2D to 3D correction factor from graph  
mu_3       = 0.025#  correction factor for sweep from graph
dc_c_f     =  0.5 # flap geometry ratio, see torenbeek book

x_cg = np.linspace(0.,1,100)

Sh_S = Sh_S_control(CL_H,CL_AH,l_h,Vh_V,x_cg,Cm_0,A,qcsweep,delta_flap,b,b_f0,b_fi,taper,c,c_f,dc_c_f,mu_2,mu_3,x_ac,CL_land,b_f,h_f,l_f,S,S_net,hcsweep,M_app,eta,CL_0)
x_as = 0.*x_cg

#plt.plot(x_cg,Sh_S,x_cg,x_as,"k")
#plt.ylim(0.,1.2)
#plt.show()







