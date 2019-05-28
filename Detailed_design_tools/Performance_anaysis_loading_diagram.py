# -*- coding: utf-8 -*-
"""
Created on Mon May 27 17:24:10 2019

@author: nikki
"""


"""Assumed that it is a jet aircraft: otherwise change equations"""
#----------------------------IMPORT MODULES-----------------------------------
import numpy as np


#------------------------------INPUTS------------------------------------------
"""Input values for the B737-800 aircraft"""
MTOW       =       78220*9.81 #Maximum take-off weight [N]
W_TO       =       MTOW       #Weight at take-off [N]
T_TO       =       107*1000   #Total static thrust of all engines at take-off [N]
T          =       91.63*1000 #Take-off thrust
CL_maxto   =        2.5
 
V          =       73.05      #Cruise velocity TAS [m/s]
#M          =       0.788      #Cruise Mach number [-]
A          =       9.44       #Aspect ratio [-]
e          =       0.85       #Oswald efficiency factor [-]

t_c        =       0.125      #Thickness over chord ration airfoil
qcsweep    =       0.4363     #Quater chord sweep in [rad]
S          =       124.60     #Surface area wing [m^2]
b_f        =       3.73       #Width of fuselage [m]
h_f        =       3.73       #Height of fuselage [m]
l_f        =       38.08      #Length of fuselage [m]

psi_TO     =       342.06     #Specific Thrust N/airflow [N/kg/s?]
bypass     =       5.5        #Bypass ratio of the engine

g          =        9.80665
#---------------------------DEFINITIONS----------------------------------------
"""ISA definitions"""

def ISA_density(h):             # enter height in m
    M = 0.0289644               #kg/mol molar mass of Earth's air
    R = 8.3144590               #universal gas constant Nm/MolK
    g = 9.80665                 #gravitationa constant 
    
    if h < 11000:
        rho0 = 1.225            #kg/m^3
        T = 288.15              #K
        h0 = 0.                 #m
        a = -0.0065             #K/m
        rho = rho0 * (T/(T + a*(h-h0)))**(1. + ((g*M)/(R*a)))
        
    if h >= 11000:
        rho0 = 0.36391          #kg/m^3
        T = 216.65              #K
        h0 = 11000.             #m
        rho = rho0*np.e**((-g*M*(h-h0))/(R*T))
        
    return rho
    
    
def ISA_temp(h):
    if h < 11000:
        T = 288.15 - 0.0065*h   #in Kelvin
        return T
    if h >= 11000:
        return 216.65           #in Kelvin

def ISA_press(h):
    R = 8.3144590               #universal gas constant Nm/MolK
    g = 9.80665                 #gravitationa constant 
    p0 = 101325.0
    T0 = 288.15 #[K]
    a = -0.0065 
    T = ISA_temp(h)
    
    if h <= 11000:
        p = p0*(T/T0)**(-g/(a*R))
    else:
        p = p0*np.exp**(((-g/(R*T))*(h)))
        
    return p
         
def Mach(V,h):                  #enter V in m/s
    gamma = 1.4                 #enter h in m
    R = 287 #J/kg/K
    a = np.sqrt(gamma*R*ISA_temp(h))
    M = (V/a)
    return M 


"""Take-off T/W jet"""
def TW_TO_jet(T,T_TO,MTOW,W_TO,M,t_c,qcsweep,S,A,l_f,b_f,h_f,psi_TO,bypass,e,rho,p0): #T/W for TO of a jet aircraft

    #r_uc       =       Correction factor for the undercarriage
    #           =       1.0 fully retracted, no fairings
    #           =       1.08 Main gear in fairings (e.g. C-130)
    #           =       1.03 Main gear retracted in nacelles of turboprop engines
    
    #r_t        =       Correction factor for the tail
    #           =       1.24
    
    #r_w        =       Correction factor for the wing drag
    #           =       1. for a cantilever wing
    #           =       1.1 for a braced wing          
                      
    #r_f        =       Correction factor for the fuselage drag
    #           =       1.0 for a cylindrical fuselage
    #           =       0.65 + 1.5*(d_f/l_f) for a fuselage with streamlined fairings
    
    #r_thr      =       Correction factor for thrust reversers
    #           =       1.0 with thrust reversers
    #           =       0.82 without thrust reversers
    
    #r_n        =       Correction factor for engine configuration
    #           =       1.5 for all engines podded
    #           =       1.65 2 podded and 1 buried in tail
    #           =       1.25 buried in nacelles on fuselage
    #           =       1.0 fully buried engines
    
    #dCD_comp   =       Correction factor for the drag coefficient due to compressibility effects
    #           =       0.0005 for long range cruise conditions
    #           =       0.0020 for high speed cruise conditions
    
    #Assumed values for now
    r_re = 1.0     #no boundary layer affected at take-off speed
    r_uc = 1.0
    r_t = 1.24
    r_w = 1.
    r_f = 1.0
    r_thr = 1.0
    r_n = 1.5
    dCD_comp = 0.   #TO conditions so not yet in compressibility region
    
    
    #Compute parameters in equation
    delta = 1.0                     #p/p0 since p0=p during take-off
    k_w = MTOW / W_TO
    
    #Drag area coefficients     
    CD_Sw = 0.0054*r_w*(1 + 3.*t_c*np.cos(qcsweep)**2)*S                    #For the wing
    CD_Sf = 0.0031*r_f*l_f*(b_f + h_f)                                      #For the fuselage
    CD_Sn = 1.72*r_n*r_thr*((5.+ bypass)/(1.+ bypass))*(T_TO/(psi_TO*p0))  #For the nacelles/engines
    
    #Variables in equation    
    d1 = r_re*r_uc*r_t*(CD_Sw/S) + dCD_comp         #Wing profile drag
    d2 = r_re*r_uc*r_t*((CD_Sf*p0)/MTOW)            #Induced and fuselage drag
    d3 = r_re*r_uc*((CD_Sn*p0)/T_TO)                #Empennage drag

    #To simplify fraction and not have an incredibly long equation = (a+b+c)/d
    a = (0.7*delta*d1*M**2)/(MTOW/(p0*S))
    b = (MTOW/(p0*S))/(0.7*k_w**2*delta*M**2*A*e*np.pi)
    c = 0.7*delta*M**2*d2
    d = (T/T_TO) - 0.7*delta*M**2*d3

    T_W_TO = (a+b+c)/d
    
    return T_W_TO
   
   
"""Take-off W/S limit: for given take-off distance"""
def WS_TO(W_TO, S, rho, CL_maxto,bypass, T_TO,A,S_to,g):
    Vs = np.sqrt((W_TO*2.)/(S*rho*1.3*CL_maxto))    
    V_LOF = 1.2*Vs
    T_mean = 0.75*((5. + bypass)/(4.+ bypass))*T_TO  
    
    gamma_LOF = 0.9*(T_mean/W_TO) - (0.3/np.sqrt(A))
    V3 = V_LOF * np.sqrt(1. + gamma_LOF*np.sqrt(2.))
    mu_dash = 0.02 + 0.01*CL_maxto
    
    h_to = 10.7 #[m] screen height of 35 ft
    f_to = 1. #or 1.15
    
    a = (S_to/f_to) - (h_to/gamma_LOF)
    b = rho*g*CL_maxto * (1. + gamma_LOF*np.sqrt(2.))
    c = (V3/Vs)**2 * ((T_mean/W_TO-mu_dash)**(-1.) + np.sqrt(2.))

    WS = a * (b/c)
    return WS
    
   
    
#------------------------------MAIN PROGRAM------------------------------------
M = Mach(V,0)
rho = ISA_density(0)
p0 = ISA_press(0)
S_to = 2700.

TW_TO = TW_TO_jet(T,T_TO,MTOW,W_TO,M,t_c,qcsweep,S,A,l_f,b_f,h_f,psi_TO,bypass,e,rho,p0)

WS_TO = WS_TO(W_TO, S, rho, CL_maxto,bypass, T_TO,A,S_to,g)

print TW_TO, WS_TO
    
    
    
    
    
    
    
    




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




