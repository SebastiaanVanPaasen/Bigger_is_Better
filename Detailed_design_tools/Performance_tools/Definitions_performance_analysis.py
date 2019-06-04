# -*- coding: utf-8 -*-
"""
Created on Mon May 27 12:41:44 2019

@author: nikki
"""
#-----------------------------IMPORT MODULES----------------------------------
import numpy as np



#---------------------------DEFINITIONS---------------------------------------


"""Take-off T/W jet"""
#Inputs needed for def jet aircraft
#MTOW       =       Maximum take-off weight [N]
#W_TO       =       Weight at take-off [N]
#T_TO       =       Total static thrust of all engines at take-off [N]

#Vcr        =       Cruise velocity TAS [m/s]
#M          =       Cruise Mach number [-]
#A          =       Aspect ratio [-]
#e          =       Oswald efficiency factor [-]

#t_c        =       Thickness over chord ration airfoil
#qcsweep    =       Quater chord sweep in [rad]
#S          =       Surface area wing [m^2]
#b_f        =       Width of fuselage [m]
#h_f        =       Height of fuselage [m]
#l_f        =       Length of fuselage [m]

#psi_TO     =       Specific Thrust N/airflow [N/kg/s?]
#phi_TO     =       Power/frontal area of the engines
#bypass     =       Bypass ratio of the engine

#p          =       Static atmospheric pressure [Pa]
#p0         =       Static atmospheric pressure at sea level ISA [Pa]
#rho        =       Atmospheric density [kg/m^3]
#mu_cr      =       Dynamic viscosity during cruise [Pa*s] or [Ns/m^2] or [kg/m*s]

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

def TW_TO_jet():                                #T/W for TO of a jet aircraft
    #Compute parameters in equation
    delta = p/p0
    k_w = MTOW / W_TO
    
    Re_f = V_cr*l_f*rho/mu_cr       #Fuselage reynolds number
    r_re = 47.*Re_f**-0.2           #Correction factor for the Reynolds number
    
    #Drag area coefficients     
    CD_Sw = 0.0054*r_w*(1 + 3.*t_c*np.cos(qcsweep)**2)*S                    #For the wing
    CD_Sf = 0.0031*r_f*l_f*(b_f + h_f)                                      #For the fuselage
    CD_Sn = 1.72*r_n*r_thr*((5.+ bypass)/(1.+ bypass))*(T_TO/(psi_TO*p0))   #For the nacelles/engines
    
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
    


"""Take-off T/W properller aircraft""" 
#(Additional) Inputs needed for def of propeller aircraft
#P_TO           =       Maximum take-off thrust [W] 
#P/P_TO         =       Power over max. TO power from manufacturer's data
#eta_prop       =       Propeller efficiency
#               =       0.85 turboprop with jet thrust
#               =       0.8 wing mounted piston engines  

#r_n_prop   =       Correction factor for the engine if turboprop
#           =       1.0 for ring type inlet
#           =       1.6 for scoop type inlet
   
def PW_TO_prop():                               #T/W for TO for a propeller aircraft
    #Compute parameters in equation
    k_w = MTOW / W_TO
    
    Re_f = V_cr*l_f*rho/mu_cr       #Fuselage reynolds number
    r_re = 47.*Re_f**-0.2           #Correction factor for the Reynolds number
    
    #Drag area coefficients     
    CD_Sw = 0.0054*r_w*(1 + 3.*t_c*np.cos(qcsweep)**2)*S            #For the wing
    CD_Sf = 0.0031*r_f*l_f*(b_f + h_f)                              #For the fuselage
    CD_Sn = 0.1*r_n_prop*(P_TO/phi_TO)                              #For the engines/nacelles
 
    #Variables in equation    
    d1 = r_re*r_uc*r_t*(CD_Sw/S) + dCD_comp         #Wing profile drag
    d2 = r_re*r_uc*r_t*((CD_Sf*p0)/MTOW)            #Induced and fuselage drag
    d3 = r_re*r_uc*((CD_Sn)/P_TO)                   #Empennage drag       
        
    #To simplify fraction and not have an incredibly long equation = (a + b)/c
    a = 0.5*rho*Vcr**2*((d1/(MTOW/S))+(d2/p0))
    b = (2.*MTOW/S)/(k_w**2*rho*Vcr**2*np.pi*A*e)
    c = eta_prop *(P/P_TO) - 0.5*rho*Vcr**3*d3
    P_W_TO = ((a + b)/c)*Vcr
    
    return P_W_TO





"""Climb T/W: for specified Rate of Climb RC for DV/Dh = 0 and n = 1"""
#Inputs needed for jet climb
#W          =       Weight during climb, dependent on flight phase
#CD0        =       Zero-lift drag
#S          =       Wing surface area [m^2]
#A          =       Aspect ratio [-]
#e          =       Oswald efficiency factor [-]

#gamma      =       Climb gradient in [rad]? or specific heats of 1.4

#a          =       Speed of sound [m/s]
#p          =       Atmospheric static pressure [Pa]

#C          =       Climb rate 
#vg_dvdh    =       Constant for the climb rate
#           =       0.5668*M**2 for const. EAS in troposphere
#           =       -0.1332*M**2 for cons. M in troposphere
#           =       0.7*M**2    for const. EAS in stratosphere
#           =       0           for const. M in stratosphere

def RC():
    C = ((T-D)*V)/(W*(1. + vg_dvdh))
    return C

def TW_steady_climb_jet(): #Steady flight climb at LOW ALTITUDE

    #To simplify fraction and not have an incredibly long equation
    x1 = ((C/a)*np.sqrt(CD0))/(np.sqrt(W/(p*S)))
    x2 = np.sqrt(W/(p*S))/((C/a)*np.sqrt(CD0))
    T_W = (3./2.)*gamma**(1./3.)*x1**(2./3.) + (2./gamma**(1./3.))*(CD0/(np.pi*A*e))*x2**(2./3.)
    
    return T_W



#Additional inputs for service ceiling thrust
#n          =       Load factor L/W
#theta      =       T/T0 relative atmospheric pressure    

def TW_ceiling_climb_jet():  #HIGH ALTITUDE indicates service ceiling thrust
    T_W = 2.*n*np.sqrt(CD0/(np.pi*A*e)) + \\
    (0.00123/np.sqrt(theta))*(((CD0*np.pi*A*e)**0.25)/(np.sqrt((n*W)/(p*S))))
    return T_W
    
#To fly at constant EAS instead of constant TAS add this term
def dTW_constEAS_jet():
    dT_W = 0.567*C/a*np.sqrt(((C/a)*(W/(p*S)))/(gamma*CD0))
    return dT_W
    
    
    
    
    
"""Climb performance of propeller aircraft"""
#Inputs needed for propeller climb
#W          =       Weight during climb, dependent on flight phase
#CD0        =       Zero-lift drag
#S          =       Wing surface area [m^2]
#A          =       Aspect ratio [-]
#e          =       Oswald efficiency factor [-]

#gamma      =       Climb gradient in [rad]? or specific heats of 1.4

#a          =       Speed of sound [m/s]
#p          =       Atmospheric static pressure [Pa]

#C          =       Climb rate 
#r_n_prop   =       Correction factor for the engine if turboprop
#           =       1.0 for ring type inlet
#           =       1.6 for scoop type inlet

#C/a0       =       0.00147 (a0 is speed of sounds at sea level)

def PW_climb_prop():
    PW = (1./eta_prop)*((C/a) + 2.217*((CD0**0.25)/(np.pi*A*e)**0.75) * np.sqrt(W/(p*S)))
    return PW
    



"""Fuel weight: cruise and reserve fuel"""
#Inputs 
#R      =       Design range in [m]
#Ct     =       Thrust specific fuel consumption [kg/N/s]
#g      =       9.81 [m/s^2]
#T      =       Static air temperature at atitude [K]
#T0     =       Static air temperature at sea level [K]
#a0     =       Speed of sound at sea level
#M      =       Mach number 
#CD     =       Drag coefficient
#CL     =       Lift coefficient 
#A      =       Aspect ratio


def Wf_Wto():       #To compute the fuel weight for the trip
    
    #Assumptions:
    #Cruise fuel is the dominant fuel weight
    #Wetted area of the fuselage is the main contributer for fuselage drag   
    
    theta = T/T0
    Wfcr_Wto = 1. - np.exp**(((-R*Ct*g/np.sqrt(theta))/(a0*M))*(CD/CL))
    Wfres_Wto = 0.18 * ((Ct*g/np.sqrt(theta))/np.sqrt(A))
    
    Wf_Wto = Wfcr_Wto + Wfres_Wto
    return Wf_Wto
    
    
"""Take-off distance: all engines functioning"""
#Inputs
#W_to       =       Take-off weight (MTOW)
#S          =       Wing surface area [m^2]
#rho0       =       Sea level density [kg/m^3]
#CL_maxto   =       CL max during take-off
#bypass     =       Bypass ratio of engine
#T_to       =       Take-off thrust [N]
#g          =       9.81 [m/s^2]
#A          =       Aspect ratio

def TO_distance():       #Required distance to pass screen height (30 ft) at a speed of 1.3*VV_stall
    Vs = np.sqrt((W_to*2.)/(S*rho0*CL_maxto))    
    V_LOF = 1.2*Vs
    T_mean = 0.75*((5. + bypass)/(4.+ bypass))*T_to    
    mu_dash = 0.02 + 0.01*CL_maxto
    
    S_run = (V_LOF**2/(2.*g))/((T_mean/W_to)-mu_dash)
    
    gamma_LOF = 0.9*(T_mean/W_to) - (0.3/np.sqrt(A))
    h_to = 10.7 #[m] screen height of 35 ft
    S_air = (V_LOF**2)/(g*np.sqrt(2)) + (h_to/gamma_LOF)
    
    S_to = S_run + S_air
    
    return S_to
    
"""Take-off W/S limit: for given take-off distance"""
#Inputs
#W_to       =       Take-off weight (MTOW)
#S          =       Wing surface area [m^2]
#rho0       =       Sea level density [kg/m^3]
#CL_maxto   =       CL max during take-off
#g          =       9.81 [m/s^2]
#A          =       Aspect ratio


def WS_TO():
    Vs = np.sqrt((W_to*2.)/(S*rho0*CL_maxto))    
    V_LOF = 1.2*Vs
    gamma_LOF = 0.9*(T_mean/W_to) - (0.3/np.sqrt(A))
    V3 = V_LOF * np.sqrt(1. + gamma_LOF*np.sqrt(2.))
    mu_dash = 0.02 + 0.01*CL_maxto
    
    h_to = 10.7 #[m] screen height of 35 ft
    f_to = 1. #or 1.15
    
    a = (S_to/f_to) - (h_to/gamma_LOF)
    b = rho0*g*CL_maxto * (1. + gamma_LOF*np.sqrt(2.))
    c = (V3/Vs)**2 * ((T_mean/W_to-mu_dash)**(-1.) + np.sqrt(2.))

    WS = a * (b/c)
    return WS
 
    
"""TO with engine failure"""
#Inputs
#W_to       =       Take-off weight (MTOW)
#T_to       =       Take-off thrust [N]
#S          =       Wing surface area [m^2]
#rho0       =       Sea level density [kg/m^3]
#CL_maxto   =       CL max during take-off
#g          =       9.81 [m/s^2]
#A          =       Aspect ratio
#CD_to      =       Average CD during take-off
#bypass     =       Bypass ratio of engine
#a_stop     =       Average decelearation usually 0.37g
#a_mean     =       Mean acceleration assume 0.25g??
#h_to       =       10.7 [m]

#gamma2_min =       Correction factor
#           =       0.024; 0.027; 0.030
#           =       2; 3; 4 number of engines

#dt         =       4.5 [s] reaction time
#mu         =       Ground friction coefficient
#           =       0.02 for concrete

def TO_eng_fail():
    Vs = np.sqrt((W_to*2.)/(S*rho0*CL_maxto))    
    V2 = 1.2*Vs
    CL_to = mu*np.pi*A*e
    
    dgamma2 = ((T_to/W_to)-(CL_to/CD_to)**(-1.)) - gamma2_min
    gamma_mean = 0.06 + dgamma2
    
    Vx = V2*( ((1. + 2.*g*h_to/V2**2)/(1. + gamma_mean / (a_mean/g)))**0.5 - ((gamma_mean*g*(dt - 1.))/V2) )
    
    S01 = Vx**2 / (2.*a_mean)                               #Distance covered before engine failure at Vx
    S12 = (1./gamma_mean)*(((V2**2-Vx**2)/(2.*g)) + h_to)   #Distance from engine failure up to save screen height at V2
    Sstop = (Vx**2/(2*a_stop)) + Vx*dt                      #If TO aborted (stop distance needed)
    
    S_continue = S01 + S12
    S_abord = S01 + Sstop
    
    return S_continue, S_abord

def BFL():      #Balanced field length: is when taking off is better than stopping in case of engine failure
    CL_to = mu*np.pi*A*e    
    dgamma2 = ((T_to/W_to)-(CL_to/CD_to)**(-1.)) - gamma2_min    
    CL2 = 0.694*CL_maxto
    T_mean = 0.75*((5. + bypass)/(4.+ bypass))*T_to    
    mu_dash = 0.02 + 0.01*CL_maxto
    
    a = 0.863/(1. + 2.3*dgamma2)    
    b = ((W_to/S)/(rho0*g*CL2)) + h_to    
    c = (1./(T_mean/W_to - mu_dash)) + 2.7
    
    BFL = a*b*c
    return BFL
    
    
"""Landing distance"""
#Inputs
#gamma_mean =       During landing is approx. 0.1
#h_land     =       15.3 [m] screen height for landing
#dn         =       0.1 (change in load factor during landing)
#a_stop/g   =       0.35-0.45 without thrust reversers
#           =       0.40 - 0.50 thrustreversers and spoilers
#           =       0.50 - 0.60 + nose wheel braking

def S_land():
    Va = 1.3*Vs
    Vtd = np.sqrt(Va**2*(1. - (gamma_mean**2/dn)))    

    S_air = (1./gamma_mean)*(((Va**2 - Vtd**2)/(2.*g)) + h_land)
    S_run = Vtd**2 / (2.*a_stop)
    
    SL = S_air + S_run
    return SL

"""ISA definitions"""

 def ISA_density(h):             # enter height in m
    M = 0.0289644               #kg/mol molar mass of Earth's air
    R = 8.3144590               #universal gas constant Nm/MolK
    
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
    p0 = 101325.0
    T0 = 293.15 #[K]
    T = ISA_temp(h)
    
    if h <= 11000:
        p = p0*(T/T0)**(-g/(a*R))
    else:
        p=p0*math.e**(((-g/(R*T))*(h)))
        
    return p
 
          
def Mach(V,h):                  #enter V in km/h
    gamma = 1.4                 #enter h in m
    R = 287 #J/kg/K
    a = np.sqrt(gamma*R*ISA_temp(h))
    M = (V/3.6)/a
    return M 







