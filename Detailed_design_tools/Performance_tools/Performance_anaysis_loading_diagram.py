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
T          =       107*1000        #Total static thrust of all engines at take-off [N]
T_TO       =       91.63*1000*2    #Take-off thrust
T_cr       =       0.2789*0.8*MTOW #Thrust during cruise, assume a weight of 80%MTOW

CL_maxto   =        2.5
CL_maxcr   =        0.9       #Max CL during cruise
CD0        =        0.02
CD         =        0.06      #Approximate CD for just after take-off, take-off  climb
 
V_to       =       73.25      #Speed during take-off [m/s]       
V_cr       =       257.22     #Cruise speed [m/s]             
A          =       9.44       #Aspect ratio [-]
e          =       0.85       #Oswald efficiency factor [-]
n          =        2.        #load factor L/W during cruise 
hcr        =       12000      #Cruise altitude [m]

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
        T = 293.15 - 0.0065*h   #in Kelvin
        return T
    if h >= 11000:
        return 216.65           #in Kelvin

def ISA_press(h):
    R = 8.3144590               #universal gas constant Nm/MolK
    M = 0.0289644               #kg/mol molar mass of Earth's air
    g = 9.80665                 #gravitationa constant 
    
    p0 = 101325.0
    T0 = 288.15 #[K]
    a = -0.0065 
    T1 = ISA_temp(11000)
    
    if h == 0:
        p = p0
    elif h <= 11000:
        p = p0*(T0/(T0 + a*h))**((g*M)/(R*a))
        #p = p0*(T/T0)**(-g/(a*R))
    else:
        p = 22632.10 *np.exp((-g*M*(h - 11000.))/(R*T1))
        
    return p
         

def Mach(V,h):                  #enter V in m/s
    gamma = 1.4                 #enter h in m
    R = 287 #J/kg/K
    a = np.sqrt(gamma*R*ISA_temp(h))
    M = (V/a)
    return M 
    
def Sos(h):
    R = 287.                    #universal gas constant J/K/kg
    gamma = 1.4                 #Specific heat constant
    T = ISA_temp(h)
    a = np.sqrt(gamma*R*T)
    return a
    


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
    
  
  
  
"""Climb T/W: for specified Rate of Climb RC for DV/Dh = 0 and n = 1"""
#C          =       Climb rate 
#vg_dvdh    =       Constant for the climb rate
#           =       0.5668*M**2 for const. EAS in troposphere
#           =       -0.1332*M**2 for cons. M in troposphere
#           =       0.7*M**2    for const. EAS in stratosphere
#           =       0           for const. M in stratosphere

def RC(T,D,V,W,vg_dvdh):
    C = ((T-D)*V)/(W*(1. + vg_dvdh))
    return C

def TW_steady_climb_jet(C,a,CD0,W,p,S,A,e): #Steady flight climb at LOW ALTITUDE
    gamma = 1.4
    #To simplify fraction and not have an incredibly long equation
    x1 = ((C/a)*np.sqrt(CD0))/(np.sqrt(W/(p*S)))
    x2 = np.sqrt(W/(p*S))/((C/a)*np.sqrt(CD0))
    T_W = (3./2.)*gamma**(1./3.)*x1**(2./3.) + (2./gamma**(1./3.))*(CD0/(np.pi*A*e))*x2**(2./3.)
    
    return T_W

def TW_ceiling_climb_jet(n,CD0,A,e,theta,W,p,S):  #HIGH ALTITUDE indicates service ceiling thrust
    T_W = 2.*n*np.sqrt(CD0/(np.pi*A*e)) + (0.00123/np.sqrt(theta))*(((CD0*np.pi*A*e)**0.25)/(np.sqrt((n*W)/(p*S))))
    return T_W
    
#To fly at constant EAS instead of constant TAS add this term
def dTW_constEAS_jet(C,a,W,p,S,CD0):
    gamma = 1.4
    dT_W = 0.567*C/a*np.sqrt(((C/a)*(W/(p*S)))/(gamma*CD0))
    return dT_W
      
    
#------------------------------MAIN PROGRAM------------------------------------
"""Climb performance"""
#For steady low climb
Vs = np.sqrt((W_TO*2.)/(S*ISA_density(0)*CL_maxto))    
V2 = 1.2*Vs                     #"Fly away" speed
D_TO = 0.5*ISA_density(0)*V2**2.*S*CD

M = Mach(V2,0)
vg_dvdh = 0.5668*M**2            #Constant EAS in tropospere

a_steady = Sos(0) 
p_steady = ISA_press(0)

C_steady = RC(T_TO,D_TO,V2,W_TO,vg_dvdh)      #Rate of climb for steady low altitude

#steady climb low altitude, at constant EAS
TW_steady_climb = TW_steady_climb_jet(C_steady,a_steady,CD0,W_TO,p_steady,S,A,e) + dTW_constEAS_jet(C_steady,a_steady,W_TO,p_steady,S,CD0) 

print ("Rate of climb steady low altitude at const. EAS: ", C_steady, "m/s")
print ("T/W for steady climb :", TW_steady_climb)
print ()

#For service ceiling climb
Dcr = 0.5*ISA_density(hcr)*V_cr**2.*S*CD
Mcr = Mach(V_cr,hcr)
vg_dvdh = 0.7*Mcr**2    
theta = ISA_temp(12000)/ISA_temp(0)                #Constant EAS in stratosphere

a_high = Sos(hcr)
p_high = ISA_press(hcr)

C_high = RC(T_cr,D_TO,V_cr,0.8*W_TO,vg_dvdh)      #Rate of climb for ceiling

#High altitude climg: service ceiling thrust at constant EAS 
TW_ceiling = TW_ceiling_climb_jet(n,CD0,A,e,theta,W_TO,p_high,S) + dTW_constEAS_jet(C_high,a_high,W_TO,p_high,S,CD0) 

print ("Rate of climb service ceiling at const. EAS: ", C_high, "m/s")
print( "T/W for ceiling climb :", TW_ceiling)
print()


"""Take-off wing loading diagram"""
rho = ISA_density(0)
p0 = ISA_press(0)

Vs = np.sqrt((W_TO*2.)/(S*rho*CL_maxto))    
V2 = 1.2*Vs
M = Mach(V2,0)

S_to = 2700.        #Specified max.take-off length to operate from major hubs

TW_TO = TW_TO_jet(T,T_TO,MTOW,W_TO,M,t_c,qcsweep,S,A,l_f,b_f,h_f,psi_TO,bypass,e,rho,p0)
WS_TO = WS_TO(W_TO, S, rho, CL_maxto,bypass, T_TO,A,S_to,g)


print ('Take-off T/W and W/S: ', TW_TO, "and", WS_TO)
    

print()
print(" * at const. EAS means that you are accelerating during climb as the density becomes less" )
  
    
    




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




