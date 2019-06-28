# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 12:18:42 2019

@author: nikki
"""
import numpy as np

"""Design inputs"""
#Wto = 1520276.626
#W_LTO = 295.12*9.81
#Wcl = Wto - W_LTO
#T0 = 432.34*100/2
#
#Vcr = 218.7110308
#Vcl = 0.599*Vcr
#
#
#S = 221.21 
#A = 15
#e = 0.6
#CD0 = 0.02
#
#mf = 0.597*0.9

#
#Wto = 1497151.235
#W_LTO = 288.222*9.81
#Wcl = Wto - W_LTO
#T0 = 428.4*100
#
#Vcr = 218.7110308
#Vcl = 0.57*Vcr
#
#
#S = 208.95
#A = 13
#e = 0.74
#CD0 = 0.0233
#
#mf = 0.5247*2
#
##
#h_start = 1000*0.3048
#h_cr = 9000



"""B737-8 MAX inputs"""
Wto = 82191.*9.81
W_LTO = 365.3*9.81
Wcl = Wto - W_LTO
W_des = 64500*9.81

T0 = 130.*1000*2.


Vcr = 233.0555
Vcl = 0.57*Vcr

mf = 0.456*2.3
S = 127
A = 10.16
e = 0.7
CD0 = 0.020

#
h_start = 1000*0.3048
h_cr = 9000


"""Verificitaion A320-200 filippone book"""
#Wto = 78000*9.81
#W_LTO = 350.
#Wcl = Wto - W_LTO
#
#T0 = 95*1000*2
#
#mf = 0.413*2
#Vcl = 137.
#S = 124.
#A = 10.3
#e = 0.8
#CD0 = 0.02
#
#h_start = 3098
#h_cr = 10058



#----------------------------DEFINITIONS--------------------------------------      
""" ISA definitions""" 
def ISA_density(h):      # enter height in m
    M = 0.0289644       #kg/mol molar mass of Earth's air
    R = 8.3144590       #universal gas constant Nm/MolK
    g = 9.81
    if h < 11000:
        rho0 = 1.225   #kg/m^3
        T = 288.15     #K
        h0 = 0.         #m
        a = -0.0065    #K/m
        rho = rho0 * (T/(T + a*(h-h0)))**(1. + ((g*M)/(R*a)))
        
    if h >= 11000:
        rho0 = 0.36391 #kg/m^3
        T = 216.65     #K
        h0 = 11000.     #m
        rho = rho0*np.e**((-g*M*(h-h0))/(R*T))
        
    return rho
    
def ISA_temp(h):
    if h < 11000:
        T = 288.15 - 0.0065*h   #in Kelvin
        return T
    if h >= 11000:
        return 216.65    #in Kelvin
    
def Mach(V,h):                  #enter V in m/s
    gamma = 1.4                 #enter h in m
    R = 287 #J/kg/K
    a = np.sqrt(gamma*R*ISA_temp(h))
    M = V/a
    return M 

def Vel(M,h):
    gamma = 1.4                 #enter h in m
    R = 287 #J/kg/K
    a = np.sqrt(gamma*R*ISA_temp(h))
    V = a*M
    return V
    

#Varying thrust with velocity AT 90% CLIMB THRUST SETTING
def T_alt(T0,h):
    m = 1.3
    T = ((T0)*(ISA_density(h)/ISA_density(0))**m)*0.9
    return T


"""Climb performance"""
#Normal rate of climb at a given altitude and airspeed (steady, no acceleration)
def RC(W,T0,V,S,A,e,CD0,h):
    k = 1. / (np.pi*A*e)
    T = T_alt(T0,h)
    RC = V*((T/W) - 0.5*ISA_density(h)*V**2 * (S/W)*CD0 - (W/S)*((2*k)/(ISA_density(h)*V**2)))
    return RC
    



def glide_range(L_D,dh):  
    R = (L_D)*dh
    return R



#print (glide_range(17.43,9000-1000*0.3048)/1000)

#---------------------------------MAIN PROGRAM--------------------------------------

"""Cimb"""
H = np.arange(h_start,h_cr,205)

RC_list = []
time = []
fuel = []


for i in range(len(H)-1):
    
    RCi = abs(RC(Wcl,T0,Vcl,S,A,e,CD0,H[i]))
    mfi = mf - 0.001*i
    ti = abs((H[i+1] - H[i])/RCi)

    #if i%2. == 0:
    time.append(ti)
    
    dW = mfi*ti

    fuel.append(dW)
    
    Wcl = Wcl - dW
    
    
print (sum(time))
print (sum(fuel))










