# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 16:46:31 2019

@author: nikki
"""

import numpy as np
import matplotlib.pyplot as plt
#----------------------------------INPUTS-------------------------------------
#t_c = thickness over chord ratio of airfoil
#Cs1 = location front spar as percentage of chord
#Cs2 = location aft spar as percentage of chord
#Cr = root chord [m]
#Ct = tip chord [m]
#D_fus = fuselage outer diameter [m]
#bf = wing span [m]
#y_eng = spanwise location of engine [m]
#Mfuel = fuel mass [kg]
#rho_fuel = density of jet fuel [kg/m^3]

t_c = 0.13
Cs1 = 0.15
Cs2 = 0.55
Cr = 5.928
Ct = 1.85
D_fus = 6.7
bf = 58.37
y_eng = 6.
Mfuel = 193990./9.81
rho_fuel = 804. #kg/m^3

#--------------------------------MAIN PROGRAM--------------------------------
def Vintegral(y):
    V = Cr**2*y + Cr*((Ct-Cr)/(bf/2.))*y**2 + (1./3.)*y**3 *((Ct-Cr)/(bf/2.))**2
    return V

V1 = (t_c*(Cs2-Cs1))*(Vintegral((y_eng - 1.)) - Vintegral((D_fus/2.)))
V2 = (t_c*(Cs2-Cs1))*(Vintegral((bf/2. - 1.)) - Vintegral((y_eng + 1.)))


#V1 = 0
#V2 = (t_c*(Cs2-Cs1))*(Vintegral((bf/2. - 1.)) - Vintegral((y_eng+1.)))

Vwing = (V1 + V2)
Vsurge = (Vwing/102.)*2.
Vfuel = (Vwing - Vsurge)*2.

Vreq = (Mfuel/rho_fuel)*1.02

print ('Wing tank volume:', Vfuel, "m^3")
print ("Required volume:", Vreq, "m^3")





#
#g = 9.81
#m =1.3
#
#
#def ISA_density(h):             # enter height in m
#    M = 0.0289644               #kg/mol molar mass of Earth's air
#    R = 8.3144590               #universal gas constant Nm/MolK
#    
#    if h < 11000:
#        rho0 = 1.225            #kg/m^3
#        T = 288.15              #K
#        h0 = 0.                 #m
#        a = -0.0065             #K/m
#        rho = rho0 * (T/(T + a*(h-h0)))**(1. + ((g*M)/(R*a)))
#        
#    if h >= 11000:
#        rho0 = 0.36391          #kg/m^3
#        T = 216.65              #K
#        h0 = 11000.             #m
#        rho = rho0*np.e**((-g*M*(h-h0))/(R*T))
#        
#    return rho
#    
#    
#def ISA_temp(h):
#    if h < 11000:
#        T = 288.15 - 0.0065*h   #in Kelvin
#        return T
#    if h >= 11000:
#        return 216.65           #in Kelvin
#
#def ISA_press(h):
#    p0 = 101325.0
#    T0 = 288.15 #[K]
#    T = ISA_temp(h)
#    a = -0.0065
#    R = 8.3144590               #universal gas constant Nm/MolK
#    
#    if h <= 11000:
#        p = p0*(T/T0)**(-g/(a*R))
#    else:
#        p=p0*np.e**(((-g/(R*T))*(h)))
#        
#    return p
#
#def T_alt(T0,h):
#    T = T0*(ISA_density(h)/ISA_density(0))**m
#    return T
#    
#
#H = np.arange(0,12000,1000)
#T0 = 200e03
#
#T_list = []
#H_list = []
#for h in H:
#    T_list.append(T_alt(T0,h))
#    H_list.append(h)
#    
#
#plt.plot(H_list,T_list)
#plt.show()
#
#print T_alt(T0,9000)/T0


    












