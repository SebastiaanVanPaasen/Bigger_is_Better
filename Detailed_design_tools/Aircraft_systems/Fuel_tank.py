# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 16:46:31 2019

@author: nikki
"""

import numpy as np

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
Cs1 = 0.2
Cs2 = 0.6
Cr = 6.0396326
Ct = 1.88879566
D_fus = 6.2877
bf = 59.463212
y_eng = bf/4.
Mfuel = 199151.8874/9.81
rho_fuel = 804. #kg/m^3

#--------------------------------MAIN PROGRAM--------------------------------
def Vintegral(y):
    V = Cr**2*y + Cr*((Ct-Cr)/(bf/2.))*y**2 + (1./3.)*y**3 *((Ct-Cr)/(bf/2.))**2
    return V

V1 = (t_c*(Cs2-Cs1))*(Vintegral((y_eng - 1.)) - Vintegral((D_fus/2.)))
V2 = (t_c*(Cs2-Cs1))*(Vintegral((bf/2. - 1.)) - Vintegral((y_eng + 1.)))

Vwing = (V1 + V2)
Vsurge = (Vwing/102.)*2.
Vfuel = (Vwing - Vsurge)*2.

Vreq = Mfuel/rho_fuel

print ('Wing tank volume:', Vfuel, "m^3")
print ("Required volume:", Vreq, "m^3")

















