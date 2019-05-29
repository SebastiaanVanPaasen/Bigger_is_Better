# -*- coding: utf-8 -*-
"""
Created on Tue May 28 09:19:12 2019

@author: mathi
"""

import loading_and_moment_diagrams as lm
from loading_and_moment_diagrams import load_diagrams
import numpy as np

def required_Izz(Cr):
    M_max = load_diagrams(100)[1]
    y_max = 0.07 * Cr
    # print(y_max)
    sigma_ult = 552 * 10 ** 6
    I_zz = M_max * y_max / sigma_ult
    return 'I_zz=', I_zz

 #8.54#7.11#6.63#6.06
print("req Izz",required_Izz(lm.Cr))
#print("req Izz",required_Izz(5.8451))

strut_length = 20.95

def required_strut_area(strut_length):
    angle = (20.78/180)*np.pi
    sigma_carbon = 1500*10**6 
    density_carbon = 1600
    cost_carbon = 100
    P_strut = lm.strutforce/np.sin(angle)
    A_req = P_strut/sigma_carbon
    strut_volume = A_req*strut_length
    strut_mass = strut_volume*density_carbon
    strut_cost = strut_mass*cost_carbon
    print(lm.strutforce)
    print(P_strut)
    return 'A_req=',A_req, 'strut_mass=',strut_mass, 'strut_cost=',strut_cost

print(required_strut_area(strut_length))