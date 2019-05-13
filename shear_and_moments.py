# -*- coding: utf-8 -*-
"""
Created on Fri May 10 16:50:26 2019

@author: mathi
"""

n_elements = 1000.

Vz = np.zeros(len(dx))
Vy = np.zeros(len(dx))
Mz = np.zeros(len(dx))
My = np.zeros(len(dx))
T = np.zeros(len(dx)

for i in range(n_elements):
    C_L = (C_L[i+1] - C_L[i])/2.
    L = 0.5*C_L*rho*dx*c
    D = 0.5*C_D*rho*dx*c
    Vz[i + 1] = Vz[i] + D
    Vy[i + 1] = Vy[i] + L