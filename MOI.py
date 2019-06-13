# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:12:22 2019

@author: floyd
"""
from matplotlib import pyplot as plt
import numpy as np
#with open('n2414.dat', 'w') as textfile:
def MOI(chordlength, airfoil_t):
    f = open("n2414.dat", "r")
    lines = f.readlines()    
    f.close
    x,y = [],[]
    for line in lines:
        a = line.split()
        x.append(float(a[0])*chordlength)
        y.append(float(a[1])*chordlength)
    x = np.asarray(x)
    y = np.asarray(y)
    x_centroid = 0.42*chordlength
    y_centroid = 1.54E-2*chordlength
    diff_x  = abs(np.diff(x))
    diff_y = abs(np.diff(y))
    x = x[:-1]
    y = y[:-1]
    diff = np.sqrt(diff_x**2 + diff_y**2)
    
    area = diff*airfoil_t
    Ixx = np.sum((y-y_centroid)**2*area)
    Iyy = np.sum((x-x_centroid)**2*area)
    Ixy = np.sum(abs(x-x_centroid)*abs(y-y_centroid)*area)
    J = Ixx + Iyy
#    print("Ixx: ",Ixx,"Iyy: " ,Iyy, "Ixy: ", Ixy,"J: ",J)
    return Ixx, Iyy, Ixy, J 
#print(MOI(1,0.001))
    
def ell_lift(L,b):
    step = 0.001
    y = np.arange(0,b+step,step) #b is semispan
    lift_dist = np.sqrt(1-(y/b)**2)
    total_lift = np.trapz(lift_dist, dx=step)
    scale_factor = total_lift/L
    lift_dist = lift_dist/scale_factor
#    print(np.trapz(lift_dist, dx=step))
#    plt.plot(y,lift_dist)
    return lift_dist,y

def MOI2(chordlength, airfoil_t):
    x_centroid = 0.42*chordlength
    y_centroid = 1.54E-2*chordlength
    total_airfoil_thickness = 0.14*chordlength
    area = 0.1*chordlength*airfoil_t
    y1 = total_airfoil_thickness-y_centroid    
    y2 = total_airfoil_thickness-y1
    x1 = 0.5*chordlength - x_centroid
    Ixx = 0.1*chordlength*airfoil_t**3/12 + (y1**2*area) +(y2**2*area)
    Iyy = 0.1*chordlength**3*airfoil_t/12 + 2*(x1**2*area)
    Ixy = (x1*y1*area + x1*y2*area)
    J = 0.1*chordlength*airfoil_t**3*(16/3 - 3.36*airfoil_t/(0.1*chordlength)*(1-1/12*airfoil_t**4/(0.1*chordlength)**4))
    print("Ixx: ",Ixx,"Iyy: " ,Iyy, "Ixy: ", Ixy,"J: ",J)
    return Ixx, Iyy, Ixy, J
MOI2(3,0.001)
