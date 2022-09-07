# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 20:14:52 2019

@author: mathi
"""
import numpy as np
#<<<<<<< HEAD
#from loading_and_moment_diagrams import c
#import loading_and_moment_diagrams as lm
#
#N= 100
#b = 60
#X_root = np.linspace(0, b / 2 - 0.00001, N)
#

nr_spars = 2
first_spar = 0.15
last_spar = 0.55

#<<<<<<< HEAD
#def c(z, Cr, b, taper):
#    Ct = Cr * taper
#    c = Cr - ((Cr - Ct) / (b / 2)) * z
#    return c
#
#def spar_loc(X_root, nr_spars, first_spar, last_spar):
#    Cr = 6.14
#    taper = 0.297 
# 
#=======
def spar_loc(N, b, nr_spars, first_spar, last_spar, c, X_root):
    spar_loc_sec = []
    for i in range(len(X_root)):
        first_spar_location = first_spar*c(X_root[i])
        last_spar_location = last_spar*c(X_root[i])
        delta_spar = (last_spar_location-first_spar_location)/(nr_spars-1)

        spar_loc = np.arange(first_spar_location,last_spar_location + delta_spar -0.001,delta_spar)

        spar_loc_sec.append(spar_loc)
        
    return spar_loc_sec, delta_spar

#<<<<<<< HEAD
#spar_loc_sec, delta_spar = spar_loc(X_root, nr_spars, first_spar, last_spar)
#print(spar_loc_sec)
#=======
#spar_loc_sec, delta_spar = spar_loc(X_root, nr_spars, first_spar, last_spar)
#>>>>>>> master
