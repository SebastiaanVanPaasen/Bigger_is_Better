# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 13:59:37 2019

@author: mathi
"""

sigma_fatigue_hoop = 350 * 10**6 # look up from graph
sigma_fatigue_long = 300 * 10**6
internal_p = 78.2 * 10**3 
external_p =30.1 * 10**3
R = 5.95/2

def t_min_pressure(sigma_fatigue_hoop, sigma_fatigue_long,R,internal_p, external_p):
    delta_p = internal_p - external_p
    t_min_hoop = delta_p*R/sigma_fatigue_hoop
    t_min_long = delta_p*R/(2*sigma_fatigue_long)
    return t_min_hoop, t_min_long

t_min_hoop,t_min_long = t_min_pressure(sigma_fatigue_hoop, sigma_fatigue_long,R,internal_p, external_p)
print(t_min_hoop)
print(t_min_long)


