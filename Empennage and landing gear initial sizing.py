# -*- coding: utf-8 -*-
"""
Created on Mon May  6 13:32:31 2019

@author: Hidde
"""
# Inputs

MAC      = 4.         # mean aerodynamic chord [m]
l_fuse   = 40.        # length fuselage [m]
x_eng    = -1.        # x-location engines w.r.t. XLEMAC [m]
l_n      = 2.         # length nacelle [m]          
xcg_OEW  = 0.25       # initial cg location OEW w.r.t. MAC [m]

# Mass fractions main components, landing gear neglected 

mass_frac_wing = 0.117
mass_frac_emp  = 0.023
mass_frac_fuse = 0.098
mass_frac_nac  = 0.018
mass_frac_prop = 0.072
mass_frac_fix  = 0.118


# xcg-locations for wing mounted engines aircraft

# Wing group - in [m] w.r.t. X-LEMAC
xcg_wing = 0.4*MAC
xcg_prop = x_eng + 0.4*l_n
xcg_nac  = x_eng + 0.4*l_n

# Fuselage group - in [m] w.r.t. nose
xcg_fuse = 0.4*l_fuse
xcg_fix  = 0.4*l_fuse
xcg_emp  = 0.9*l_fuse

# Averaged mass fractions and cg locations

# Wing group

mass_frac_winggroup = mass_frac_wing + mass_frac_prop + mass_frac_nac
xcg_winggroup       = (xcg_wing*mass_frac_wing + xcg_prop*mass_frac_prop + xcg_nac*mass_frac_nac)/mass_frac_winggroup

# Fuselage group

mass_frac_fusegroup = mass_frac_fuse + mass_frac_fix + mass_frac_emp
xcg_fusegroup       = (xcg_fuse*mass_frac_fuse + xcg_fix*mass_frac_fix + xcg_emp*mass_frac_emp)/mass_frac_fusegroup

# X location of the leading edge of the mean aerodynamic chord

X_LEMAC = xcg_fusegroup + MAC*((xcg_winggroup/MAC)*(mass_frac_winggroup/mass_frac_fusegroup) - xcg_OEW*(1 + (mass_frac_winggroup/mass_frac_fusegroup)))

print(X_LEMAC)