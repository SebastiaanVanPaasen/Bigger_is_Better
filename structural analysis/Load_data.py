# -*- coding: utf-8 -*-
"""
Created on Mon May 27 17:20:51 2019

@author: mathi
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 11:10:28 2019

@author: Max
"""
import os

def load(file_name):
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
#    rel_path = "design_results/"+str(file_name)
    rel_path = str(file_name)
    abs_file_path = os.path.join(script_dir, rel_path)


    f = open(abs_file_path, 'r')
    lines = f.readlines()

    f.close
    data = []
    for line in lines:
        x = line.split()
        data.append(x)  
    
    inputs = []
        
    for i in range(len(data)-3):
        inputs.append(float(data[i][-1]))
        
    Mach = inputs[0]
    cruise_alt = inputs[1]
    W = inputs[2]
    OEW = inputs[3]
    Wing_W_frac = inputs[4] 
    Fuselage_weight = inputs[5]
    Empennage_weight =inputs[6] 
    Nacelle_weight = inputs[7]
    Engine_weight = inputs[8]
    Landing_weight = inputs[9]
    Fixed_equip_weight =inputs[10] 
    Prop_weight =inputs[11]
    Payload_weight = inputs[12]
    Fuel_W_tot = inputs[13]
    total_thrust = inputs[14]
    Cl_cruise = inputs[15]
    Cd_cruise= inputs[16]
    rho = inputs[17] 
    CD0 = inputs[18]
    Lift_over_Drag = inputs[19]
    V_cruise = inputs[20]
    n_ult = inputs[21]
    V = inputs[22]
    b = inputs[23]
    S = inputs[24]
    Fuselage_length =inputs[25] 
    Cr =inputs[26]
    Sweep0 = inputs[27]
    taper = inputs[28]
    AR = inputs[29]
    e = inputs[30]
    Wing_W = inputs[31]
    n_engines = inputs[32]
    engine_weigth = inputs[33]
    
    
    return Mach,cruise_alt,W,OEW,Wing_W_frac,Fuselage_weight,Empennage_weight,Nacelle_weight ,Engine_weight,Landing_weight ,Fixed_equip_weight ,Prop_weight,Payload_weight,Fuel_W_tot,total_thrust,Cl_cruise,Cd_cruise,rho,CD0,Lift_over_Drag,V_cruise,n_ult,V,b,S,Fuselage_length,Cr,Sweep0,taper,AR,e,Wing_W,n_engines,engine_weigth 
    
inputs = load('aerodynamic_concept')
print (inputs)

#file_name = 'aerodynamic_concept'
#
#rel_path = "design_results/"+str(file_name)
#
#
#script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
#abs_file_path = os.path.join(script_dir, rel_path)
#print (rel_path,abs_file_path )