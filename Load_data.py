#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 11:10:28 2019

@author: Max
"""
import os

def load(file_name):
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    rel_path = "design_results/"+str(file_name)
    abs_file_path = os.path.join(script_dir, rel_path)


    f = open(abs_file_path, 'r')
    lines = f.readlines()

    f.close
    data = []
    for line in lines:
        x = line.split()
        data.append(x)  
    
    inputs = []
        
    for i in range(len(data)-1):
        inputs.append(float(data[i][-1]))
        
    return inputs
    
inputs = load('Aerodynamic_concept')
print (inputs)

#file_name = 'Aerodynamic_concept'
#
#rel_path = "design_results/"+str(file_name)
#
#
#script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
#abs_file_path = os.path.join(script_dir, rel_path)
#print (rel_path,abs_file_path )