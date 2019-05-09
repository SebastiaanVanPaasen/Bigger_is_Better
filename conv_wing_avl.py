# -*- coding: utf-8 -*-
"""
Created on Thu May  9 13:38:52 2019

@author: floyd
"""
import numpy as np

M_cruise = 0.75      #inputs from different part 
surface_area = 200  #inputs from different part 
aspect_ratio = 9    #inputs from different part 
CL_cruise = 0.5
CD0 = 0.010
root_chord = 4 #from planform file
tip_chord = 2#planform
Angle = 0.0 #Incidence angle
leading_edge_sweep = 25 #planform file
dihedral = 0.5

span = np.sqrt(surface_area*aspect_ratio)
chords= [root_chord, tip_chord]
Ainc = [0.0, 0.0]
Nspanwise = [0, 0]
Sspace = [0, 0]
x_loc_LE = [0, span/2*np.tan(np.radians(leading_edge_sweep))]
y_loc_LE = [0, span/2]
z_loc_LE = [0,span/2*np.tan(dihedral)] #change for dihedral angle low/high wing


with open("conv_wing.avl", "w") as text_file:
    print("Boxed Wing" +"\n"
        "#Mach" +"\n" + 
        str(M_cruise) +"\n"
        "#IYsym IZsym Zsym" +"\n"
        "0  0  0" +"\n" 
        "#Sref  Cref  Bref" +"\n" +
        str(surface_area), str(root_chord), str(span), "\n"          
        "#Xref Yref Zref" +"\n"
        "0.3  0.0  0.0806" + "\n"
        "CDdp" + "\n" +
        str(CD0), "\n"
        "\n" + "SURFACE" +"\n" +
        "Wing", "\n"
        "8  1.0  12  -2.0"+"\n"
        "YDUPLICATE"+"\n" +
        str(0.0), "\n" +
        "ANGLE"+"\n"+
        str(Angle), file=text_file)
    for i in range(2):
        print("SECTION", file=text_file)            
        print(x_loc_LE[i],y_loc_LE[i],z_loc_LE[i],chords[i],Ainc[i],Nspanwise[i], Sspace[i], file=text_file)        
            