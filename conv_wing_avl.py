# -*- coding: utf-8 -*-
"""
Created on Thu May  9 13:38:52 2019

@author: floyd
"""
import numpy as np
from inputs import *
from planform import wing_parameters
from drag_estimation import *
import subprocess
import os


#MANUALLY PUT IN PARAMETERS from latest iteration
surface_area = 200  #Once this is included in inputs file then this should be deleted 
d_fuselage, l_nosecone, l_cabin, l_tailcone = 5, 10, 50, 15 #same as surface area, should come from iteration.
#GET PLANFORM PARAMETERS
wing_geo  = wing_parameters(M_cruise, CL_cruise, surface_area, A) #planform parameters: 
#0 quarter_chord_sweep, 1 leading_edge_sweep, 2 taper_ratio, 3 span, 4 chord_root, 5 chord_tip, 6dihedral, 7thickness_over_chord)


#OBTAIN ZERO LIFT DRAG
wing_wet = Wing_wetted_area(wing_geo[4],wing_geo[5], d_fuselage, wing_geo[3], surface_area,0)
h_wet = H_tail_wetted_area(root_chord_h, taper_ratio_h, span_h)
v_wet = V_tail_wetted_area(root_chord_v, taper_ratio_v, span_v)
fus_wet = Fus_wetted_area(d_fuselage, l_nosecone, l_cabin, l_tailcone)
CD0 = Zero_Lift_Drag_est(surface_area, wing_wet, h_wet, v_wet, fus_wet)

#CD0 = 0.010 
Angle = 0.0 #Incidence angle

chords= [wing_geo[4], wing_geo[5]]
Ainc = [0.0, 0.0]
Nspanwise = [0, 0]
Sspace = [0, 0]
x_loc_LE = [0, wing_geo[3]/2*np.tan(wing_geo[1])]
y_loc_LE = [0, wing_geo[3]]
z_loc_LE = [0,wing_geo[3]*np.tan(wing_geo[6])] #change for dihedral angle low/high wing


with open("conv_wing.avl", "w") as text_file:
    print("Wing" +"\n"
        "#Mach" +"\n" + 
        str(0.70) +"\n"
        "#IYsym IZsym Zsym" +"\n"
        "0  0  0" +"\n" 
        "#Sref  Cref  Bref" +"\n" +
        str(surface_area), str(round(wing_geo[4],3)), str(round(wing_geo[3])), "\n"          
        "#Xref Yref Zref" +"\n"
        "0.3  0.0  0.0806" + "\n"
        "#CDcp" + "\n" +
        str(round(CD0,3)), "\n"
        "\n" + "SURFACE" +"\n" +
        "Wing", "\n"
        "8  1.0  12  -2.0"+"\n"
        "YDUPLICATE"+"\n" +
        str(0.0), "\n" +
        "ANGLE"+"\n"+
        str(Angle), file=text_file)
    for i in range(2):
        print("SECTION", file=text_file)            
        print(round(x_loc_LE[i],3),round(y_loc_LE[i],3),round(z_loc_LE[i],3),round(chords[i],3),Ainc[i],Nspanwise[i], Sspace[i], file=text_file)        


p = subprocess.Popen(r"C:\Users\floyd\Desktop\avl.exe", stdin=subprocess.PIPE, universal_newlines=True)
p.communicate(os.linesep.join(["load", "testwing", "oper", "W", "test"]))            

#
#p = subprocess.Popen(r"C:\Users\floyd\Desktop\avl.exe", stdin=subprocess.PIPE, universal_newlines=True)
#p.communicate(os.linesep.join(["load", "testwing", "oper","a a 2", "x","ft", "results"]))   