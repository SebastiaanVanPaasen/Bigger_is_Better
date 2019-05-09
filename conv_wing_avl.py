# -*- coding: utf-8 -*-
"""
Created on Thu May  9 13:38:52 2019

@author: floyd
"""
import numpy as np
from inputs import *
from planform import wing_parameters


surface_area = 200  #Once this is included in inputs file then this should be deleted 
wing_geo  = wing_parameters(M_cruise, CL_cruise, surface_area, A) #planform parameters: 
#0 quarter_chord_sweep, 1 leading_edge_sweep, 2 taper_ratio, 3 span, 4 chord_root, 5 chord_tip, 6dihedral, 7thickness_over_chord)

CD0 = 0.010 
Angle = 0.0 #Incidence angle

chords= [wing_geo[4], wing_geo[5]]
Ainc = [0.0, 0.0]
Nspanwise = [0, 0]
Sspace = [0, 0]
x_loc_LE = [0, wing_geo[3]/2*np.tan(wing_geo[1])]
y_loc_LE = [0, wing_geo[3]]
z_loc_LE = [0,wing_geo[3]*np.tan(wing_geo[6])] #change for dihedral angle low/high wing


with open("conv_wing.avl", "w") as text_file:
    print("Boxed Wing" +"\n"
        "#Mach" +"\n" + 
        str(M_cruise) +"\n"
        "#IYsym IZsym Zsym" +"\n"
        "0  0  0" +"\n" 
        "#Sref  Cref  Bref" +"\n" +
        str(surface_area), str(wing_geo[4]), str(wing_geo[3]), "\n"          
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
            