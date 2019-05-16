# -*- coding: utf-8 -*-
"""
Created on Wed May 15 20:25:01 2019

@author: floyd

find CL_CD polar of the wing 
"""
import numpy as np
import os
import subprocess
from run_conditions import define_run_condition
ROOT_DIR = os.path.dirname(os.path.abspath("CL_CD"))
from matplotlib import pyplot as plt
#inputs from wingplanform
##reference 777
S = 427.80
span = 60.90
MAC = 8.75
AR = 8.67
taper = 0.149
qc_sweep = np.radians(31.60)
dihedral = 0
Cr = (2*S)/((1+taper)*span)
Ct = Cr*taper
CD_0 = 0.01
M_cruise = 0.0
spanwise_discretize_points = 20
chordwise_discretize_point = 8
LE_sweep = np.arctan((0.25*Cr + span/2*np.tan(qc_sweep) - 0.25*Ct)/(span/2))

def make_avl_file(root_chord, tip_chord, span, LE_sweep, dihedral, S, CD_0, M_cruise, spanwise_discretize_points,
                  chordwise_discretize_point):  # for wingtype: input "T-tail" in case
    # of T-tail else random input
    Angle = 0.0  # Incidence angle
    chords = [root_chord, tip_chord]
    Ainc = [0.0, 0.0]
    Nspanwise = [0, 0]
    Sspace = [0, 0]
    x_loc_LE = [0, span / 2 * np.tan(LE_sweep)]
    y_loc_LE = [0, span / 2]
    z_loc_LE = [0, span / 2 * np.tan(dihedral)]  # change for dihedral angle low/high wing

    with open("CL_CD.avl", "w") as text_file:
        print("Wing" + "\n"
                       "#Mach" + "\n" +
              str(M_cruise) + "\n"
                              "#IYsym IZsym Zsym" + "\n"
                                                    "0  0  0" + "\n"
                                                                "#Sref  Cref  Bref" + "\n" +
              str(S), str(round(root_chord, 3)), str(round(span)), "\n"
                                                                   "#Xref Yref Zref" + "\n"
                                                                                       "0  0.0  0" + "\n"
                                                                                                     "#CDcp" + "\n" +
              str(round(CD_0, 3)), "\n"
                                   "\n" + "SURFACE" + "\n" +
              "Wing", "\n" +
              str(chordwise_discretize_point), "1.0 " + str(spanwise_discretize_points), "-2.0" + "\n"
                                                                                                  "YDUPLICATE" + "\n" +
              str(0.0), "\n" +
              "ANGLE" + "\n" +
              str(Angle), file=text_file)
        for i in range(2):
            print("SECTION", file=text_file)
            print(round(x_loc_LE[i], 3), round(y_loc_LE[i], 3), round(z_loc_LE[i], 3), round(chords[i], 3), Ainc[i],
                  Nspanwise[i], Sspace[i], file=text_file)
            print("CDCL" + "\n"
                  "-0.5  1.200  0.0  1.00  1.5  1.200", file=text_file)
        print("AFILE" + "\n"
                        "n3414.dat", file=text_file)
make_avl_file(Cr,Ct,span,LE_sweep,dihedral,S,CD_0, M_cruise,spanwise_discretize_points,chordwise_discretize_point)


step = .05
CL = np.arange(-2,2+step,step)
CD = np.zeros(len(CL))
for i in range(len(CL)):
    define_run_condition(M_cruise,CD_0)
    p = subprocess.Popen(str(ROOT_DIR) + "/avl.exe", stdin=subprocess.PIPE, universal_newlines=True)
    set_CL = "a c " + str(CL[i])
    p.communicate(os.linesep.join(["load", "CL_CD","case","mach" + str(M_cruise) + ".run", "oper", set_CL, "x","ft", "endresult"]))          
    lines = [line.rstrip('\n') for line in open('endresult')]
    CD[i] = float(lines[24].split()[2])
    
    os.remove("endresult")
    os.remove("mach" + str(M_cruise) + ".run")
plt.plot(CD,CL)
plt.grid()    
plt.xlabel("CD")
plt.ylabel("CL")
