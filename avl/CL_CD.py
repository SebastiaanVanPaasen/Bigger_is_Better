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
from xfoil import xfoil_shit
from matplotlib import pyplot as plt
#inputs from wingplanform
##reference 777
#S = 427.80
#span = 60.900
#MAC = 8.750
#AR = 8.67
#taper = 0.149
#qc_sweep = np.radians(31.60)
#dihedral = 0
#Cr = (2*S)/((1+taper)*span)
#Ct = Cr*taper
#CD_0 = 0.01
#M_cruise = 0.0
#spanwise_discretize_points = 20
#chordwise_discretize_point = 8
#LE_sweep = np.arctan((0.25*Cr + span/2*np.tan(qc_sweep) - 0.25*Ct)/(span/2))


#AERO CONCEPT USED FOR VERIFICATION
#S = 223
#span = 55.87
#AR = 14
#taper = 0.31
#dihedral = 0.04
#Cr = 6.1
#Ct = 1.89
#CD_0 = 0.0262
#M_cruise = 0.75
#spanwise_discretize_points = 16
#chordwise_discretize_point = 8
#LE_sweep = np.radians(27.69)

# HIGH BYPASS ENGINE USED FOR NO SWEEP VERIFICATION
#S = 291.42
#span = 48.28
#AR = 8
#taper = 0.4
#dihedral = 0.09
#Cr = 8.62
#Ct = 3.45
#CD_0 = 0.0219
#M_cruise = 0.7
#spanwise_discretize_points = 16
#chordwise_discretize_point = 8
#LE_sweep = np.radians(3.07)
S = np.array([222.99,193.23,186.97,291.42,208.79])
b = np.array([55.87,41.70,39.87,48.28,61.30])
A = np.array([14,9,8.5,8,18])
taper = np.array([0.31,0.4,0.41,0.4,0.4])
dihedral = np.array([np.radians(2.4),np.radians(5),np.radians(2.4),np.radians(1.0),np.radians(1.0)])
Cr = np.array([6.10,6.62,7.16,8.62,4.87])
Ct= np.array([1.89,2.65,2.22,3.45,1.95])
CD_0 = np.array([0.0262,0.0202,0.0205,0.0219,0.0276])
M_cruise = np.array([0.75,0.7,0.75,0.7,0.7])
spanwise_discretize_points = np.array([16,16,16,16,16])
chordwise_discretize_point = np.array([8,8,8,8,8])
LE_sweep = np.array([np.radians(27.69),np.radians(2.73),np.radians(28.77),np.radians(3.07),np.radians(1.36)])
#figurenumber = [1,2,3,4,5]
#%%
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
                  "-1.54 0.002 0.69 0.0054 1.8 0.002", file=text_file)
        print("AFILE" + "\n"
                        "n3414.dat", file=text_file)
#make_avl_file(Cr,Ct,span,LE_sweep,dihedral,S,CD_0, M_cruise,spanwise_discretize_points,chordwise_discretize_point)

for j in range(len(S)):
    make_avl_file(Cr[j],Ct[j],b[j],LE_sweep[j],dihedral[j],S[j],CD_0[j], M_cruise[j],spanwise_discretize_points[j],chordwise_discretize_point[j])
    step = .05
    CL = np.arange(-2,2.2+step,step)
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
    plt.figure(1)
    plt.plot(CD,CL)
    plt.grid()    
    plt.xlabel("CD [-]")
    plt.ylabel("CL [-]")
    plt.show()
plt.legend(["AERO", "DD2E","DD4E", "HBPE","STRW"])
#%%
#CL2 = np.arange(-2,2.2+step,step)
##e = 4.61*(1 - 0.045*AR**0.68) *(np.cos(LE_sweep))**0.15 - 3.1
#e = 1/(1.05 + 0.007* np.pi *AR)
#print(e)
#k = 1/(np.pi*AR*e)
#CD2 = CD_0 + k*CL**2
#plt.plot(CD2,CL2)
#plt.legend(["AVL","Obert Approx."])
#
#Re = 7.5
#data = xfoil_shit(7.5)
#plt.plot(data[1],data[0])


