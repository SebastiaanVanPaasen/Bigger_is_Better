# -*- coding: utf-8 -*-
"""
Created on Thu May  9 13:38:52 2019

@author: floyd
"""
import numpy as np
from inputs import *
from planform import wing_parameters
from drag_estimation import *
from subprocess import Popen, DEVNULL, PIPE
import os

"""
Comments about this file:
    avl.exe, mach0.7.run and the python file should be in the same file location. 
    Also airfoil files should be in the same file location
    The filepath in run_avl function should be changed to your filepath. 
"""


def make_avl_file(root_chord, tip_chord, span, LE_sweep, dihedral, S, CD_0, M_cruise, tailarm_h, span_h, \
                  root_chord_h, quarter_chord_sweep_h, taper_ratio_h, tailarm_v, span_v, root_chord_v, \
                  quarter_chord_sweep_v, taper_ratio_v, wingtype):
    Angle = 0  # Incidence angle
    chords = [root_chord, tip_chord]
    Ainc = [0.0, 0.0]
    Nspanwise = [0, 0]
    Sspace = [0, 0]
    x_loc_LE = [0, span / 2 * np.tan(LE_sweep)]
    y_loc_LE = [0, span / 2]
    z_loc_LE = [0, span / 2 * np.tan(dihedral)]  # change for dihedral angle low/high wing

    dx = 0.25 * root_chord_h * (1 - taper_ratio_h) + span_h / 2 * np.tan(quarter_chord_sweep_h)
    y_loc_LE_h = [0, span_h / 2]
    if wingtype == 1:
        dx_sweep = 0.25 * root_chord_v * (1 - taper_ratio_v) + span_v / 2 * np.tan(quarter_chord_sweep_v)
        z_loc_LE_h = [span_v / 2, span_v / 2]
        x_loc_LE_h = [tailarm_h + dx_sweep, tailarm_h + dx + dx_sweep]
    else:
        z_loc_LE_h = [0, 0]
        x_loc_LE_h = [tailarm_h, tailarm_h + dx]
    chords_h = [root_chord_h, root_chord_h * taper_ratio_h]
    Angle_h = 0.0  # Incidence angle horizontal stabilizer

    dx = 0.25 * root_chord_v * (1 - taper_ratio_v) + span_v / 2 * np.tan(quarter_chord_sweep_v)
    x_loc_LE_v = [tailarm_v, tailarm_v + dx]
    y_loc_LE_v = [0, 0]
    z_loc_LE_v = [0, span_v / 2]
    chords_v = [root_chord_v, root_chord_v * taper_ratio_v]
    Angle_v = 0.0
    with open("conv_wing.avl", "w") as text_file:
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
              "Wing", "\n"
                      "8  1.0  12  -2.0" + "\n"
                                           "YDUPLICATE" + "\n" +
              str(0.0), "\n" +
              "ANGLE" + "\n" +
              str(Angle), file=text_file)
        for i in range(2):
            print("SECTION", file=text_file)
            print(round(x_loc_LE[i], 3), round(y_loc_LE[i], 3), round(z_loc_LE[i], 3), round(chords[i], 3), Ainc[i],
                  Nspanwise[i], Sspace[i], file=text_file)
        print("AFILE" + "\n"
                        "n2414.dat.txt", file=text_file)

        # HORIZONTAL TAIL surface
        print("\n" + "SURFACE" + "\n" +
              "Stab", "\n"
                      "6 1.0 15 -1.1" + "\n"
                                        "YDUPLICATE" + "\n" +
              str(0.0), "\n" +
              "ANGLE" + "\n" +
              str(Angle_h), file=text_file)
        for i in range(2):
            print("SECTION", file=text_file)
            print(round(x_loc_LE_h[i], 3), round(y_loc_LE_h[i], 3), round(z_loc_LE_h[i], 3), round(chords_h[i], 3),
                  Ainc[i], Nspanwise[i], Sspace[i], file=text_file)
        print("AFILE" + "\n"
                        "n2414.dat.txt", file=text_file)

        # Vertical tail surface
        print("\n" + "SURFACE" + "\n" +
              "FIN", "\n"
                     "7 1.0 11 1.0" + "\n"
                                      "YDUPLICATE" + "\n" +
              str(0.0), "\n" +
              "ANGLE" + "\n" +
              str(Angle_v), file=text_file)
        for i in range(2):
            print("SECTION", file=text_file)
            print(round(x_loc_LE_v[i], 3), round(y_loc_LE_v[i], 3), round(z_loc_LE_v[i], 3), round(chords_v[i], 3),
                  Ainc[i], Nspanwise[i], Sspace[i], file=text_file)
        print("AFILE" + "\n"
                        "n0010.dat.txt", file=text_file)
    return


make_avl_file(12.22, 12.22 * 0.149, 60.90, 0.56, 0.03, 427.8, 0.015, 0.7, 27.90, 21.35, 5, np.radians(35), 0.5, 26.6,
              18.48, 4, np.radians(46), 0.5, 1)


def run_avl(cl_cruise):
    p = Popen(r"C:\Users\sebas\OneDrive\Documents\DSE\Bigger_is_Better\avl.exe", stdin=PIPE,
              stdout=DEVNULL, universal_newlines=True)
    set_cl_cruise = "a c " + str(cl_cruise)
    p.communicate(
        os.linesep.join(["load", "conv_wing", "case", "mach0.7", "oper", set_cl_cruise, "x", "ft", "endresult"]))
    lines = [line.rstrip('\n') for line in open('endresult')]
    alpha = float(lines[15].split()[2])
    CD = float(lines[24].split()[2])
    e = float(lines[27].split()[5])
    print("Alpha is ", alpha)
    print("CD is ", CD)
    print("eff factor is ", e)
    os.remove("endresult")
    return e


run_avl(0.7)


# path for floyd: "C:\Users\floyd\Desktop\avl.exe"
# path for sebas: "C:\Users\sebas\OneDrive\Documents\DSE\Bigger_is_Better\avl.exe"
def find_clalpha():
    alpha_range = [0, 5]
    CL_range = []
    for j in range(len(alpha_range)):
        p = Popen(r"C:\Users\sebas\OneDrive\Documents\DSE\Bigger_is_Better\avl.exe", stdin=PIPE,
                  stdout=DEVNULL, universal_newlines=True)
        set_alpha = "a a " + str(alpha_range[j])
        p.communicate(
            os.linesep.join(["load", "conv_wing", "case", "mach0.7", "oper", set_alpha, "x", "ft", "endresult"]))
        lines = [line.rstrip('\n') for line in open('endresult')]

        CL = float(lines[23].split()[2])
        CL_range.append(CL)
        os.remove("endresult")
    CL_alpha = round((CL_range[1] - CL_range[0]) / (alpha_range[1] - alpha_range[0]), 3)
    return CL_alpha


print("CLa is ", find_clalpha())
