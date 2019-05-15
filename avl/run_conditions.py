# -*- coding: utf-8 -*-
"""
Created on Mon May 13 11:13:44 2019

@author: floyd
"""


def define_run_condition(Mach, CD_0):
    fname = "mach" + str(Mach) + ".run"
    with open(fname, "w") as text_file:
        print("\n---------------------------------------------\n"
              "\n"
              "Run case  1:   -unnamed-\n"
              "\n"
              "alpha        ->  alpha       =   0.00000\n"
              "beta         ->  beta        =   0.00000\n"
              "pb/2V        ->  pb/2V       =   0.00000\n"
              "qc/2V        ->  qc/2V       =   0.00000\n"
              "rb/2V        ->  rb/2V       =   0.00000\n"
              "\n"
              "alpha     =   0.00000     deg\n"
              "beta      =   0.00000     deg\n"
              "pb/2V     =   0.00000\n"
              "qc/2V     =   0.00000\n"
              "rb/2V     =   0.00000\n"
              "CL        =   0.00000\n"
              "CDo       =   " + str(CD_0) + "\n"
                                             "bank      =   0.00000     deg\n"
                                             "elevation =   0.00000     deg\n"
                                             "heading   =   0.00000     deg\n"
                                             "Mach      =   " + str(Mach) + "\n"
                                                                            "velocity  =   0.00000     Lunit/Tunit\n"
                                                                            "density   =   1.00000     Munit/Lunit^3\n"
                                                                            "grav.acc. =   1.00000     Lunit/Tunit^2\n"
                                                                            "turn_rad. =   0.00000     Lunit\n"
                                                                            "load_fac. =   0.00000\n"
                                                                            "X_cg      =   1.00000     Lunit\n"
                                                                            "Y_cg      =   0.00000     Lunit\n"
                                                                            "Z_cg      =   0.00000     Lunit\n"
                                                                            "mass      =   1.00000     Munit\n"
                                                                            "Ixx       =   1.00000     Munit-Lunit^2\n"
                                                                            "Iyy       =   1.00000     Munit-Lunit^2\n"
                                                                            "Izz       =   1.00000     Munit-Lunit^2\n"
                                                                            "Ixy       =   0.00000     Munit-Lunit^2\n"
                                                                            "Iyz       =   0.00000     Munit-Lunit^2\n"
                                                                            "Izx       =   0.00000     Munit-Lunit^2\n"
                                                                            "visc CL_a =   0.00000\n"
                                                                            "visc CL_u =   0.00000\n"
                                                                            "visc CM_a =   0.00000\n"
                                                                            "visc CM_u =   0.00000", file=text_file)
# define_run_condition(0.4,0.020)
