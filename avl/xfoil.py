# -*- coding: utf-8 -*-
"""
Created on Thu May 16 09:56:03 2019

@author: floyd
"""
import os
import subprocess
from matplotlib import pyplot as plt
from conv_wing_avl import * 
#from CL_CD import *
ROOT_DIR = os.path.dirname(os.path.abspath("XFOIL"))
Re = 7.5 #in millions
def xfoil_shit(Re):
    p = subprocess.Popen(str(ROOT_DIR) + "/xfoil.exe", stdin=subprocess.PIPE, universal_newlines=True)
    p.communicate(os.linesep.join(["load", "sc20616.dat", "pane","oper", "iter 100", "visc " + str(Re) + "e6", "seqp", "pacc", "polar2", "dump", "cseq -1  3  .1"]))
    lines = [line.rstrip('\n') for line in open('polar2')]
    CL = []
    CD = []
    alpha = []
    for i in range(12,len(lines)-1):
        CL.append(float(lines[i].split()[1]))
        CD.append(float(lines[i].split()[2])) 
        alpha.append(float(lines[i].split()[0]))
    os.remove("polar2")
    os.remove("dump")
    return(CL,CD,alpha)
#%%
data = xfoil_shit(31.5 /4.33*6.04)
#
##make_avl_file(12.22, 12.22 * 0.149, 60.90, 0.56, 0.03, 427.8, 0.015, 0.6, 27.90, 21.35, 5, np.radians(35), 0.5, 26.6,18.48, 4, np.radians(46), 0.5, " ", 18, 5)
##run_avl(0.5, 0.85, 0.020)
###print("CLa is ", find_clalpha(0.7, 0.020))
#
#
##plt.figure(10)
##Drag polar
#plt.subplot(221)
#plt.plot(data[1],data[0])
#plt.grid()
##%%
##CLalpha curves
#
##cruise conditions
#mach = 0.75
#CD_0 = 0.0262
#filename = "aero_avl"
#
##alpha_range for avl
#x = np.arange(-5,20,1)
##find zero lift AoA
CL_array = np.asarray(data[0])
alpha_arr = np.asarray(data[2])
index = np.where(CL_array == 0.0)
zero_lift_angle = alpha_arr[index]
print(zero_lift_angle)
#CLa =  find_clalpha(mach, CD_0, filename)
#y = CLa*(x-zero_lift_angle) 
#
##DATCOM
#
##A = 14
##chord_tip = 1.89
##chord_root = 6.1 
##span = 55.87
##LE_sweep = np.radians(27.69)
##beta = np.sqrt(1-mach**2)
##hc_sweep = np.arctan(
##        ((chord_tip / 2 + span / 2 * np.tan(LE_sweep)) - (chord_root / 2)) / (span / 2))
##CLa_DATCOM = ((2*np.pi*A)/(2+np.sqrt(4+(A*beta/0.95)**2*(1+(np.tan(hc_sweep))**2/beta**2))))/57.3
##print(CLa)
##print(CLa_DATCOM)
###percentage = (CLa- CLa_DATCOM)/CLa_DATCOM*100
##y2 = CLa_DATCOM*(x-zero_lift_angle) 
#
##plt.subplot(222)
##plt.plot(data[2],data[0])
#plt.plot(x,y)
#plt.plot(x,y2)
#plt.xlim((-10,25))
#plt.xlabel("Angle of Attack [deg]")
#plt.ylabel("CL [-]")
#plt.legend([ "AVL", "DATCOM"])
#
##y = np.arange(-1.5,1.7,0.01)
##x = 0.001*(y-0.7)**2 + 0.0054
##plt.plot(x,y)
##plt.plot(CL,CD)
#plt.grid()
##plt.axis('scaled')
