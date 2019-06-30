# -*- coding: utf-8 -*-
"""
Created on Sun Jun 30 09:05:42 2019

@author: mathi
"""
from matplotlib import pyplot as plt

strut_pos = [9, 11, 13, 15, 17, 19,21]
C_QS = [40713.11525, 31827.36436,24147.58158,17670.91814,13472.24446, 15500.25665, 20714.41559]
ALU =[41677.45441, 33149, 25905, 19933, 16304, 18977, 24865]
C_UD = [40043, 31006, 23114, 16379, 11879, 13318, 18201.1]

strut_pos_2 = [17,17.5,18,18.5,19]
C_UD_2 = [11879, 11352, 10881, 12054, 13252.36223]

plt.rcParams.update({'font.size': 20})        

plt.figure()
plt.plot(strut_pos, C_QS, label = " QI carbon fibre")
plt.plot(strut_pos, C_UD, label = " UD carbon fibre")
plt.plot(strut_pos, ALU, label = "Aluminium 2024")
plt.xlabel("Strut position along the halfspan [m]")
plt.ylabel("Comined stiffener and strut [m]")
ax = plt.gca()
ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

plt.legend()

plt.figure()
plt.plot(strut_pos_2, C_UD_2, label = "UD carbon fibre")
plt.xlabel("Strut position along the halfspan [m]")
plt.ylabel("Comined stiffener and strut [m]")
ax = plt.gca()
ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

plt.legend()
plt.show()