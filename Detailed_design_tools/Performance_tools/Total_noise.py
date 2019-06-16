# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 16:51:44 2019

@author: nikki
"""

from Airframe_noise import SPL, SPL_cor, freq
from Engine_noise import SPL_fan_list, SPL_fan_effects
import numpy as np
import matplotlib.pyplot as plt


"""Combined results"""
SPL_overall_list = []
SPL_coroverall_list = []

#To compute spectral SPL
for i in range(len(SPL)):
    SPL_overall = 10.*np.log10(10.**(SPL[i]/10.) + 10.**(SPL_fan_list[i]/10.))
    SPL_coroverall = 10.*np.log10(10.**(SPL_cor[i]/10.) + 10.**(SPL_fan_effects[i]/10.))   

    SPL_overall_list.append(SPL_overall)
    SPL_coroverall_list.append(SPL_coroverall)

#To compute the OASPL
SPL_sum = []
SPL_sum_cor = []


for j in range(len(SPL_overall_list)):
    SPL_sum.append(10.**(SPL_overall_list[j]/10.))
    SPL_sum_cor.append(10.**(SPL_coroverall_list[j]/10.)) 
    
OASPL = 10.*np.log10(sum(SPL_sum)) 
OASPL_cor = 10.*np.log10(sum(SPL_sum_cor))

print ("OASPL model:", OASPL)
print ("OASPL corrected model:", OASPL_cor)


plt.figure(3)
plt.subplot(121)  
plt.plot(freq,SPL_overall_list)
plt.title("Spectral SPL levels overall aircraft") 
plt.grid(True)
plt.xlabel("1/3 Octave Band Centre Frequency [Hz]")
plt.ylabel("1/3 Octave Band SPL [dB]")
plt.xlim(50,10000)
plt.ylim(40,120)
plt.xscale('log')


plt.subplot(122)  
plt.plot(freq,SPL_coroverall_list)
plt.title("Spectral SPL corrected levels of overall aircraft") 
plt.grid(True)
plt.xlabel("1/3 Octave Band Centre Frequency [Hz]")
plt.ylabel("1/3 Octave Band SPL [dB]")
plt.xlim(50,10000)
plt.xscale('log')


plt.show()














    
    






