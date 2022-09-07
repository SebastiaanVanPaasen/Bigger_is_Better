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
A_weight = []
A_weight_cor = []

#To compute spectral SPL
for i in range(len(SPL)):
    if freq[i] < 1000:
        SPL_overall = 10.*np.log10(10.**(SPL[i]/10.) )#+ 10.**(SPL[i]/10.))
        SPL_coroverall = 10.*np.log10(10.**(SPL_cor[i]/10.))# + 10.**(SPL_cor[i]/10.))   
        
        dLa = -145.528 + 98.262*np.log10(freq[i]) - 19.509*np.log10(freq[i])**2 + 0.975*np.log10(freq[i]**3)
        A_weight_i = 10.*np.log10(10.**((SPL_overall + dLa)/10.))
        A_weight_icor = 10.*np.log10(10.**((SPL_coroverall + dLa)/10.))
        
    elif freq[i] > 1000:
        SPL_overall = 10.*np.log10(10.**(SPL[i]/10.) + 10.**(SPL_fan_list[i]/10.))
        SPL_coroverall = 10.*np.log10(10.**(SPL_cor[i]/10.)+ 10.**(SPL_fan_effects[i]/10.))   

        dLa = -145.528 + 98.262*np.log10(freq[i]) - 19.509*np.log10(freq[i])**2 + 0.975*np.log10(freq[i])**3
        A_weight_i = 10.*np.log10(10.**((SPL_overall + dLa)/10.))
        A_weight_icor = 10.*np.log10(10.**((SPL_coroverall + dLa)/10.))
        
    SPL_overall_list.append(SPL_overall)
    SPL_coroverall_list.append(SPL_coroverall)
    A_weight.append(A_weight_i)
    A_weight_cor.append(A_weight_icor)

#To compute the OASPL
SPL_sum = []
SPL_sum_cor = []
A_weight_sum = []
A_weight_corsum = []

for j in range(len(SPL_overall_list)):
    SPL_sum.append(10.**(SPL_overall_list[j]/10.))
    SPL_sum_cor.append(10.**(SPL_coroverall_list[j]/10.)) 
    
    A_weight_sum.append(10.**(A_weight[j]/10.)) 
    A_weight_corsum.append(10.**(A_weight_cor[j]/10.)) 
    
OASPL = 10.*np.log10(sum(SPL_sum)) 
OASPL_cor = 10.*np.log10(sum(SPL_sum_cor))
A_OASPL = 10.*np.log10(sum(A_weight_sum))
A_OASPL_cor = 10.*np.log10(sum(A_weight_corsum))

print ("OASPL model:", OASPL)
print ("OASPL corrected model:", OASPL_cor)
print("A weighted:", A_OASPL)
print ("A weighted corrected:", A_OASPL_cor)


plt.figure(3)  
plt.plot(freq,SPL_overall_list)
#plt.title("Spectral SPL levels overall aircraft") 
plt.grid(True)
plt.xlabel("1/3 Octave Band Centre Frequency [Hz]",fontsize ='x-large')
plt.ylabel("1/3 Octave Band SPL [dB]",fontsize ='x-large')
plt.xlim(50,10000)
#plt.ylim(50,80)
plt.xscale('log')

plt.figure(4)
 
plt.plot(freq,SPL_coroverall_list)
#plt.title("Spectral SPL corrected levels of overall aircraft") 
plt.grid(True)
plt.xlabel("1/3 Octave Band Centre Frequency [Hz]",fontsize ='x-large')
plt.ylabel("1/3 Octave Band SPL [dB]",fontsize ='x-large')
plt.xlim(50,10000)
plt.ylim(30,80)
plt.xscale('log')


#plt.figure(5) 
##plt.plot(freq,A_weight)
#plt.title("A weighted levels overall aircraft") 
#plt.grid(True)
#plt.xlabel("1/3 Octave Band Centre Frequency [Hz]",fontsize ='x-large')
#plt.ylabel("1/3 Octave Band La [dBA]",fontsize ='x-large')
#plt.xlim(50,10000)
##plt.ylim(40,120)
#plt.xscale('log')
#
#
#plt.figure(6) 
#plt.plot(freq,A_weight_cor)
##plt.title("A weighted levels corrected levels of overall aircraft") 
#plt.grid(True)
#plt.xlabel("1/3 Octave Band Centre Frequency [Hz]",fontsize ='x-large')
#plt.ylabel("1/3 Octave Band La [dBA]",fontsize ='x-large')
#plt.xlim(50,10000)
#plt.xscale('log')
#



plt.show()














    
    






