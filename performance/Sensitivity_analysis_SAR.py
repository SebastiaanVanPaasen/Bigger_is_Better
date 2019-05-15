# -*- coding: utf-8 -*-
"""
Created on Wed May 15 11:43:05 2019

@author: nikki
"""
import matplotlib.pyplot as plt

#---------------------SENSITIVITY ANALYSIS RESULTS---------------------------
"""Aspect ratio"""
A = 1.
A1 = 0.8            #decrease of 20% in A
A2 = 1.2            #increase of 20% in A

H = range(7000,12500,500)    #range of cruise altitudes

#Difference in min SAR due to decrease
minSAR_A1 = [10.20619665379791, 10.906268391494073, 11.530168098395135, 11.797703244608851, 11.803741179981559, 11.803270401059375, 11.804412818054905, 11.802160864343852, 11.803870777743739, 11.802601351980233, 11.804421041547872]
#Difference in V at min SAR due to decrease
V_minSAR_A1 = [0.0, 0.0, 2.0202020202020332, 5.3872053872053955, 5.2459016393442699, 5.7507987220447276, 6.1919504643962755, 6.0060060060060145, 5.7971014492753703, 5.5710306406685319, 5.8981233243967672]
#Overall difference in SAR due to decrease 
SAR_A1= [5.0403294502889064, 5.4904506565581936, 5.9737711190444882, 6.4909837381687385, 7.042423076320067, 7.628000392827504, 8.2471402811850982, 8.8987217982262887, 9.575954508606813, 10.448918480421542, 11.341707426159623]

#Difference in min SAR due to increase
minSAR_A2 = [-6.8041311025318825, -7.2708455943294155, -7.7509424650911622, -8.2416758479037249, -8.6319678292529307, -8.7127515793493533, -8.7133628039944799, -8.7139105923248366, -8.7133415816786748, -8.7129151588411258, -8.7118852152922202]
#Difference in V at min SAR due to increase
V_minSAR_A2 = [0.0, 0.0, 0.0, 0.0, -2.6229508196721349, -4.4728434504792487, -4.334365325077397, -4.2042042042041841, -4.6376811594202794, -4.4568245125348263, -4.2895442359249385]
#Overall difference in SAR due to increase 
SAR_A2 = [-3.3602196335259347, -3.6603004377054651, -3.9825140793629976, -4.3273224921124953, -4.6949487175467119, -5.0853335952183354, -5.4980935207900652, -5.9324811988175243, -6.3839696724045414, -6.9659456536143667, -7.561138284106419]

##plt.subplot(212)
#plt.plot(H,minSAR_A1,"r",label = "%change in min SAR ")
#plt.plot(H,V_minSAR_A1,"b",label = "%change in V at min SAR")
#plt.plot(H,SAR_A1,"g",label = "%average change in SAR")
#plt.title("Changes due to a 20% decrease in A")
#plt.xlabel("Altitude [m]")
#plt.ylabel("Change in %")
#plt.legend(loc = "lower right")
##plt.show()
#
##plt.subplot(211)
#plt.plot(H,minSAR_A2,"r",label = "%change in min SAR ")
#plt.plot(H,V_minSAR_A2,"b",label = "%change in V at min SAR")
#plt.plot(H,SAR_A2,"g",label = "%average change in SAR")
#plt.title("Changes due to a 20% increase in A")
#plt.xlabel("Altitude [m]")
#plt.ylabel("Change in %")
#plt.legend()
#plt.show()


"""Zero lift drag CD0"""
CD0 = 1.
CD01 = 0.8            #decrease of 20% in CD0
CD02 = 1.2            #increase of 20% in CD0

H = range(7000,12500,500)    #range of cruise altitudes

#Difference in min SAR due to decrease
minSAR_CD01 =[-11.83504267696167, -11.27498528680472, -10.775865521283894, -10.561837404312907, -10.557007056014752, -10.557383679152501, -10.556469745556115, -10.558271308524954, -10.556903377805019, -10.55791891841581, -10.556463166761681]
#Difference in V at min SAR due to decrease
V_minSAR_CD01 = [0.0, 0.0, 2.0202020202020332, 5.3872053872053955, 5.2459016393442699, 5.7507987220447276, 6.1919504643962755, 6.0060060060060145, 5.7971014492753703, 5.5710306406685319, 5.8981233243967672]
#Overall difference in SAR due to decrease 
SAR_CD01= [-15.967736439768869, -15.607639474753439, -15.220983104764409, -14.807213009465011, -14.366061538943946, -13.897599685737994, -13.402287775051924, -12.88102256141897, -12.33923639311455, -11.640865215662771, -10.926634059072301]

#Difference in min SAR due to increase
minSAR_CD02 = [11.835042676961731, 11.274985286804732, 10.698869041890635, 10.109988982515533, 9.6416386048965084, 9.5446981047807764, 9.5439646352065779, 9.5433072892101904, 9.5439901019855782, 9.5445018093906739, 9.5457377416492957]
#Difference in V at min SAR due to increase
V_minSAR_CD02 = [0.0, 0.0, 0.0, 0.0, -2.6229508196721349, -4.4728434504792487, -4.334365325077397, -4.2042042042041841, -4.6376811594202794, -4.4568245125348263, -4.2895442359249385]
#Overall difference in SAR due to increase 
SAR_CD02 = [15.967736439768869, 15.60763947475345, 15.220983104764407, 14.807213009465004, 14.366061538943944, 13.897599685738001, 13.402287775051931, 12.881022561418972, 12.339236393114552, 11.640865215662764, 10.926634059072299]

##plt.subplot(212)
#plt.plot(H,minSAR_CD01,"r",label = "%change in min SAR ")
#plt.plot(H,V_minSAR_CD01,"b",label = "%change in V at min SAR")
#plt.plot(H,SAR_CD01,"g",label = "%average change in SAR")
#plt.title("Changes due to a 20% decrease in CD0")
#plt.xlabel("Altitude [m]")
#plt.ylabel("Change in %")
#plt.legend(loc = "upper left")
#plt.show()
#
##plt.subplot(211)
#plt.plot(H,minSAR_CD02,"r",label = "%change in min SAR ")
#plt.plot(H,V_minSAR_CD02,"b",label = "%change in V at min SAR")
#plt.plot(H,SAR_CD02,"g",label = "%average change in SAR")
#plt.title("Changes due to a 20% increase in CD0")
#plt.xlabel("Altitude [m]")
#plt.ylabel("Change in %")
#plt.legend()
#plt.show()

