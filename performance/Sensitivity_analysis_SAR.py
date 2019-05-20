# -*- coding: utf-8 -*-
"""
Created on Wed May 15 11:43:05 2019

@author: nikki
"""
import matplotlib.pyplot as plt
#---------------------------SENSITIVITY ANALYSIS CRUISE ALTITUDE---------------------
H = range(7500,12500,500)    #range of cruise altitudes

#Difference in min SAR due to decrease
minSAR_H = [-0.88924196295702018, -0.57301721587421872, -0.23847334349236435, -0.0051125378857262462, 0.00010912561194382052, -0.00052209251888409741, 0.00059026811977999801, -0.00061617253561938566, 0.00023998881178734962, -0.00062762048901803468]
#Difference in V at min SAR due to decrease
V_minSAR_H = [0.0, 0.0, 0.0, 2.6936026936026947, 2.6229508196721318, 3.1948881789137387, 3.0959752321981431, 3.6036036036036041, 4.0579710144927548, 3.8997214484679672]


plt.plot(H,minSAR_H,"r",label = "%change in min SAR ")
plt.plot(H,V_minSAR_H,"b",label = "%change in V at min SAR")
#plt.title("Changes with respect to the previous values due to step changes in altitude")
plt.xlabel("Altitude [m]")
plt.ylabel("Change in % ")
plt.legend()
plt.show()



#---------------------SENSITIVITY ANALYSIS DESIGN PARAMETERS---------------------------
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

plt.subplot(121)
plt.plot(H,minSAR_A1,"r",label = "%change in min SAR ")
plt.plot(H,V_minSAR_A1,"b",label = "%change in V at min SAR")
plt.plot(H,SAR_A1,"g",label = "%average change in SAR")
#plt.title("Changes due to a 20% decrease in A")
plt.xlabel("Altitude [m]")
plt.ylabel("Change in %")
plt.legend(loc = "lower right")


plt.subplot(122)
plt.plot(H,minSAR_A2,"r",label = "%change in min SAR ")
plt.plot(H,V_minSAR_A2,"b",label = "%change in V at min SAR")
plt.plot(H,SAR_A2,"g",label = "%average change in SAR")
#plt.title("Changes due to a 20% increase in A")
plt.xlabel("Altitude [m]")
plt.ylabel("Change in %")
plt.legend()
plt.show()

"""Surface Area"""
S = 1.
S1 = 0.8            #decrease of 20% in A
S2 = 1.2            #increase of 20% in A

H = range(7000,12500,500)    #range of cruise altitudes

#Difference in min SAR due to decrease
minSAR_S1 = [-1.6978101578155467, -0.81572803622943013, -0.24424964177565084, -0.0049779927413986878, -9.3022150099715372e-07, -0.00094055144435837331, 0.00024796977405798466, -0.0010646330086639892, -0.00048014350284827136, -0.00054249428787121116, -1.8156931194226044e-05]
#Difference in V at min SAR due to decrease
V_minSAR_S1 = [2.0202020202020332, 4.7138047138047305, 8.080808080808092, 11.447811447811455, 11.803278688524589, 12.140575079872194, 11.764705882352949, 12.012012012012029, 11.594202898550725, 11.69916434540391, 11.79624664879355]
#Overall difference in SAR due to decrease
SAR_S1= [-10.927406989479964, -10.117188818195247, -9.2472119857199111, -8.3162292712962653, -7.3236384626238795, -6.2695992929104944, -5.155147493866826, -3.9823007631926806, -2.7632818845077325, -1.1919467352412285, 0.41507336708732795]

#Difference in min SAR due to increase
minSAR_S2 =[5.0309115744297621, 4.0041396924753521, 2.9479265767994844, 1.8683131346117952, 0.89135898180397566, 0.26410772949638595, 0.0057994623428088447, -0.00092551971025761517, 0.00012639199792893114, -0.00072180017626833357, 0.00086364110749315777]
#Difference in V at min SAR due to increase
V_minSAR_S2 = [0.0, 0.0, 0.0, 0.0, -2.6229508196721349, -5.1118210862619886, -8.0495356037151726, -8.4084084084084036, -8.6956521739130466, -8.9136490250696347, -8.579088471849877]
#Overall difference in SAR due to increase
SAR_S2 = [12.607516806242936, 11.947339037047987, 11.238469025401418, 10.479890517352514, 9.6711128213972373, 8.8122660905196657, 7.9041942542618688, 6.9485413626014489, 5.9552667207100134, 4.674919562048407, 3.3654957749658867]

plt.subplot(121)
plt.plot(H,minSAR_S1,"r",label = "%change in min SAR ")
plt.plot(H,V_minSAR_S1,"b",label = "%change in V at min SAR")
plt.plot(H,SAR_S1,"g",label = "%average change in SAR")
#plt.title("Changes due to a 20% decrease in S")
plt.xlabel("Altitude [m]")
plt.ylabel("Change in %")
plt.legend(loc = "lower right")


plt.subplot(122)
plt.plot(H,minSAR_S2,"r",label = "%change in min SAR ")
plt.plot(H,V_minSAR_S2,"b",label = "%change in V at min SAR")
plt.plot(H,SAR_S2,"g",label = "%average change in SAR")
#plt.title("Changes due to a 20% increase in S")
plt.xlabel("Altitude [m]")
plt.ylabel("Change in %")
plt.legend()


plt.show()


"""Cruise weight"""
W = 1.
W1 = 0.8            #decrease of 20% in A
W2 = 1.2            #increase of 20% in A

H = range(7000,12500,500)    #range of cruise altitudes

#Difference in min SAR due to decrease
minSAR_W1 = [-17.860844144146284, -19.085969685114641, -20.346223970864283, -21.634399100747327, -22.851172216778277, -23.807264316773704, -24.489932124150329, -24.890885266010677, -24.999567554265766, -24.999805943051424, -25.000017837264693]
#Difference in V at min SAR due to decrease
V_minSAR_W1 = [0.0, 0.0, 0.0, 0.0, -2.6229508196721349, -5.1118210862619886, -8.0495356037151726, -10.810810810810811, -13.333333333333345, -13.370473537604443, -13.404825737265428]
#Overall difference in SAR due to decrease
SAR_W1=[-8.8205765380055912, -9.6082886489768438, -10.454099458327864, -11.359221541795307, -12.324240383560126, -13.349000687448145, -14.432495492073922, -15.572763146896008, -16.757920390061937, -18.285607340737705, -19.847987995779356]

#Difference in min SAR due to increase
minSAR_W2 = [22.877737302730555, 23.980339954713138, 24.694687947780398, 24.993777509073212, 24.999998837223103, 24.998824310694573, 25.00030996221756, 24.998669208739141, 24.999399820621427, 24.999321882140162, 24.99997730383599]
#Difference in V at min SAR due to increase
V_minSAR_W2 = [2.0202020202020332, 4.7138047138047305, 8.080808080808092, 11.447811447811455, 11.803278688524589, 12.140575079872194, 11.764705882352949, 12.012012012012029, 11.594202898550725, 11.69916434540391, 11.79624664879355]
#Overall difference in SAR due to increase
SAR_W2 = [11.340741263150033, 12.353513977255927, 13.440985017850098, 14.604713410879658, 15.845451921720132, 17.163000883861866, 18.556065632666446, 20.022124046009129, 21.545897644365308, 23.510066580948454, 25.518841708859146]

# plt.subplot(121)
plt.figure(1)
plt.plot(H,minSAR_W1,"r",label = "%change in min SAR ")
#plt.plot(H,V_minSAR_W1,"b",label = "%change in V at min SAR")
plt.plot(H,SAR_W1,"g",label = "%average change in SAR")
#plt.title("Changes due to a 20% decrease in W")
plt.xlabel("Altitude [m]")
plt.ylabel("Change in %")
plt.legend(loc = "upper right")

# plt.subplot(122)
plt.figure(2)
plt.plot(H,minSAR_W2,"r",label = "%change in min SAR ")
#plt.plot(H,V_minSAR_W2,"b",label = "%change in V at min SAR")
plt.plot(H,SAR_W2,"g",label = "%average change in SAR")
#plt.title("Changes due to a 20% increase in W")
plt.xlabel("Altitude [m]")
plt.ylabel("Change in %")
plt.legend(loc = "lower right")
plt.show()

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

plt.subplot(121)
plt.plot(H,minSAR_CD01,"r",label = "%change in min SAR ")
plt.plot(H,V_minSAR_CD01,"b",label = "%change in V at min SAR")
plt.plot(H,SAR_CD01,"g",label = "%average change in SAR")
plt.title("Changes due to a 20% decrease in CD0")
plt.xlabel("Altitude [m]")
plt.ylabel("Change in %")
plt.legend(loc = "upper right")


plt.subplot(122)
plt.plot(H,minSAR_CD02,"r",label = "%change in min SAR ")
plt.plot(H,V_minSAR_CD02,"b",label = "%change in V at min SAR")
plt.plot(H,SAR_CD02,"g",label = "%average change in SAR")
plt.title("Changes due to a 20% increase in CD0")
plt.xlabel("Altitude [m]")
plt.ylabel("Change in %")
plt.legend()

plt.show()




"""Nominal Thrust specific fuel consumption Ct0"""
#To not forget this one since the plot is not interesting but the results are!!!
# print
print("changing Ct0 has a uniform effect on the parameters")
print("Decreasing Ct0 by 20% decreases the min SAR and overal SAR by 20%")
print("The V at which min SAR is obtained is not affected")
# print()
print("Increasing Ct0 by 20% increases the min SAR and overal SAR by 20%")
print("Again the V at which min SAR is obtained is not affected")




Ct0 = 1.
Ct01 = 0.8            #decrease of 20% in CD0
Ct02 = 1.2            #increase of 20% in CD0

H = range(7000,12500,500)    #range of cruise altitudes

#Difference in min SAR due to decrease
minSAR_Ct01 =[-19.999999999999996, -19.999999999999982, -19.999999999999989, -19.999999999999979, -19.999999999999979, -19.999999999999982, -19.999999999999986, -19.999999999999993, -19.999999999999979, -19.999999999999979, -19.999999999999986]
#Difference in V at min SAR due to decrease
V_minSAR_Ct01 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#Overall difference in SAR due to decrease
SAR_Ct01= [-19.999999999999996, -19.999999999999996, -19.999999999999996, -19.999999999999996, -20.0, -19.999999999999996, -19.999999999999996, -19.999999999999996, -19.999999999999996, -19.999999999999996, -19.999999999999996]

#Difference in min SAR due to increase
minSAR_Ct02 = [20.000000000000021, 20.000000000000021, 20.000000000000014, 19.999999999999979, 20.000000000000004, 19.999999999999996, 20.000000000000011, 20.000000000000004, 20.000000000000004, 20.000000000000014, 19.999999999999986]
#Difference in V at min SAR due to increase
V_minSAR_Ct02 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#Overall difference in SAR due to increase
SAR_Ct02 = [20.000000000000004, 20.000000000000004, 20.000000000000004, 20.0, 20.000000000000004, 20.000000000000004, 20.000000000000004, 20.000000000000004, 20.0, 20.0, 20.000000000000004]

##plt.subplot(212)
#plt.plot(H,minSAR_Ct01,"r",label = "%change in min SAR ")
#plt.plot(H,V_minSAR_Ct01,"b",label = "%change in V at min SAR")
#plt.plot(H,SAR_Ct01,"g",label = "%average change in SAR")
#plt.title("Changes due to a 20% decrease in Ct0")
#plt.xlabel("Altitude [m]")
#plt.ylabel("Change in %")
#plt.legend(loc = "upper left")
#plt.show()
#
##plt.subplot(211)
#plt.plot(H,minSAR_Ct02,"r",label = "%change in min SAR ")
#plt.plot(H,V_minSAR_Ct02,"b",label = "%change in V at min SAR")
#plt.plot(H,SAR_Ct02,"g",label = "%average change in SAR")
#plt.title("Changes due to a 20% increase in Ct0")
#plt.xlabel("Altitude [m]")
#plt.ylabel("Change in %")
#plt.legend()
#plt.show()









