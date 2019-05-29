# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:00:27 2019

@author: nikki
"""

import matplotlib.pyplot as plt


#High wing twin engine double decker
hcr = [8000,8500,8500,8500]
sar = [0.008875167,0.008652922,	0.008546917,0.008471038]
Mcr = [0.75,0.75,0.75,	0.75]

sar_ref = 0.00927126508106
h_ref = 11500
Mref = 0.79


for i in range(len(sar)):
    #plt.subplot(231)
    plt.figure(1)
    plt.plot(Mcr[i],sar[i],'o',label = "%.0f m and %.2f Mach" %(hcr[i],Mcr[i]))

plt.plot(Mref,sar_ref,'o',label = "Ref. aircraft")
plt.hlines(sar_ref*0.9,Mcr[0]-0.05,Mcr[-1]+0.05,"gray","--")    
plt.legend(loc = "upper left")
plt.title("HIGH 2E DD")   
plt.xlabel("Mach number")
plt.ylabel("Fuel consumption [kg/km-pax]")


#High wing twin engine strutted double decker
hcr = [11000,11000,11000,11500,11500,11500,12000,12000,11000,11000,11000,11500,11500,11500,11000,11000,11000]
sar = [0.007326261,0.007257847,0.00718606,0.007336119,0.007268737,0.007216544,0.007356067,0.007273407,0.007359841,0.007297933,0.00720945,0.007410994,0.007304704,0.007226114,0.007440889,0.007351348,0.00728756]
Mcr = [0.75,0.75,0.75,0.775,0.775,0.775,0.8,0.8,0.775,0.775,0.775,0.8,0.8,0.8,0.8,0.8,0.8]

sar_ref = 0.00927126508106
h_ref = 11500
Mref = 0.79

for i in range(len(sar)):
    #plt.subplot(232)
    plt.figure(2)
    plt.plot(Mcr[i],sar[i],'o',label = "%.0f m and %.2f Mach" %(hcr[i],Mcr[i]))
    
plt.plot(Mref,sar_ref,'o',label = "Ref. aircraft")
plt.hlines(sar_ref*0.9,Mcr[0]-0.02,Mcr[-1]+0.02,"gray","--")  
plt.xlim(min(Mcr)-0.02,max(Mcr)+0.02) 
plt.legend(loc = "upper left") 
plt.title("HIGH STRUT 2E DD") 
plt.xlabel("Mach number")
plt.ylabel("Fuel consumption [kg/km-pax]")

#High wing twin engine strutted semi double decker
hcr = [9000,9000,9000,10000,10000,10000,9000,9000,9000]
sar = [0.008491673,0.008392962,0.008323609,0.008360158,0.008277767,0.00819461,0.008282552,0.00818993,0.008189164]
Mcr = [0.75,0.75,0.75,0.8,0.8,0.8,0.75,0.75,0.75]

sar_ref = 0.00927126508106
h_ref = 11500
Mref = 0.79

for i in range(len(sar)):
    #plt.subplot(232)
    plt.figure(3)
    plt.plot(Mcr[i],sar[i],'o',label = "%.0f m and %.2f Mach" %(hcr[i],Mcr[i]))
    
plt.plot(Mref,sar_ref,'o',label = "Ref. aircraft")
plt.hlines(sar_ref*0.9,min(Mcr)-0.02,max(Mcr)+0.02,"gray","--")  
plt.xlim(min(Mcr)-0.02,max(Mcr)+0.02)  
plt.legend(loc = "upper left") 
plt.title("HIGH STRUT 2E SEMI DD") 
plt.xlabel("Mach number")
plt.ylabel("Fuel consumption [kg/km-pax]")

#High wing twin engine semi double decker
hcr = [9000,9000,9000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000]
sar = [0.009058995,0.008958824,0.008862814,0.008446061,0.008281633,0.008135173,0.008652075,0.008529623,0.008421779,0.008592276,0.008472932,0.008317172,0.008881601,0.008725746,0.008586266]
Mcr = [0.75,0.75,0.75,0.75,0.75,0.75,0.775,0.775,0.775,0.775,0.775,0.775,0.8,0.8,0.8]

sar_ref = 0.00927126508106
h_ref = 11500
Mref = 0.79

for i in range(len(sar)):
    #plt.subplot(232)
    plt.figure(4)
    plt.plot(Mcr[i],sar[i],'o',label = "%.0f m and %.2f Mach" %(hcr[i],Mcr[i]))
    
plt.plot(Mref,sar_ref,'o',label = "Ref. aircraft")
plt.hlines(sar_ref*0.9,min(Mcr)-0.02,max(Mcr)+0.02,"gray","--")  
plt.xlim(min(Mcr)-0.02,max(Mcr)+0.02)  
plt.legend(loc = "upper left") 
plt.title("HIGH 2E SEMI DD") 
plt.xlabel("Mach number")
plt.ylabel("Fuel consumption [kg/km-pax]")

#Low wing 4 engines semi double decker
hcr = [10000,10000,10000,10000,10000,10000,10000,10000,10000,9000,9000,9000,9000,9000,9000,10000,10000,10000,10000,10000,10000]
sar = [0.008521515,0.008363761,0.008230474,0.008758282,0.008627132,0.008496162,0.009028616,0.008931609,0.008835845,0.009003854,0.008871894,0.008788272,0.009327135,0.009205896,0.009136345,0.008701435,0.008513077,0.008383098,0.008911064,0.008775528,0.008659411]
Mcr = [0.75,0.75,0.75,0.775,0.775,0.775,0.8,0.8,0.8,0.75,0.75,0.75,0.775,0.775,0.775,0.775,0.775,0.775,0.8,0.8,0.8]

sar_ref = 0.00927126508106
h_ref = 11500
Mref = 0.79

for i in range(len(sar)):
    #plt.subplot(232)
    plt.figure(5)
    plt.plot(Mcr[i],sar[i],'o',label = "%.0f m and %.2f Mach" %(hcr[i],Mcr[i]))
    
plt.plot(Mref,sar_ref,'o',label = "Ref. aircraft")
plt.hlines(sar_ref*0.9,min(Mcr)-0.02,max(Mcr)+0.02,"gray","--")  
plt.xlim(min(Mcr)-0.02,max(Mcr)+0.02)  
plt.legend(loc = "upper left") 
plt.title("LOW 2E SEMI DD") 
plt.xlabel("Mach number")
plt.ylabel("Fuel consumption [kg/km-pax]")

plt.show()
