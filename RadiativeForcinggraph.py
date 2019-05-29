import numpy as np
#import matplotlib.pyplot as plt

def FuelVRFDOC(h,M,Fnew):

    #Machspeed(h,M):
    #Radiative(h,Fincini):
    #Costs(Vini, Fincini)
    F737 = 0.0094

    F_dif= ((Fnew - F737) / F737 ) 
    
    Vini = Machspeed(h,M)
    
    RFreduction = Radiative(h,F_dif)
    DOCreduction = Costs(Vini, F_dif)
    
    F_dif = F_dif * 100
    #print Fincini
    #print Vini
    #print RFreduction
    #print DOCreduction
    print(DOCreduction)
    return F_dif, RFreduction, DOCreduction


def Machspeed(h,M):
    def ISA_temp(h):
        if h < 11000:
            T = 288.15 - 0.0065 * h  # in Kelvin
            return T
        if h >= 11000:
            return 216.65  # in Kelvin
        

    def Velocity(M, h):
        gamma = 1.4
        R = 287
        a = np.sqrt(gamma * R * ISA_temp(h))
        V = M * a
        return V
    return Velocity(M,h)


def Costs(Vini, Fincini):
     V737 = 828.479
     Vini = Vini * 3.6
     
     Vincs = [ 0.9, 0.95, 1, 1.05, 1.1, ] #
     Fincs = [ 1.,1., 1., 1., 1., ]
     TOCres = [0] * len(Vincs)
     Doc =[0] * len(Vincs)
     Doct =[0] * len(Vincs)
     TOCcor= [0] * len(Vincs)
     for i in range(len(Vincs)):
     #Changeable variables
      Vnew = Vincs[i] * Vini # new speed
      Finc = Fincs[i] # fuel increase
      R = 1400. # range
      Takelanding = 25.
      Turnaround = 80.

     #Direct operative costs
     #constants, from literature
      DOCnames =  ['Fuel', 'Crew', 'Depreciation', 'Maintenance']

      DOCs737 = [40., 20., 20., 20.]
      DOC737 = sum(DOCs737) / 100

      DOCsini = [40., 13., 20., 20.]
      DOCini = sum(DOCsini) / 100
      DOCs = [_ for _ in DOCsini]
      DOCs[0] = DOCs[0]
      DOC = sum(DOCs) / 100 / DOCini# new relative DOC due to increase in fuel

     #Total operative costs
      TOCnames =  ['DOCT', 'IOC', 'OHC']

      TOCs737 = [55., 20., 25.] 
      TOC737 = sum(DOCs737) / 100
      #TOCsini = [55., 20., 25.] 
      #TOCini = sum(DOCsini) / 100

    #Total costs increase due to increase in just fuel consumption
    #DOCinc = DOC/ DOCini
    #TOCs = [_ for _ in TOCsini]
    #TOCs[0] = TOCs[0] * DOCinc
    #TOC = sum(TOCs) /100

    #Flight trajectory
      Vini = Vini #km/h speed

      FTnames = ['Cruise','Takelanding' , 'Turnaround'] # times of procedures
      FT737 = [ R/V737, Takelanding /60 ,  Turnaround /60]
      Fltime7 = FT737[0] + FT737[1]
      TT737 = sum(FT737) #total time
      TOCh737 = [0] * 3

      
      TOCh737[0] = TOCs737[0] / Fltime7 # costs per hour per flight time
      TOCh737[1] = TOCs737[1] / FT737[2] # per 
      TOCh737[2] = TOCs737[2] / TT737

      FTini = [ R/Vini, Takelanding /60 ,  Turnaround /60]
      Fltime = FTini[0] + FTini[1]
      TTini = sum(FTini) #total time
      TOCsini = [0] * 3
      TOChini = [_ for _ in TOCh737]
      TOChini[0] = TOChini[0] * DOCini # costs per hour
     
     
      TOCsini[0] = TOChini[0] * Fltime # costs per hour per thing
      TOCsini[0] = TOCsini[0] - (DOCsini[0]/ DOCini  * (1 - Fincini)) # reduction of fuel consumption per km
      TOCsini[1] = TOChini[1] * FTini[2]
      TOCsini[2] = TOChini[2] * TTini
      TOCini = sum(TOCsini)
     
     
      TOCh = [_ for _ in TOChini]
      TOCh[0] = TOCh[0] * DOC # costs per hour
      FT = [_ for _ in FTini]
      FT[0] = R/Vnew
      TT = sum(FT)

      Flt = FT[0] + FT[1]

     
      DOCs[0] = DOCs[0] * Finc * (Fincini + 1)
      #print DOCs[0]

      #for i in range(1,4):
        #  print DOCs[i] * Flt / Fltime7
    
      TOCs = [0] * 3
      TOCs[0] = TOCh[0] * Flt
      TOCs[0] = TOCs[0] - (DOCs[0]/ DOCini  * (1-(Finc*(1 + Fincini))))
      TOCs[1] = TOCh[1] * FT[2]
      TOCs[2] = TOCh[2] * TT
     
      TOC = sum(TOCs)

      TOCres[i] = TOC
      Corf = TT/ TT737
      TOCcor[i] = TOCres[i] * Corf
      Doc[i] = DOC
      Doct[i] = ((TOCs[0] - TOCs737[0]) / TOCs737[0]) * 100 

     Dif = (TOCini / TOC737) /100
     #print Dif
     NTOC = [x /1 for x in TOCres]
     NDoc  = [x / DOC737 * 0.9 *100 for x in Doc]
     Vnews  = [x * Vini for x in Vincs]
     D737  = [x * Dif for x in NTOC]
     F100  = [x *100. for x in Fincs]
     
     #gr1 = plt.plot( TOCres  , Vnews,label='TOC compared to 737' )
     #gr5 = plt.plot( TOCcor  , Vnews,label='TOC with correction' )
     #gr2 = plt.plot(TOCcor, Vnews, )
     #gr3 = plt.plot(NDoc , Vnews, label='DOC per hour' )
     #gr3 = plt.plot(Doct , Vnews, label='DOC cost per flight' )
     #gr4 = plt.plot(D737, Vnews, label='TOC compared to initial 737')
     #gr4 = plt.plot(F100, Vnews, label='Fuelincrease')
     #plt.axvline(x=100)
     #plt.axhline(y=Vini)
     #plt.axvline(x=90)
     #plt.legend()
     #print TOCres[2]
     #plt.legend([gr1, gr2], ['TOC differnce', 'With time correction'])
     #plt.axis([0, 6, 0, 20])
     #plt.show()
     
     pt90 = 745
     pt88 = 733
     pt85 = 715
     
     return Doct[2]
 
    
def Radiative(h,Fincini):
    H737 = 12000
    Hini = 11000
    h = h
     
    
    trfi = [0] * 1
    for i in range(1):


     TotalRFI = 100.
     CO2RFI = 52.
     NOXRFI = 23.4
     CONTRAILRFI = 21.9
     H2ORFI = 5.2
     OtherRFI = -2.6

     #Contrail reduction, starting at 11000km
     Dif11km = Hini- h
     #Effect of difference
     
     deltacontrail = 0.000246
     if Dif11km <= 3000:
      conrfi  =  (1-( deltacontrail * Dif11km)) * CONTRAILRFI
      h2rfi =( 1 - (deltacontrail * Dif11km)) * H2ORFI
      NOXrfi = 23.4
     elif Dif11km <= 4064:
      conrfi  =  (1-( deltacontrail * Dif11km)) * CONTRAILRFI
      h2rfi =( 1 - (deltacontrail * Dif11km)) * H2ORFI
      NOXrfi = 23.4 * 0.9
     elif Dif11km <= 5000:
      conrfi  =  0
      h2rfi = 0
      NOXrfi = 23.4 * 0.5
     elif Dif11km <= 6000:
      conrfi  =  0
      h2rfi = 0
      NOXrfi = 23.4 * 0.3
     else :
      conrfi  =  0
      h2rfi = 0
      NOXrfi = 23.4 * 0.2

    
     totrfi = CO2RFI + NOXrfi + conrfi + h2rfi + OtherRFI
     trfi = totrfi * (Fincini + 1)


     Dif11km = Hini- H737
     if Dif11km <= 3000:
      conrfi7  =  (1-( deltacontrail * Dif11km)) * CONTRAILRFI
      h2rfi7 =( 1 - (deltacontrail * Dif11km)) * H2ORFI
      NOXrfi7 = 23.4
    totrfi7 = CO2RFI + NOXrfi7 + conrfi7 + h2rfi7 + OtherRFI
    
    #plt.plot( h, trfiinc)
    #plt.plot( h, Finc) 
    #plt.axis([0, 6, 0, 20])
    #plt.show()
    Dif737 = ((trfi - totrfi7) / totrfi7 ) * 100
    #print TotalRFI * Fincini
    #print CO2RFI * Fincini
    #print NOXrfi * Fincini
    #print conrfi / conrfi7 * CONTRAILRFI * Fincini
    #print h2rfi / h2rfi7 * H2ORFI * Fincini
#    print(OtherRFI * Fincini)
    #print Dif737 
    #print Fincini
    
    return Dif737
#
#hkm= [ 12.,11.,10.,9.,8.,7.,6.,5.,]
#h = [ x * 1000. for x in hkm]
#M10= [7.,7.,7.,7.,7.,7.,7.,7.]
#M = [ x / 10 for x in M10]
#F1000 = [80.,76.,75.,76.,78.,83.,88.,97.]
#Fnew = [ x / 10000 for x in F1000]
#values = [0] * len(h)
#RF = [0] * len(h)
#F = [0] * len(h)
#both = [0] * len(h)
#
#for i in range(len(h)):
#    values[i] = FuelVRFDOC(h[i],M[i],Fnew[i])
#    RF[i] = FuelVRFDOC(h[i],M[i],Fnew[i])[2] / 0.7 * 100 -100
#    F[i] = FuelVRFDOC(h[i],M[i],Fnew[i])[0] / 0.8 * 100 -100
#print(values)
#
#plt.xlabel('Altitude [km]')
#plt.ylabel('Difference [%]')
#show = plt.plot(hkm,RF, label='Radiative Forcing')
#show2 = plt.plot(hkm,F, label='Fuel Consumption')
#plt.grid('on')
#plt.legend()
#plt.show()
#
#result = [0] * len(h)
#for i in range(len(h)):
#    result[i] = values[i][0] * 2 + values[i][2] * 5 + values[i][3] * 5 /100

#print result
