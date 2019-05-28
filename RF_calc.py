
import numpy as np

def ISA_temp(h):
    if h<11000:
        Temp = 288.15 - 0.0065*h
    else: 
        Temp = 216.65
    return Temp
        
def Velocity(M,Temp):
    gamma = 1.4
    speed_of_sound = np.sqrt(gamma*287*Temp)
    V = M*speed_of_sound
    return(V)


def Radiative(h, Fnew):  #Fnew is the MTOW fuel fraction of our design
    H737 = 12000 #Reference height of B737
    Hini = 11000 #Reference height for contrail shit
    
#    TotalRFI = 100.  #Impact of the various RF's
    CO2RFI = 52.
    NOXRFI = 23.4
    CONTRAILRFI = 21.9
    H2ORFI = 5.2
    OtherRFI = -2.6
    NOXrfi =23.4
    deltacontrail = 0.000246
    Dif11km = Hini- h
    
    # RF's contributions setup
    conrfi = 0
    h2rfi = 0
    if Dif11km <= 4064:
        conrfi = (1-( deltacontrail * Dif11km)) * CONTRAILRFI
        h2rfi = ( 1 - (deltacontrail * Dif11km)) * H2ORFI
    
    
    if Dif11km <= 3000 :
        NOXrfi = NOXRFI
    elif Dif11km <=4064:
        NOXrfi = NOXrfi * 0.9
    elif Dif11km <= 5000:
        NOXrfi = NOXrfi*0.5
    elif Dif11km <= 6000:
        NOXrfi = NOXrfi*0.3
    elif Dif11km >6000:
        NOXrfi = NOXrfi *0.2
        
    F737 = 0.0094
    F_frac = Fnew/F737                                              #Below is the RF of our design
    totrfi = (CO2RFI + NOXrfi + conrfi + h2rfi + OtherRFI)*F_frac #Total de/increase in RFi due to delta in fuel weight 
    
    conrfi7  =  (1-( deltacontrail * (Hini-H737))) * CONTRAILRFI
    h2rfi7 =( 1 - (deltacontrail * (Hini-H737))) * H2ORFI
    NOXrfi7 = 23.4    
    totrfi7 = CO2RFI + NOXrfi7 + conrfi7 + h2rfi7 + OtherRFI #737 Radiative forces
    
    Dif737 = (1-(totrfi/totrfi7))*100
    return Dif737, F_frac
#print(Radiative(13000,0.0094))
    
Vini = Velocity(.78,ISA_temp(12000))
#print("Velocity " + str(Vini))

def Costs(Vini, Fnew):
    F737 = 0.0094
    F_frac = Fnew/F737
    
    V737 = 828.479 #km/h
    Vini = Vini * 3.6 #m/s->km/h
#    DOC = .55
    DOC_frac_crew = 0.93
    Range = 1400 #km
    TO_landing = 25 #minutes
    flight_time_737 = Range/V737 +TO_landing/60
    flight_time = Range/Vini + TO_landing/60
    flight_time_fraction = flight_time/flight_time_737 
#    print(flight_time_fraction)
    DOC_reduction_flighttime = (1-flight_time_fraction)
        
#    cost_per_hour737 =DOC/flight_time_737
#    cost_per_hour = DOC/flight_time
    
    DOC_reduction_fuel = 0.4*(1-F_frac)
    DOC_reduction_crew = (1-DOC_frac_crew)
    total_DOC_reduction = (DOC_reduction_flighttime + DOC_reduction_fuel + DOC_reduction_crew)*100
    return(total_DOC_reduction)
#print(Costs(Vini,0.0094))
    
    
    
