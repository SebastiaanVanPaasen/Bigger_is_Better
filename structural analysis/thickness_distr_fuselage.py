import numpy as np
from matplotlib import pyplot as plt


def t_fus_due_to_moment(sigma, Mx, My):
    
    
#    file = open("C:/Users/sebas/OneDrive/Documents/DSE/Bigger_is_Better/structural analysis/fus_coor.txt", 'r')
    file = open("C:/Users/mathi/Documents/DSE/Bigger_is_Better/structural analysis/fus_coor.txt", 'r')

    x_pos = []
    y_pos = []
    
    for line in file:
        x_pos.append(float(line.split()[0]))
        y_pos.append(float(line.split()[1]))
        
    x_pos = np.array(x_pos[::-1])
    y_pos = -np.array(y_pos[::-1])


    R_2 = ((x_pos**2 + y_pos**2)**0.5) / 1000
    
#    print(R_2)
#    print("ypos", y_pos)
    alpha_segment = []
    
    #Iteration 1 of the variable thickness
    for i in range(len(y_pos)):
        
        if y_pos[i]>0 and x_pos[i]>0:
            alpha_segment.append(np.arctan(y_pos[i]/x_pos[i]))
            
        elif y_pos[i]>0 and x_pos[i]<0:
            alpha_segment.append(np.arctan(-x_pos[i]/y_pos[i]) + np.pi/2)
            
        elif y_pos[i]<0 and x_pos[i]<0:
            alpha_segment.append(np.arctan(-y_pos[i]/-x_pos[i]) + np.pi)
            
        elif y_pos[i]<0 and x_pos[i]>0:
            alpha_segment.append(np.arctan(x_pos[i]/-y_pos[i]) + np.pi*3/2)#*180/np.pi#*np.pi/180
    
    alpha_segment = np.array(alpha_segment)
  
    print(alpha_segment)
    
    theta = list(np.arange(30, 35, 10) * np.pi / 180) #list(np.array([-30, 0, 30])* np.pi / 180) # #angle between the moment and the z-axis
    t_circ = np.zeros((len(theta), len(alpha_segment)))
    t_max_all = []
    
    for i in range(len(theta)):
        
        for j in range(len(alpha_segment)):
            
            part_mx = (Mx * np.cos(theta[i]) + My * np.sin(theta[i])) * np.sin(alpha_segment[j])
            part_my = (My * np.cos(theta[i]) - Mx * np.sin(theta[i])) * np.cos(alpha_segment[j])
            
            t = (part_mx + part_my) / (np.pi * (R_2[j] ** 2) * sigma)
            
            t_circ[i][j] = abs(t)
     
    
    for i in range(len(alpha_segment)):
        segment_values = []
        
        for j in range(len(theta)):
            
            segment_values.append(t_circ[j][i])
        
        t_max = max(segment_values)
        t_max_all.append(t_max)
   
    
    segment_length = np.zeros(len(alpha_segment))
    
    for i in range(len(R_2)-1):
        
        segment_length[i] = ((alpha_segment[i + 1] - alpha_segment[i]) / (2 * np.pi)) * (2 * np.pi * (R_2[i+1] + R_2[i]) / 2)

      
        
    segment_length[-1] = (alpha_segment[0] + 2 * np.pi - alpha_segment[-1]) * (R_2[0] + R_2[-1]) / 2
    
#    plt.rcParams.update({'font.size': 20})        
#    plt.figure(1)
#    plt.plot(alpha_segment*180/np.pi,t_circ[0]*1000,label = 'theta = -30')
#    plt.plot(alpha_segment*180/np.pi,t_circ[1]*1000, label = 'theta = -20')
#    plt.plot(alpha_segment*180/np.pi,t_circ[2]*1000, label = 'theta = -10')
#    plt.plot(alpha_segment*180/np.pi,t_circ[3]*1000, label = 'theta = 0')
#    plt.plot(alpha_segment*180/np.pi,t_circ[4]*1000, label = 'theta = 10')
#    plt.plot(alpha_segment*180/np.pi,t_circ[5]*1000, label = 'theta = 20')
#    plt.plot(alpha_segment*180/np.pi,t_circ[6]*1000, label = 'theta = 30')
#    plt.plot(alpha_segment*180/np.pi,np.array(t_max_all)*1000, label = 'maximum thickness')
#
#    plt.xlabel("\u03B1 [degrees]")
#    plt.ylabel("Thickness distribution [mm]")
#
#    plt.legend(bbox_to_anchor=(1.03,1), loc="upper left")    
#    plt.show()
    
#    print(t_max_all)
    
    return(t_max_all, alpha_segment, R_2, segment_length, theta, t_circ)
        

    
def max_stress(theta, alpha_segment, I_xx, I_yy, Mx, My, sigma):
    
    stress_ratios = np.zeros((len(alpha_segment), len(theta)))
    max_stress_ratios = []
    
    for i in range(len(theta)):
        
        for j in range(len(alpha_segment)):  
            
            part_mx = (Mx * np.cos(theta[i]) + My * np.sin(theta[i])) * R_2[j] * np.sin(alpha_segment[j]) / I_xx
            part_my = (My * np.cos(theta[i]) - Mx * np.sin(theta[i])) * R_2[j] * np.cos(alpha_segment[j]) / I_yy
            
            
#            print(part_mx)
#            print()
#            
#            print(part_my)
#            print()
            
            stress_ratio = (part_mx + part_my) / sigma
            stress_ratios[j][i] = stress_ratio
    
    
    for i in range(len(stress_ratios)):
        
#        print(stress_ratios[i])
        
        max_stress_ratios.append(max(abs(stress_ratios[i])))
       
    max_stress_ratios = np.asarray(max_stress_ratios)
    
#    print(max_stress_ratios)

    
    return max_stress_ratios

finalthicknessreached = False
counter = 0

sigma = 105 * (10 ** 6)   #ultimate stress of aluminium 
Mx = -37799412
My = 216.17 * (10 ** 3) * 7

t_iter, alpha_segment, R_2, segment_length, theta, t_circ = t_fus_due_to_moment(sigma, Mx, My)


while finalthicknessreached == False:  # and counter < 100: 
#    print(t_iter)
    
#    I_xx = t_circ[2] * (segment_length ** 3) * (np.sin(alpha_segment + 0.5 * np.pi) ** 2) / 12 + t_circ[2] * segment_length * (R_2 ** 2) * np.sin(alpha_segment)**2
#    I_yy = t_circ[2] * (segment_length ** 3) * (np.cos(alpha_segment + 0.5 * np.pi) ** 2) / 12 + t_circ[2] * segment_length * (R_2 ** 2) * np.cos(alpha_segment)**2
#    
#    
#    I_xx_tot = sum(I_xx)
#    I_yy_tot = sum(I_yy) 

    
    
    I_xx = t_iter * (segment_length ** 3) * (np.sin(alpha_segment + 0.5 * np.pi) ** 2) / 12 + t_iter * segment_length * (R_2 ** 2) * np.sin(alpha_segment)**2
    I_yy = t_iter * (segment_length ** 3) * (np.cos(alpha_segment + 0.5 * np.pi) ** 2) / 12 + t_iter * segment_length * (R_2 ** 2) * np.cos(alpha_segment)**2
    
    
    I_xx_tot = sum(I_xx)
    I_yy_tot = sum(I_yy) 
#    
#    print()
#    print(I_xx_tot)
#    print(I_yy_tot)
#    print()
##    
    stress_ratio = max_stress(theta, alpha_segment, I_xx_tot, I_yy_tot, Mx, My, sigma)
    
    #factor = 1/(counter+1)
    
#    t_iter = t_new
    
    
    counter += 1
    
    
#    print(counter)
#    print(stress_ratio)
    
    if counter > 49:
        break
##    
#    
    for i in range(len(stress_ratio)):
#        print(i)
        
        if 0.9 < stress_ratio[i] <= 0.99 or t_iter[i] < 0.001:
#            finalthicknessreached = True
#            t_circ[2][i] = t_circ[2][i]
            t_iter[i] = t_iter[i]
#            print(t_iter[1])
        else:
#            finalthicknessreached = False
#            print()
#            print(t_iter[i])
#            t_circ[2][i] = t_circ[2][i] * stress_ratio[i]
            t_iter[i] = t_iter[i] * stress_ratio[i]  #*(1+0.2*factor)
#            print(stress_ratio[i])
#            print(t_iter[i])
#            print()
            

            
plt.figure(2)
plt.plot(alpha_segment * 180 / np.pi, np.array(t_iter) * 1000, label = 'Iteration = ' + str(counter) + ", theta = 30")
plt.xlabel("\u03B1 [degrees]")
plt.ylabel("Thickness distribution [mm]")
plt.legend(bbox_to_anchor=(1.03,1), loc="upper left")    
plt.show()




sigma_fatigue_hoop = 114 * 10**6 # look up from graph
sigma_fatigue_long = 105 * 10**6
internal_p = 78.2 * 10**3 
external_p = 30.1 * 10**3


def t_fus_due_to_pressure(sigma_fatigue_hoop, sigma_fatigue_long,R,internal_p, external_p):
    delta_p = internal_p - external_p
    t_min_hoop = []
    t_min_long = []
    for i in range(len(R)):
        t_min_hoop.append(delta_p*R[i]/sigma_fatigue_hoop)
        t_min_long.append(delta_p*R[i]/(2*sigma_fatigue_long))
    return t_min_hoop, t_min_long

t_min_hoop,t_min_long = t_fus_due_to_pressure(sigma_fatigue_hoop, sigma_fatigue_long, R_2, internal_p, external_p)
print("tmin hoop", max(t_min_hoop))
print("tmin long", max(t_min_long))


#max_stress_ratios = np.zeros((1,len(t_iter)))

#for i in range(len(alpha_segment)):
#    I_xx += t_iter[i]*segment_length[0][i]**3*np.sin(alpha_segment[i]+0.5*np.pi)**2/12 + t_iter[i]*segment_length[0][i]*R**2*np.sin(alpha_segment[i])**2
#    I_yy += t_iter[i]*segment_length[0][i]**3*np.cos(alpha_segment[i]+0.5*np.pi)**2/12 + t_iter[i]*segment_length[0][i]*R**2*np.cos(alpha_segment[i])**2
#
#    stress_ratios = np.zeros((len(alpha_segment), len(theta)))
#        
#for i in range(len(theta)): 
#    for j in range(len(alpha_segment)):    
#        stress_ratio = (M*np.cos(theta[i])*R*np.sin(alpha_segment[j])/I_xx + -M*np.sin(theta[i])*R*np.cos(alpha_segment[j])/I_yy)/sigma
#        stress_ratios[j][i] = stress_ratio
#    #    print(np.shape(stress_ratios))
#    #    print(stress_ratios)
##max_stress_ratios=[]
##    print(t_old)
#print(t_iter[0:5])
#
#for i in range(len(stress_ratios)):
#    max_stress_ratios.append(max(abs(stress_ratios[i])))
#
#    t_new = t_iter[i]*max(abs(stress_ratios[i]))    
#    t_iter[i] = t_new
#    
#print(max_stress_ratios[0:5])   
#print(t_iter[0:5])
#finalthicknessreached = False


#for i in range(len(max_stress_ratios)):
#    I_xx = 0
#    I_yy = 0
#    
#    while finalthicknessreached == False:
#        
#        for p in range(len(alpha_segment)):
#            I_xx += t_iter[p]*segment_length[0][p]**3*np.sin(alpha_segment[p]+0.5*np.pi)**2/12 + t_iter[p]*segment_length[0][p]*R**2*np.sin(alpha_segment[p])**2
#            I_yy += t_iter[p]*segment_length[0][p]**3*np.cos(alpha_segment[p]+0.5*np.pi)**2/12 + t_iter[p]*segment_length[0][p]*R**2*np.cos(alpha_segment[p])**2
#
#        if 0.92>abs(max_stress_ratios[i]) or abs(max_stress_ratios[i])>0.99:
##            stress_ratios = np.zeros((len(alpha_segment), len(theta)))
##            for k in range(len(theta)): 
#            for j in range(len(alpha_segment)):#len(alpha_segment)):    
#                    stress_ratio = (M*np.cos(theta[k])*R*np.sin(alpha_segment[j])/I_xx + -M*np.sin(theta[k])*R*np.cos(alpha_segment[j])/I_yy)/sigma
#                    stress_ratios[j][k] = stress_ratio
#    
#            max_stress_ratios[i]= max(abs(stress_ratios[i]))
#            t_new = t_iter[i]*max(abs(stress_ratios[i]))    
#            t_iter[i] = t_new
#            print(max_stress_ratios[i])
#            print("thickness",t_iter[i])
#        
#        elif abs(max_stress_ratios[i])>= 0.92 and abs(max_stress_ratios[i])<= 0.99:
#            finalthicknessreached = True