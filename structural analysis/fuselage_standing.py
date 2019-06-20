import numpy as np

x_MTOW_for = 23.83
x_MTOW_aft = 26.66
MTOW = 1555183
nose_gear = 0.08*MTOW
main_gear = 0.92*MTOW
x_nose_gear = 6.6
x_main_gear = 28.4

x_MTOW = (x_MTOW_for+x_MTOW_aft)/2
E = 70*10**9
I = 0.005
tot_deflection_fuse = []
x = np.arange(x_nose_gear,x_main_gear,0.01)
for i in x:
    deflection_fuse = (((-1*MTOW)*(x_main_gear-x_MTOW)*i)/(6*E*I*(x_main_gear-x_nose_gear)))*((x_main_gear-x_nose_gear)**2-((x_main_gear-x_MTOW)**2)-(i**2))
    tot_deflection_fuse.append(deflection_fuse)
    
plt.plot(x, tot_deflection_fuse)
plt.show()    