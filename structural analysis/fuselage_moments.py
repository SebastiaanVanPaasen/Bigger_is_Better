import numpy as np
from math import *
## define forces
fuse_W = -268536
Vertical_tail = -10075
fuse_L = 
wing_component = 150000
h_tail_component = 298056
v_tail_component = 0.0
strut = 

## define positions
fuse_length = 49.3
wing_pos = 20
h_tail_pos = 51
v_tail_pos = 48
cg_pos = 
strut_pos = 
main_lg_pos = 
nose_lg_pos =
empennage_pos = 
fuse_L_pos = 


##angles --> define the angle of attack t which the forces are in equilibrium during flight!!
aoa = np.arange(0.,40.,0.0001)
for i in aoa:
    #print(i)
    fuse_L = sin(90-i)*fuse_L
    wing_component = sin(90-i)*wing_component
    h_tail_component = sin(90-i)*h_tail_component
    v_tail_component = sin(90-i)*v_tail_component
    strut = sin(90-i)*strut
    neg_forces = fuse_W + empennage_W
    pos_forces = fuse_L+wing_component+h_tail_component+v_tail_component+strut
    #print(neg_forces)
    diff = neg_forces-pos_forces
    if abs(diff)<0.1:
        print(i)
        break


##Moments
#Moment around the cg, right hand rule positive:
dx = 0.01
length = []
moment = []
for i in range(0,fuse_L,dx):
    M = fuse_w*(cg_pos-i)+empennage_W*(empennage_pos-i)+fuse_L*(fuse_L_pos-i)+wing_component*(wing_pos-i)+h_tail_component*(h_tail_pos-i)+v_tail_component*(v_tail_pos-i)+strut*(strut_pos-i)
    length.append(i)
    moment.append(M)
##during ground operations

