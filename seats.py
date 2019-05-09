# -*- coding: utf-8 -*-
"""
Created on Thu May  9 15:08:57 2019

@author: Hidde
"""
import numpy as np

seatcounter = np.array([[3, 4, 3], [3, 0, 3]])

def Seats(seatcounter):  
    Window = 0
    Aisle  = 0
    Middle = 0
    
    if seatcounter[0][0] != 0:
        Window = Window + 1
        seatcounter[0][0] = seatcounter[0][0]-1   
    if seatcounter[0][2] != 0:
        Window = Window + 1
        seatcounter[0][2] = seatcounter[0][2]-1
    if seatcounter[1][0] != 0:
        Window = Window + 1
        seatcounter[1][0] = seatcounter[1][0]-1
    if seatcounter[1][2] != 0:
        Window = Window + 1
        seatcounter[1][2] = seatcounter[1][2]-1
        
    if seatcounter[0][0] != 0:
        Aisle = Aisle + 1
        seatcounter[0][0] = seatcounter[0][0]-1
    if seatcounter[0][2] != 0:
        Aisle = Aisle + 1
        seatcounter[0][2] = seatcounter[0][2]-1
    if seatcounter[1][0] != 0:
        Aisle = Aisle + 1
        seatcounter[1][0] = seatcounter[1][0]-1
    if seatcounter[1][2] != 0:
        Aisle = Aisle + 1
        seatcounter[1][2] = seatcounter[1][2]-1
    if seatcounter[0][1] != 0:
        Aisle = Aisle + 2
        seatcounter[0][1] = seatcounter[0][1]-2
    if seatcounter[1][1] != 0:
        Aisle = Aisle + 2
        seatcounter[1][1] = seatcounter[1][1]-2
    
    Middle = sum(seatcounter[0])+sum(seatcounter[1])
        
    return(Window, Aisle, Middle)
print(Seats(seatcounter))