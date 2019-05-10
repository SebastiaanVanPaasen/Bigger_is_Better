# -*- coding: utf-8 -*-
"""
Created on Thu May  9 15:08:57 2019

@author: Hidde
"""
import numpy as np
from fuselage_cross_section import *

def Seats(seats_ab):  
    Window = 0
    Aisle  = 0
    Middle = 0
    seatcounter = seats_ab*1
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
    Window = float(Window)
    Aisle = float(Aisle)
    Middle = sum(seatcounter[0])+sum(seatcounter[1])
        
    return(Window, Aisle, Middle)





def cg_seats(Pseat, tot_seating_abreast,N_rows_above,N_rows_below,lpax_below,lpax_above,l_nosecone_range,tot_seating_abreast_last_row):
    l_nosecone = (l_nosecone_range[0]+l_nosecone_range[1])/2.
    location_last_row_below = l_nosecone+lpax_below
    location_last_row_above = l_nosecone+lpax_above
    xcg_seats = [l_nosecone]
    for i in range(max([N_rows_below,N_rows_above])-1):
        xcg_seats.append(xcg_seats[i]+Pseat)
    return(xcg_seats)


def W_seats(tot_seating_abreast,N_rows_above,N_rows_below,lpax_below,lpax_above,l_nosecone_range,tot_seating_abreast_last_row):
    #print(tot_seating_abreast)
    W_window = min(N_rows_below, N_rows_above)*[Seats(tot_seating_abreast)[0]]
    W_aisle = min(N_rows_below, N_rows_above)*[Seats(tot_seating_abreast)[1]]
    W_middle = min(N_rows_below, N_rows_above)*[Seats(tot_seating_abreast)[2]]
    
    remain_rows=abs(N_rows_below-N_rows_above)
    
    if remain_rows > 0:
        if N_rows_below > N_rows_above:
            tot_seating_abreast[0] = np.array([0.,0.,0.])
        if N_rows_below < N_rows_above:
            tot_seating_abreast[1] = np.array([0.,0.,0.])
        
        if (sum(tot_seating_abreast_last_row[0])+sum(tot_seating_abreast_last_row[1])) > 0:
            W_window = W_window+((remain_rows-1)*[Seats(tot_seating_abreast)[0]])
            W_window = W_window+([Seats(tot_seating_abreast_last_row)[0]])
            W_aisle = W_aisle+((remain_rows-1)*[Seats(tot_seating_abreast)[1]])
            W_aisle = W_aisle+([Seats(tot_seating_abreast_last_row)[1]])
            W_middle = W_middle+((remain_rows-1)*[Seats(tot_seating_abreast)[2]])
            W_middle = W_middle+([Seats(tot_seating_abreast_last_row)[2]])
        else:
            W_window = W_window+((remain_rows)*[Seats(tot_seating_abreast)[0]])
            W_aisle = W_aisle+((remain_rows)*[Seats(tot_seating_abreast)[1]])
            W_middle = W_middle+((remain_rows)*[Seats(tot_seating_abreast)[2]])
    
    else: 
        W_window[-1] = Seats(tot_seating_abreast_last_row)[0]
        W_aisle[-1] = Seats(tot_seating_abreast_last_row)[1]
        W_middle[-1] = Seats(tot_seating_abreast_last_row)[2]
    

    return(W_window, W_aisle, W_middle)

