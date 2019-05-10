# -*- coding: utf-8 -*-
"""
Created on Thu May  9 15:08:57 2019

@author: Hidde
"""
import numpy as np
from fuselage_cross_section import *

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
Npax = 450
Npax_below = 300
d_inner, d_outer, lcabin, lcabin_below, lcabin_above, l_tailcone_range, l_nosecone_range, l_fuselage,tot_seating_abreast, N_aisle_above, N_aisle_below,N_rows_above,N_rows_below,lpax_below,lpax_above,tot_seating_abreast_last_row,Pseat=fuselage_cross_section(Npax,Npax_below)


def cg_seats(tot_seating_abreast,N_rows_above,N_rows_below,lpax_below,lpax_above,l_nosecone_range,tot_seating_abreast_last_row):
    l_nosecone = (l_nosecone_range[0]+l_nosecone_range[1])/2.
    location_last_row_below = l_nosecone+lpax_below
    location_last_row_above = l_nosecone+lpax_above
    xcg_seats = [l_nosecone]
    for i in range(max([N_rows_below,N_rows_above])-1):
        xcg_seats.append(xcg_seats[i]+Pseat)
    return(xcg_seats)


def W_seats(tot_seating_abreast,N_rows_above,N_rows_below,lpax_below,lpax_above,l_nosecone_range,tot_seating_abreast_last_row):
    print(tot_seating_abreast)
    W_window = min(N_rows_below, N_rows_above)*[Seats(tot_seating_abreast)[0]]
    W_aisle = min(N_rows_below, N_rows_above)*[Seats(tot_seating_abreast)[1]]
    W_middle = min(N_rows_below, N_rows_above)*[Seats(tot_seating_abreast)[2]]
    print(W_window, W_aisle, W_middle)
    return

W_seats(tot_seating_abreast,N_rows_above,N_rows_below,lpax_below,lpax_above,l_nosecone_range,tot_seating_abreast_last_row)
    