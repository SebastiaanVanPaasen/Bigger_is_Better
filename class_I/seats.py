import numpy as np


def Seats(seats_ab):
    Window = 0.
    Aisle = 0.

    seatcounter = seats_ab * 1
    if seatcounter[0][0] != 0:
        Window = Window + 1
        seatcounter[0][0] = seatcounter[0][0] - 1
    if seatcounter[0][2] != 0:
        Window = Window + 1
        seatcounter[0][2] = seatcounter[0][2] - 1
    if seatcounter[1][0] != 0:
        Window = Window + 1
        seatcounter[1][0] = seatcounter[1][0] - 1
    if seatcounter[1][2] != 0:
        Window = Window + 1
        seatcounter[1][2] = seatcounter[1][2] - 1

    if seatcounter[0][0] != 0:
        Aisle = Aisle + 1
        seatcounter[0][0] = seatcounter[0][0] - 1
    if seatcounter[0][2] != 0:
        Aisle = Aisle + 1
        seatcounter[0][2] = seatcounter[0][2] - 1
    if seatcounter[1][0] != 0:
        Aisle = Aisle + 1
        seatcounter[1][0] = seatcounter[1][0] - 1
    if seatcounter[1][2] != 0:
        Aisle = Aisle + 1
        seatcounter[1][2] = seatcounter[1][2] - 1
    if seatcounter[0][1] != 0:
        Aisle = Aisle + 2
        seatcounter[0][1] = seatcounter[0][1] - 2
    if seatcounter[1][1] != 0:
        Aisle = Aisle + 2
        seatcounter[1][1] = seatcounter[1][1] - 2

    Middle = sum(seatcounter[0]) + sum(seatcounter[1])

    return (Window, Aisle, Middle)


def cg_seats(Pseat, N_rows_above, N_rows_below, wing_start, start_pax, space_below_2):
    n_rows1 = int((wing_start - start_pax)/Pseat)
    loc_front = n_rows1
    xcg_seats = [start_pax]
    for i in range(max([N_rows_below, N_rows_above]) - 2):
        
        xcg_seats.append(xcg_seats[i] + Pseat)
        if i == n_rows1:
           xcg_seats[-1] = xcg_seats[-1] + space_below_2
    
#    print(xcg_seats)
#    print(space_below_2)
#    print(n_rows1)
    return (xcg_seats, loc_front)


def W_seats(tot_seating_abreast, N_rows_above, N_rows_below, tot_seating_abreast_last_row, loc_front):
    # print(tot_seating_abreast)
    W_window = min(N_rows_below, N_rows_above) * [Seats(tot_seating_abreast)[0]]
    W_aisle = min(N_rows_below, N_rows_above) * [Seats(tot_seating_abreast)[1]]
    W_middle = min(N_rows_below, N_rows_above) * [Seats(tot_seating_abreast)[2]]
    tot_seating_abreast_last_row = np.array([[2.,0.,2.], [2.,2.,2.]])
    remain_rows = abs(N_rows_below - N_rows_above)
#    print(tot_seating_abreast)
    if remain_rows > 0:
        
        if N_rows_below > N_rows_above:
            tot_seating_abreast[0] = np.array([0., 0., 0.])
        if N_rows_below < N_rows_above:
            tot_seating_abreast[1] = np.array([0., 0., 0.])
#        print(tot_seating_abreast, tot_seating_abreast_last_row)
        if (sum(tot_seating_abreast_last_row[0]) + sum(tot_seating_abreast_last_row[1])) != 0:
#            W_window = W_window + ((remain_rows - 1) * [Seats(tot_seating_abreast)[0]])
#            W_window = W_window + ([Seats(tot_seating_abreast_last_row)[0]])
#            W_aisle = W_aisle + ((remain_rows - 1) * [Seats(tot_seating_abreast)[1]]) 
#            W_aisle = W_aisle + ([Seats(tot_seating_abreast_last_row)[1]])
#            W_middle = W_middle + ((remain_rows - 1) * [Seats(tot_seating_abreast)[2]])
#            W_middle = W_middle + ([Seats(tot_seating_abreast_last_row)[2]])
#            print(tot_seating_abreast_last_row)
             tot_seating_abreast_last_row = np.array([[2.,0.,2.], [2.,2.,2.]])
#            W_window = W_window[0:loc_front] + ((remain_rows-1) * [Seats(tot_seating_abreast)[0]]) + ([Seats(tot_seating_abreast_last_row)[0]]) + W_window[loc_front:]
#            W_aisle = W_aisle[0:loc_front] + ((remain_rows-1) * [Seats(tot_seating_abreast)[1]]) + ([Seats(tot_seating_abreast_last_row)[1]]) + W_aisle[loc_front:]
#            W_middle = W_middle[0:loc_front] + ((remain_rows-1) * [Seats(tot_seating_abreast)[2]]) + ([Seats(tot_seating_abreast_last_row)[2]]) + W_middle[loc_front:]
             W_window = W_window[0:loc_front] + ((remain_rows) * [Seats(tot_seating_abreast)[0]]) + W_window[loc_front:-2] + 2*([Seats(tot_seating_abreast_last_row)[0]]) 
             W_aisle = W_aisle[0:loc_front] + ((remain_rows) * [Seats(tot_seating_abreast)[1]])  + W_aisle[loc_front:-2] + 2*([Seats(tot_seating_abreast_last_row)[1]])
             W_middle = W_middle[0:loc_front] + ((remain_rows) * [Seats(tot_seating_abreast)[2]]) + W_middle[loc_front:-2] + 2*([Seats(tot_seating_abreast_last_row)[2]])
            
        else:
#            W_window = W_window + ((remain_rows) * [Seats(tot_seating_abreast)[0]])
#            W_aisle = W_aisle + ((remain_rows) * [Seats(tot_seating_abreast)[1]])
#            W_middle = W_middle + ((remain_rows) * [Seats(tot_seating_abreast)[2]])
            
            W_window = W_window[0:loc_front] + ((remain_rows) * [Seats(tot_seating_abreast)[0]]) + W_window[loc_front:]
            W_aisle = W_aisle[0:loc_front] + ((remain_rows) * [Seats(tot_seating_abreast)[1]]) + W_aisle[loc_front:]
            W_middle = W_middle[0:loc_front] + ((remain_rows) * [Seats(tot_seating_abreast)[2]]) + W_middle[loc_front:]
    else:
        W_window[-1] = Seats(tot_seating_abreast_last_row)[0]
        W_aisle[-1] = Seats(tot_seating_abreast_last_row)[1]
        W_middle[-1] = Seats(tot_seating_abreast_last_row)[2]
        
#    print(W_window, W_aisle, W_middle)
    print(sum(W_window) + sum(W_aisle) + sum(W_middle))
    return (W_window, W_aisle, W_middle)
