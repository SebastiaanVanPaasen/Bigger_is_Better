from math import *
from create_circle import *
import numpy as np


def fuselage_cross_section(Npax, Npax_below):
    Haisle = 1.95  # Height of the aisle [m] />60
    Waisle = 15 * 2.54 * 0.01  # Width of the aisle  [m] />15
    Wseat = 20 * 2.54 * 0.01  # Width of the seat   [m] /16-18
    Pseat = 32 * 2.54 * 0.01  # Seat Pitch          [m] /30-32
    Warmrest = 2 * 2.54 * 0.01  # Width of the armrest[m]
    Sclearance = 2 * 0.01  #
    Nos_lat = 2  # Number of lateral overhead storage space
    Nos_ce = 1  # Number of central overhead storage space
    Aos_lat = 0.20  # Area of the lateral overhead storage space [m^2]
    Aos_ce = 0.24  # Area of the central overhead storage space [m^2]
    Kos = 0.74  # Ration of the amount of storage places in the cabin
    rho_luggage = 170  # The typical density of the luggage  [kg/m^3]
    rho_cargo = 160  # the typical density of the crago    [kg/m^3]
    Ratio_tail_d_outer = 1.6  # Ratio of the tail length over outer diameter of the aircraft (for a typical aircraft)
    m_luggage = 40 * 0.45359237 * Npax
    m_cargo = 100
    H_shoulder = 1.10
    H_headroom = 1.65
    H_cargo = 64 * 2.54 * 0.01
    B_cargo = 96 * 2.54 * 0.01

    tot_seating_abreast = np.zeros([2, 3])
    # single floor configuration
    if Npax_below == Npax:
        N_sa_below = round(0.45 * sqrt(Npax))
        N_sa_above = 0
        # Check number of aisle and create the seating abreast
        if N_sa_below <= 6:
            k_cabin = 1.08
            N_aisles = 1
            N_right = ceil(N_sa_below / 2)
            N_left = N_sa_below - N_right
            N_middle = 0
        elif N_sa_below >= 6 and N_sa_below <= 12:
            k_cabin = 1.17
            N_aisles = 2
            if N_sa_below < 9:
                N_left = 2
                N_right = 2
            else:
                N_left = 3
                N_right = 3
            N_middle = N_sa_below - N_left - N_right

        else:
            print("This aircraft has more than 2 aisle")
            N_right = 0
            N_left = 0
            N_middle = 0
            N_aisles = 0
            N_sa_below = 0
            k_cabin = 0

        # Place the seating configuration in a array
        tot_seating_abreast[1, 0] = N_right
        tot_seating_abreast[1, 1] = N_middle
        tot_seating_abreast[1, 2] = N_left

        # Compute length of compartments
        lcabin = (Npax / N_sa_below) * k_cabin
        W_cabin = N_sa_below * Wseat + (N_sa_below + N_aisles + 1) * Warmrest + N_aisles * Waisle + 2 * Sclearance
        W_headroom = W_cabin - 2 * (Warmrest + Sclearance) - Wseat
        V_os = (Nos_lat * Aos_lat + Nos_ce * Aos_ce) * lcabin * Kos

        # create points to plot the layout if neccessary, these points are used to create a circle
        W_cabin_left_x = -W_cabin / 2
        W_headroom_left_x = -W_headroom / 2
        W_cabin_right_x = W_cabin / 2
        W_headroom_right_x = W_headroom / 2
        Waisle_left_x = -Waisle / 2
        Waisle_right_x = Waisle / 2

        cabin_left_point = (W_cabin_left_x, 0)
        cabin_right_point = (W_cabin_right_x, 0)
        shoulder_left_point = (W_cabin_left_x, H_shoulder)
        shoulder_right_point = (W_cabin_right_x, H_shoulder)
        headroom_left_point = (W_headroom_left_x, H_headroom)
        headroom_right_point = (W_headroom_right_x, H_headroom)
        aisle_left_point = (Waisle_left_x, Haisle)
        aisle_right_point = (Waisle_right_x, Haisle)

        # Function that compute the most inner circle of the points
        r = make_circle(
            [cabin_left_point, cabin_right_point, headroom_left_point, headroom_right_point, aisle_left_point,
             aisle_right_point])
        d_inner = 2 * r[2]
        d_outer = 1.045 * d_inner + 0.084

        # plot the cross_section
        #        plot_array_x = [W_cabin_left_x,W_cabin_right_x,W_cabin_left_x,W_cabin_right_x,W_headroom_left_x,W_headroom_right_x,Waisle_left_x,Waisle_right_x]
        #        plot_array_y = [0,0,H_shoulder,H_shoulder,H_headroom,H_headroom,Haisle,Haisle]
        #        plt.plot(plot_array_x,plot_array_y,"ro")
        #        circle = plt.Circle((r[0],r[1]), d-inner, color='b', fill=False)
        #        plt.show()

        # Several length of the fuselage
        lcockpit = 4
        ltail = Ratio_tail_d_outer * d_inner
        l_fuselage = lcabin + lcockpit + ltail
        lcabin_below = lcabin
        lcabin_above = 0
        N_aisle_above = 0
        N_aisle_below = N_aisles
        N_rows_above = 0
        N_rows_below = ceil(Npax_below/N_sa_below)

    # Same applies for a double decker configuration
    else:
        tot_seating_abreast = np.zeros([2, 3])
        Npax_above = 450 - Npax_below
        N_sa_below = round(0.45 * sqrt(Npax_below))
        N_sa_above = round(0.45 * sqrt(Npax_above))
        N_sa = [N_sa_below, N_sa_above]
        for i in range(len(N_sa)):
            if N_sa[i] <= 6:
                k_cabin = 1.08
                N_aisles = 1
                N_right = ceil(N_sa[i] / 2)
                N_left = N_sa[i] - N_right
                N_middle = 0

            elif N_sa[i] >= 6 and N_sa[i] <= 12:
                k_cabin = 1.17
                N_aisles = 2

                if N_sa[i] < 9:
                    N_left = 2
                    N_right = 2

                else:
                    N_left = 3
                    N_right = 3
                N_middle = N_sa[i] - N_left - N_right

            else:
                print("This aircraft have more than 2 aisle")
                N_right = 0
                N_left = 0
                N_middle = 0
                k_cabin = 0
                N_sa = 0
                N_aisles = 0

            if i == 0:
                tot_seating_abreast[1, 0] = N_right
                tot_seating_abreast[1, 1] = N_middle
                tot_seating_abreast[1, 2] = N_left

                if N_middle == 0:
                    N_aisle_below = 1
                    Nos_ce_below = 0

                elif N_middle != 0:
                    N_aisle_below = 2
                    Nos_ce_below = 1

            elif i == 1:
                tot_seating_abreast[0, 0] = N_right
                tot_seating_abreast[0, 1] = N_middle
                tot_seating_abreast[0, 2] = N_left

                if N_middle == 0:
                    N_aisle_above = 1
                    Nos_ce_above = 0

                elif N_middle != 0:
                    N_aisle_above = 2
                    Nos_ce_above = 1

        # Compute the length of the different compartments
        lcabin_below = (Npax_above / N_sa_above) * k_cabin
        lcabin_above = (Npax_below / N_sa_below) * k_cabin
        #compute Number of rows
        N_rows_above = ceil(Npax_above / N_sa_above)
        N_rows_below = ceil(Npax_below/N_sa_below)
        lcabin = max(lcabin_below, lcabin_above)
        W_cabin_above = N_sa_above * Wseat + (
                    N_sa_above + N_aisle_above + 1) * Warmrest + N_aisle_above * Waisle + 2 * Sclearance
        W_cabin_below = N_sa_below * Wseat + (
                    N_sa_below + N_aisle_below + 1) * Warmrest + N_aisle_below * Waisle + 2 * Sclearance
        W_cabin = max(W_cabin_above, W_cabin_below)
        W_headroom_above = W_cabin_above - 2 * (Warmrest + Sclearance) - Wseat
        W_headroom_below = W_cabin_below - 2 * (Warmrest + Sclearance) - Wseat
        V_os_above = (Nos_lat * Aos_lat + Nos_ce_above * Aos_ce) * lcabin_above * Kos
        V_os_below = (Nos_lat * Aos_lat + Nos_ce_below * Aos_ce) * lcabin_below * Kos
        V_os = V_os_above + V_os_below
        # cabin points below
        W_cabin_left_below_x = -W_cabin_below / 2
        W_headroom_left_below_x = -W_headroom_below / 2
        W_cabin_right_below_x = W_cabin_below / 2
        W_headroom_right_below_x = W_headroom_below / 2
        Waisle_left_below_x = -Waisle / 2
        Waisle_right_below_x = Waisle / 2

        cabin_left_point_below = (W_cabin_left_below_x, 0)
        cabin_right_point_below = (W_cabin_right_below_x, 0)
        headroom_left_point_below = (W_headroom_left_below_x, H_headroom)
        headroom_right_point_below = (W_headroom_right_below_x, H_headroom)
        aisle_left_point_below = (Waisle_left_below_x, Haisle)
        aisle_right_point_below = (Waisle_right_below_x, Haisle)
        shoulder_left_point_below = (W_cabin_left_below_x, H_shoulder)
        shoulder_right_point_below = (W_cabin_right_below_x, H_shoulder)

        # cabin points above
        W_cabin_left_above_x = -W_cabin_above / 2
        W_headroom_left_above_x = -W_headroom_above / 2
        W_cabin_right_above_x = W_cabin_above / 2
        W_headroom_right_above_x = W_headroom_above / 2
        Waisle_left_above_x = -Waisle / 2
        Waisle_right_above_x = Waisle / 2
        floor_height_above = Haisle + 0.05 * Haisle

        cabin_left_point_above = (W_cabin_left_above_x, 0 + floor_height_above)
        cabin_right_point_above = (W_cabin_right_above_x, 0 + floor_height_above)
        headroom_left_point_above = (W_headroom_left_above_x, H_headroom + floor_height_above)
        headroom_right_point_above = (W_headroom_right_above_x, H_headroom + floor_height_above)
        aisle_left_point_above = (Waisle_left_above_x, Haisle + floor_height_above)
        aisle_right_point_above = (Waisle_right_above_x, Haisle + floor_height_above)
        shoulder_left_point_above = (W_cabin_left_above_x, H_shoulder + floor_height_above)
        shoulder_right_point_above = (W_cabin_right_above_x, H_shoulder + floor_height_above)
        
        print(cabin_left_point_above, cabin_left_point_below, cabin_right_point_below, cabin_right_point_above,
             headroom_left_point_below, headroom_left_point_above, headroom_right_point_below,
             headroom_right_point_above, aisle_left_point_below, aisle_left_point_above, aisle_right_point_below,
             aisle_right_point_above, shoulder_left_point_below, shoulder_right_point_below, shoulder_left_point_above,
             shoulder_right_point_above)
        # function that creates the most inner circle
        r = make_circle(
            [cabin_left_point_above, cabin_left_point_below, cabin_right_point_below, cabin_right_point_above,
             headroom_left_point_below, headroom_left_point_above, headroom_right_point_below,
             headroom_right_point_above, aisle_left_point_below, aisle_left_point_above, aisle_right_point_below,
             aisle_right_point_above, shoulder_left_point_below, shoulder_right_point_below, shoulder_left_point_above,
             shoulder_right_point_above])
        print(r)
        d_inner = 2 * r[2]
        d_outer = 1.045 * d_inner + 0.084
        # other dimension of the total fuselage
        lcockpit = 4
        ltail = Ratio_tail_d_outer * d_inner
        l_fuselage = lcabin + lcockpit + ltail

    l_tailcone_range = np.array([2.6 * d_outer, 4 * d_outer])
    l_nosecone_range = np.array([1.2 * d_outer, 2.5 * d_outer])
    print(W_cabin_below)
    print(W_cabin_above)
    ##plot fuselage cross-section

    V_luggage = m_luggage / rho_luggage
    V_cargo = m_cargo / rho_cargo
    Vcc = V_cargo + (V_luggage - V_os)
    lpax_below = N_rows_below * Pseat
    lpax_above = N_rows_above * Pseat
    Npax_tot_seats = N_rows_below * N_sa_below+ N_rows_above * N_sa_above
    diff = Npax_tot_seats - Npax
    
    tot_seating_abreast_last_row = tot_seating_abreast*1
    if diff>0:
        tot_seating_abreast_last_row = np.zeros([2, 3])
        if N_sa_below>=N_sa_above:
        
            N_sa_below_last_row = N_sa_below-diff
            
            if N_sa_below_last_row <= 6:
                k_cabin = 1.08
                N_aisles = 1
                N_right = ceil(N_sa_below_last_row / 2)
                N_left = N_sa_below_last_row - N_right
                N_middle = 0
    
            elif N_sa_below_last_row >= 6 and N_sa_below_last_row <= 12:
                k_cabin = 1.17
                N_aisles = 2
    
                if N_sa_below_last_row < 9:
                    N_left = 2
                    N_right = 2
    
                else:
                    N_left = 3
                    N_right = 3
                    N_middle = N_sa_below_last_row - N_left - N_right
            tot_seating_abreast_last_row[1, 0] = N_right
            tot_seating_abreast_last_row[1, 1] = N_middle
            tot_seating_abreast_last_row[1, 2] = N_left
        else:
            N_sa_above_last_row = N_sa_above-diff
            if N_sa_above_last_row <= 6:
                k_cabin = 1.08
                N_aisles = 1
                N_right = ceil(N_sa_above_last_row / 2)
                N_left = N_sa_above_last_row - N_right
                N_middle = 0
    
            elif N_sa_above_last_row >= 6 and N_sa_above_last_row <= 12:
                k_cabin = 1.17
                N_aisles = 2
    
                if N_sa_above_last_row < 9:
                    N_left = 2
                    N_right = 2
    
                else:
                    N_left = 3
                    N_right = 3
                    N_middle = N_sa_above_last_row - N_left - N_right
            tot_seating_abreast_last_row[0, 0] = N_right
            tot_seating_abreast_last_row[0, 1] = N_middle
            tot_seating_abreast_last_row[0, 2] = N_left
    
    # return output function:
    # inner diameter         [0]
    # outer diameter         [1]
    # length cabin           [2]
    # Length cabin below     [3]
    # Length cabin above     [4]
    # Length tailcone range  [5]
    # Length nosecone range  [6]
    # length fuselage        [7]
    # tot seating abreast    [8]
    # Number of aisles above [9]
    # Number of aisles below [10]
    # Number of rows above   [11]
    # Number of rows below   [12]
    # Length of the passengers below [13]
    # Length of the passengers above [14]
    # tot seating abreast last row [15]
    # Seat pitch                    [16]
    return (d_inner, d_outer, lcabin, lcabin_below, lcabin_above, l_tailcone_range, l_nosecone_range, l_fuselage,
            tot_seating_abreast, N_aisle_above, N_aisle_below,N_rows_above,N_rows_below,lpax_below,lpax_above,tot_seating_abreast_last_row,Pseat)


d_inner, d_outer, lcabin, lcabin_below, lcabin_above, l_tailcone_range, l_nosecone_range, l_fuselage, tot_seating_abreast, N_aisle_above, N_aisle_below,N_rows_above,N_rows_below,lpax_below,lpax_above,tot_seating_abreast_last_row,Pseat=fuselage_cross_section(450, 242)


