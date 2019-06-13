import matplotlib.pyplot as plt
import numpy as np
#from class_I.class_I_empennage_landinggear import *
#from class_I.fuselage_cross_section import *
#from class_I.seats import *
#from input_files.conventional_double_decker import *
# ---------------INPUTS--------------------------------------------------------

#W_payload = 450000  # [N]
#W_person = 90  # mass of passenger and hand baggage kg
#
#cargo_fwdfrac = 0.6  # Fraction of the amount the front compartment holds, if N_cargo = 1 then the value is 1
#
## Weights from class II in [N]
#W_fuse = 80000 * 5
#W_nlg = 1500 * 5
#W_vt = 2500 * 5
#W_ht = 3500 * 5
#W_fix = 50000 * 5
#W_wing = 70000 * 5
#W_nac = 8000 * 5
#W_prop = 15000 * 5
#W_mlg = 6000 * 5
#
#W_fuel = 938780
#
#xcg_fuel = 0.55  # % MAC
#
#Safety_margin = 0.02  # for cg range
#
#d_inner, d_outer, lcabin, lcabin_below, lcabin_above, l_tailcone_range, l_nosecone_range, l_fuselage, tot_seating_abreast, N_aisle_above, N_aisle_below, N_rows_above, N_rows_below, lpax_below, lpax_above, tot_seating_abreast_last_row, Pseat = fuselage_cross_section(
#    N_pas, N_pas_below)
#xcg_seats = cg_seats(Pseat, N_rows_above, N_rows_below,l_nosecone_range,)
#W_window, W_aisle, W_middle = W_seats(tot_seating_abreast, N_rows_above, N_rows_below, tot_seating_abreast_last_row)
#D_fuse = d_outer
#l_nosecone = (l_nosecone_range[0] + l_nosecone_range[1]) / 2.
#xcg_payload = l_fuselage * 0.5  # cg-location payload w.r.t. nose [m]
## Import cg values from class I sizing
#
##cg_locations, S_h_frac_S, S_v_frac_S, X_LEMAC, x_h, x_v = class_I_empennage(MAC, l_fuselage, x_eng, l_n, xcg_OEW_MAC,
##                                                                       mass_frac_OEW, xcg_payload, mass_frac_payload,
##                                                                       xcg_fuel, mass_frac_fuel, D_fuse, b)
#
#xcg_nlg = l_nosecone  # assumption, use landing gear function !!!
#xcg_mlg = 0.75 * MAC  # assumption, use landing gear function !!!


# Inputs from fuselage dimensions

def potato(l_nosecone, W_TO, x_eng, l_n, mass_fractions, tail_h, tail_v, Safety_margin, cg_locations, X_LEMAC, constants, MAC):
    
    N_cargo, l_fuselage, cargo_fwdfrac, xcg_seats, W_window, W_aisle, W_middle, n_pax, W_person = constants
    
    x_le_h, sweep_LE_h, y_MAC_h, MAC_h, l_h, c_root_h, c_tip_h, b_h = tail_h
    x_le_v, sweep_LE_v, y_MAC_v, MAC_v, c_root_v, c_tip_v, b_v = tail_v
    
    W_wing = mass_fractions[0]*W_TO
    W_ht = W_vt = mass_fractions[1]*W_TO
    W_fuse = mass_fractions[2]*W_TO
    W_nac = mass_fractions[3]*W_TO
    W_prop = mass_fractions[4]*W_TO
    W_fix = mass_fractions[5]*W_TO
    W_mlg   = 0.8*mass_fractions[6]*W_TO
    W_nlg   = mass_fractions[6]*W_TO - W_mlg
    
    W_payload = mass_fractions[8]*W_TO
    W_fuel = mass_fractions[9]*W_TO
    
    xcg_wing = 0.4 * MAC
    xcg_prop = x_eng + 0.4 * l_n
    xcg_nac = x_eng + 0.4 * l_n

    # Fuselage group - in [m] w.r.t. nose
    xcg_fuse = 0.5 * l_fuselage
    xcg_fix = 0.5 * l_fuselage
    xcg_emp = 0.9 * l_fuselage

    xcg_nlg = l_nosecone  # assumption, use landing gear function !!!
    xcg_mlg = 0.75 * MAC   # assumption, use landing gear function !!!
    xcg_fuel = 0.55  # mac %
    
    
    cg_locations[0] = xcg_fuse 
    cg_locations[1] = xcg_emp  
    cg_locations[2] = xcg_fix  
    cg_locations[3] = xcg_nac  
    cg_locations[4] = xcg_prop  
    cg_locations[5] = xcg_wing
    cg_locations[8] = xcg_mlg
    cg_locations[9] = xcg_nlg
    
    # Fuselage group weight and cg
    W_fusegroup = W_fuse + W_nlg + W_vt + W_ht + W_fix
    xcg_vt = x_le_v + np.tan(sweep_LE_v)*y_MAC_v + 0.42*MAC_v
    xcg_ht = x_le_h + np.tan(sweep_LE_h)*y_MAC_h + 0.42*MAC_h
    xcg_fusegroup = (xcg_fuse * W_fuse + xcg_nlg * W_nlg + xcg_vt * W_vt + xcg_ht * W_ht + xcg_fix * W_fix) / W_fusegroup

    # Wing group weight and cg
    W_winggroup = W_wing + W_nac + W_prop + W_mlg
    xcg_winggroup = (xcg_wing * W_wing + xcg_nac * W_nac + xcg_prop * W_prop + xcg_mlg * W_mlg) / W_winggroup
    xcg_wg_MAC = xcg_winggroup / MAC
    plt.close()
    min_cg = []
    max_cg = []
    X_LEMAC_range = []
    xcg_OEW_MAC = 0
    # print(X_LEMAC/l_fuselage)
    for j in np.arange((X_LEMAC / l_fuselage) - 0.2, (X_LEMAC / l_fuselage) + 0.2, 0.001):
#    for j in np.arange(26.235227232293713/l_fuselage, 0.5, 100):
    
        X_LEMAC_range.append(j) 
        X_LEMAC = j * l_fuselage

        # Aircraft OEW weight and cg
        # xcg_OEW_MAC = (xcg_fusegroup + (W_winggroup/W_fusegroup)*xcg_wg_MAC - X_LEMAC)/(1. + W_winggroup/W_fusegroup)
        W_OEW = W_fusegroup + W_winggroup
        xcg_OEW_MAC = ((xcg_fusegroup - X_LEMAC) / MAC * W_fusegroup + xcg_wg_MAC * W_winggroup) / W_OEW

        # ------------Cargo loading-------------------------------------------------
        W_pax = n_pax * W_person
        W_cargo = W_payload - W_pax

        if N_cargo == 1:
            xcg_cargo1_MAC = (0.45 * l_fuselage - X_LEMAC) / MAC  # cg cargo as percentage of MAC
            xcg_cargo2_MAC = 0
            W_cargo1 = W_cargo
            W_cargo2 = 0

        else:
            xcg_cargo1_MAC = (0.3 * l_fuselage - X_LEMAC) / MAC  # cg cargo 1 as percentage of MAC
            xcg_cargo2_MAC = (0.75 * l_fuselage - X_LEMAC) / MAC  # cg cargo 2 as percentage of MAC
            W_cargo1 = W_cargo * (cargo_fwdfrac)
            W_cargo2 = W_cargo * (1 - cargo_fwdfrac)

            # --------------- Backward loading------------------------------------------

        # Weight order
        W_Orderback = [W_OEW] + [W_cargo2, W_cargo1] + list((np.array(W_window[::-1]) * W_person)) + list(
            (np.array(W_aisle[::-1]) * W_person)) + list((np.array(W_middle[::-1]) * W_person)) + [W_fuel]

        # cg order
        xcg_Orderback = [xcg_OEW_MAC] + [xcg_cargo2_MAC, xcg_cargo1_MAC] + 3 * list(
            (np.array(xcg_seats[::-1]) - X_LEMAC) / MAC) + [xcg_fuel]

        # Weight actual values
        W_Backloading = [W_OEW]
        for i in range(1, len(W_Orderback)):
            W_Backloading.append(W_Backloading[i - 1] + W_Orderback[i])

        # cg actual values
        xcg_Backloading = [xcg_OEW_MAC]
        for i in range(1, len(xcg_Orderback)):
            xcg_Backloading.append(
                (xcg_Backloading[i - 1] * W_Backloading[i - 1] + xcg_Orderback[i] * W_Orderback[i]) / W_Backloading[i])

        # ----------------Forward loading-------------------------------------------

        # Weight order
        W_Orderfwd = [W_OEW] + [W_cargo1, W_cargo2] + list((np.array(W_window) * W_person)) + list(
            (np.array(W_aisle) * W_person)) + list((np.array(W_middle) * W_person)) + [W_fuel]

        # cg order
        xcg_Orderfwd = [xcg_OEW_MAC] + [xcg_cargo1_MAC, xcg_cargo2_MAC] + 3 * list(
            (np.array(xcg_seats) - X_LEMAC) / MAC) + [xcg_fuel]

        # Weight actual values
        W_Fwdloading = [W_OEW]
        for i in range(1, len(W_Orderfwd)):
            W_Fwdloading.append(W_Fwdloading[i - 1] + W_Orderfwd[i])

        # cg actual values
        xcg_Fwdloading = [xcg_OEW_MAC]
        for i in range(1, len(xcg_Orderfwd)):
            xcg_Fwdloading.append(
                (xcg_Fwdloading[i - 1] * W_Fwdloading[i - 1] + xcg_Orderfwd[i] * W_Orderfwd[i]) / W_Fwdloading[i])
        
        W_fwd_cargo = W_Fwdloading[0:3]
        W_fwd_window = W_Fwdloading[2:(3+len(xcg_seats))]
        W_fwd_aisle = W_Fwdloading[(2+len(xcg_seats)):(3+2*len(xcg_seats))]
        W_fwd_middle = W_Fwdloading[(2+2*len(xcg_seats)):(3+3*len(xcg_seats))]
        W_fwd_fuel = W_Fwdloading[(2+3*len(xcg_seats)):(4+3*len(xcg_seats))]
        
        W_bwd_cargo = W_Backloading[0:3]
        W_bwd_window = W_Backloading[2:(3+len(xcg_seats))]
        W_bwd_aisle = W_Backloading[(2+len(xcg_seats)):(3+2*len(xcg_seats))]
        W_bwd_middle = W_Backloading[(2+2*len(xcg_seats)):(3+3*len(xcg_seats))]
#        W_bwd_fuel = W_Backloading[(2+3*len(xcg_seats)):(4+3*len(xcg_seats))]
        
        xcg_fwd_cargo = xcg_Fwdloading[0:3]
        xcg_fwd_window = xcg_Fwdloading[2:(3+len(xcg_seats))]
        xcg_fwd_aisle = xcg_Fwdloading[(2+len(xcg_seats)):(3+2*len(xcg_seats))]
        xcg_fwd_middle = xcg_Fwdloading[(2+2*len(xcg_seats)):(3+3*len(xcg_seats))]
        xcg_fwd_fuel = xcg_Fwdloading[(2+3*len(xcg_seats)):(4+3*len(xcg_seats))]
        
        xcg_bwd_cargo = xcg_Backloading[0:3]
        xcg_bwd_window = xcg_Backloading[2:(3+len(xcg_seats))]
        xcg_bwd_aisle = xcg_Backloading[(2+len(xcg_seats)):(3+2*len(xcg_seats))]
        xcg_bwd_middle = xcg_Backloading[(2+2*len(xcg_seats)):(3+3*len(xcg_seats))]
#        xcg_bwd_fuel = xcg_Backloading[(2+3*len(xcg_seats)):(4+3*len(xcg_seats))]
        
        minlist = [min(xcg_Fwdloading) - Safety_margin, min(xcg_Fwdloading) - Safety_margin]
        maxlist = [max(xcg_Backloading) + Safety_margin, max(xcg_Backloading) + Safety_margin]
#        print((max(xcg_Backloading) + Safety_margin) - (min(xcg_Fwdloading) - Safety_margin))
        W_range = [W_Fwdloading[0], W_Fwdloading[-1]]
            
        
        min_cg.append(min(xcg_Fwdloading) - Safety_margin)
        max_cg.append(max(xcg_Backloading) + Safety_margin)
        
        
#        plt.ylabel("$W$ [N]")
#        plt.xlabel("$x_{cg}$/$MAC$ [-]")
#        plt.plot(xcg_fwd_cargo, W_fwd_cargo, label="Cargo fwd")
#        plt.plot(xcg_fwd_window, W_fwd_window, label="Window fwd")
#        plt.plot(xcg_fwd_aisle, W_fwd_aisle, label="Aisle fwd")
#        plt.plot(xcg_fwd_middle, W_fwd_middle, label="Middle fwd")
#        plt.plot(xcg_bwd_cargo, W_bwd_cargo, label="Cargo bwd")
#        plt.plot(xcg_bwd_window, W_bwd_window, label="Window bwd")
#        plt.plot(xcg_bwd_aisle, W_bwd_aisle, label="Aisle bwd")
#        plt.plot(xcg_bwd_middle, W_bwd_middle, label="Middle bwd")
#        plt.plot(xcg_fwd_fuel, W_fwd_fuel, label = "Fuel")
#        plt.plot(minlist, W_range, "k")
#        plt.plot(maxlist, W_range, "k")
#        plt.ylim(round(min(W_range)-50000.,-5), round(max(W_range)+50000.,-5))
#        plt.yticks(np.arange(round(min(W_range)-50000.,-5), round(max(W_range)+160000.,-5), 100000.))
#        plt.xlim(-0.1, 0.6)
#        plt.xticks(np.arange(-0.1, 0.7, 0.1))
#        plt.legend(loc="upper left", prop={'size': 7}, bbox_to_anchor=(.2, .4))
#        plt.grid(True)
#        plt.show()
        
    difflist = []
    for j in range(len(min_cg)):
        difflist.append(abs(min_cg[j]-max_cg[j]))
        
    min_cg_range = min(difflist)

#    plt.plot(min_cg, X_LEMAC_range)
#    plt.plot(max_cg, X_LEMAC_range)
#    plt.show()
    
    return min_cg, max_cg, X_LEMAC_range, min_cg_range


#potato(Safety_margin, xcg_fuel, W_fuel, W_fuse, W_nlg, W_vt, W_ht, W_fix, W_wing, W_nac, W_prop, W_mlg, cg_locations,
#       xcg_nlg, xcg_mlg, X_LEMAC, W_payload, N_cargo, l_fuselage, MAC, cargo_fwdfrac, xcg_seats,
#       W_window, W_aisle, W_middle)
