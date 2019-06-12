import pandas
import os

def read_input(design):
    ROOT_DIR =os.path.dirname(os.path.abspath("read_csv_input.py"))
    input_path = ROOT_DIR + '/input_files'
#    print(input_path)
    df = pandas.read_csv(input_path+'/csv_input_HIGH_2E_DD_STRUT_2.csv', delimiter=',', index_col="Parameter")
#    print(df)
    
    coef_series = df[design][0:4]
#    print(coef_series)
    char_series = df[design][4:9]
#    print(char_series)
    cruise_series = df[design][9:12]
#    print(cruise_series)
    engine_series = df[design][12:16]
    opt_series = df[design][16:35]
    tail_series = df[design][35:42]
    
    choices = [pandas.Series.tolist(opt_series[0:6]), pandas.Series.tolist(opt_series[6:8]),
               [opt_series[8]], [opt_series[9]], pandas.Series.tolist(opt_series[10:12]), [opt_series[12]], 
               [opt_series[13]], [opt_series[14]], pandas.Series.tolist(opt_series[15:17]), [opt_series[17]],
               [opt_series[18]]]
    
#    print(engine_series)
#    print(opt_series)
#    print(tail_series)
    
    coefficients = pandas.Series.tolist(coef_series)
    ac_characteristics = pandas.Series.tolist(char_series)
    cruise_conditions = pandas.Series.tolist(cruise_series)
    engine = pandas.Series.tolist(engine_series)
    tail = pandas.Series.tolist(tail_series)

    return(coefficients, ac_characteristics,cruise_conditions, engine,choices,tail)  

import sys


def read_output(design):
#    sys.path.append("C:/Users/Hidde/Documents/Aerospace Engineering Bachelor Year 3/Design Synthesis Exercise/Bigger_is_Better")
#    sys.path.append("C:/Users/mathi/Documents/DSE/Bigger_is_Better")
    sys.path.append("C:/Users/sebas/OneDrive/Documents/DSE/Bigger_is_Better")
#    ROOT_DIR = os.path.dirname(os.path.abspath("excel_results"))
    input_path = sys.path[-1] + '/excel_results'
#    print(input_path)
    
    df = pandas.read_csv(input_path+'/combined_results.csv', delimiter=',', index_col="Parameter")
#    print(df)
    
    weights = {"W_TO" : df[design][0], "W_F" : df[design][2], "W_W" : df[design][4], "W_N" : df[design][7], "W_E" : df[design][38]}
#    print(weights)
    
    wing = {"A" : df[design][15], "S" : df[design][16], "b" : df[design][17], "C_root" : df[design][18], "C_tip" : df[design][19], "Sweep" : df[design][20], "Taper" : df[design][21]}
#    print(wing)
    
    cruise_conditions = {"T_TO" : df[design][29], "V_cr" : df[design][31], "H_cr" : df[design][32], "CD_0" : df[design][37]}
#    print(cruise_conditions)

    return weights, wing, cruise_conditions

