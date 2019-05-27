import pandas
import os

def read_input(design):
    ROOT_DIR =os.path.dirname(os.path.abspath("read_csv_input.py"))
    input_path = ROOT_DIR + '/input_files'
#    print(input_path)
    df = pandas.read_csv(input_path+'/csv_input_HIGH_2E.csv', delimiter=',', index_col="Parameter")
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


