import pandas
import os
def read_input(design):
    ROOT_DIR =os.path.dirname(os.path.abspath("read_csv_input.py"))
    input_path = ROOT_DIR + '/input_files'
    df = pandas.read_csv(input_path+'/csv_input.csv', delimiter=';', index_col="parameter")
    #print(df)
    coef_series = df[design][0:4]
#    print(coef_series)
    char_series = df[design][4:9]
#    print(char_series)
    cruise_series = df[design][9:12]
#    print(cruise_series)
    engine_series = df[design][12:16]
    opt_series = df[design][16:27]
    tail_series = df[design][27:33]
#    print(engine_series)
#    print(opt_series)
#    print(tail_series)
    
    coefficients = pandas.Series.tolist(coef_series)
    ac_characteristics = pandas.Series.tolist(char_series)
    cruise_conditions = pandas.Series.tolist(cruise_series)
    engine = pandas.Series.tolist(engine_series)
    options = pandas.Series.tolist(opt_series)
    tail = pandas.Series.tolist(tail_series)
    return(coefficients, ac_characteristics,cruise_conditions, engine,options,tail)
#print(read_input('Design 1'))

