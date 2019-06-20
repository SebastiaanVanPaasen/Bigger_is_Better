#import csv
#
#f = open('CL1.csv')
#content = csv.reader(f)
#numbers = []
#for num in content:
#    


import pandas as pd
def importpan(textfile):
    df = pd.read_csv('CL1.csv')
    df.columns = ["Values"]

    tot_values = []
    for i in range(len(df)):
        values = df.iloc[i]['Values']
        tot_values.append(values)
    return(tot_values)
    