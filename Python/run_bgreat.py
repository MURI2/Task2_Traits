from __future__ import division
import os, re
import numpy as np
import pandas as pd

mydir = os.path.expanduser("~/GitHub/Task2_Traits/")

def checkTemp(df):
    temp_min = min(df['Temp_C'])
    temp_max = max(df['Temp_C'])
    temp_diff = temp_max - temp_min
    if temp_diff > 3:
        print "Temperature difference greater than 3C, check for temperature effects"

# clean Data
def clean_data(path_IN, path_OUT, wells = 48):
    IN = open(path_IN, 'r')
    OUT = open(path_OUT, 'w')
    path_name = path_IN[:-4] + '_wellNames' + path_IN[-4:]
    IN_wells = pd.read_csv(path_name, sep = '\t')
    IN_wells["SampleRep"] = IN_wells['Sample'] + '_' + IN_wells["Rep"].map(str)
    names_dict = dict([(i, a) for i, a in zip(IN_wells.Well, IN_wells.SampleRep)])
    for line in IN:
        line_clean = line.strip().split('\t')
        if len(line_clean) == wells + 2:
            line_clean
            if line_clean[0] == 'Time':
                line_clean[1] = 'Temp_C'
                for j, item in enumerate(line_clean):
                    if item in names_dict:
                        line_clean[j] = names_dict[item]
            print>> OUT, '\t'.join(line_clean)
    OUT.close()

    df = pd.read_csv(path_OUT, sep = '\t')
    cols = [c for c in df.columns if ('empty_' not in c) and ('blank_media_' not in c) and ('blank_water_' not in c)]
    df = df[cols]
    df.to_csv(path_OUT, sep = '\t', index = False)

    #checkTemp

#def merge_clean_data():

def run_analysis():
    raw_data_path = mydir + 'data/uMax/raw_data/'
    clean_data_path = mydir + 'data/uMax/clean_data/'
    dfs = []
    for x in os.listdir(raw_data_path):
        if x == '.DS_Store':
            continue
        #if x != 'Task2_48hr_48well_discon_170928_135419':
        #    continue
        path_in = raw_data_path + x + '/' + x + '.txt'
        path_out = clean_data_path + x + '.txt'
        #clean_data(path_in, path_out, wells = 48)
        df = pd.read_csv(path_out, sep = '\t')
        # round to Minutes
        time_split = [map(int, x.split(':'))for x in df['Time'].values]
        time_hours = []
        # divide by 60 (get hours)
        for x in time_split:
            x_hour = ((x[0] * 60) + x[1]) / 60
            time_hours.append(round(x_hour, 2))
        t =  np.asarray(time_hours)
        df['Time'] = t
        df = df.drop('Temp_C', 1)
        dfs.append(df)

    data = reduce(lambda x, y: pd.merge(x, y, on = 'Time', how='outer'), dfs)
    #  data (n x p)
    print data

    # merge metadata



# merge cleaned data


# get just B and S lines

# run analysis

run_analysis()
