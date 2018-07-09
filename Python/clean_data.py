from __future__ import division
import os, re, pathlib
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
import traits_tools as tt


def clean_data(path_IN, path_OUT, wells = 48):
    # damn plate reader exports non utf characters
    IN = open(path_IN, 'r', errors='ignore')
    OUT = open(path_OUT, 'w')
    path_name = path_IN[:-4] + '_wellNames' + path_IN[-4:]
    IN_wells = pd.read_csv(path_name, sep = '\t')
    IN_wells["SampleRep"] = IN_wells['Sample'] + '_' + IN_wells["Bio_Rep"].map(str) + '_' + IN_wells["Tech_rep"].map(str)
    names_dict = dict([(i, a) for i, a in zip(IN_wells.Well, IN_wells.SampleRep)])
    for line in IN:
        line_clean = line.strip().split('\t')
        if len(line_clean) == wells + 2:
            if line_clean[0] == 'Time':
                line_clean[1] = 'Temp_C'
                for j, item in enumerate(line_clean):
                    if item in names_dict:
                        line_clean[j] = names_dict[item]
            #print>> OUT, '\t'.join(line_clean)
            OUT.write('\t'.join(line_clean) + '\n')
    OUT.close()


def checkTemp(df):
    temp_min = min(df['Temp_C'])
    temp_max = max(df['Temp_C'])
    temp_diff = temp_max - temp_min
    if temp_diff > 3:
        print("Temperature difference greater than 3C, check for temperature effects")





def plot_data(clean_data_path):
    out_path = tt.get_path() + '/figures/uMax/data_plots/' + re.split('[./]' , clean_data_path)[-2]
    pathlib.Path(out_path).mkdir(exist_ok=True)
    df = pd.read_csv(clean_data_path, sep = '\t')
    time_hours = get_time(df)
    for column in df:
        if ('blank' in column) or ('Time' in column) or ('Temp_C' in column):
            continue
        od = df[column].values
        fig = plt.figure()
        plt.scatter(time_hours, od, c='#87CEEB', marker='o', label='_nolegend_', s = 60)
        plt.title(column, fontsize = 24)
        plt.xlabel('Time (hours)', fontsize = 18)
        plt.ylabel('OD 600', fontsize = 18)
        fig_name = out_path + '/' + column + '.png'
        fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4)#, dpi = 600)
        plt.close()




in_path = tt.get_path() + '/data/uMax/raw_data/Task2_48hr_48well_discon_180614_140037/Task2_48hr_48well_discon_180614_140037.txt'
out_path = tt.get_path() + '/data/uMax/clean_data/Task2_48hr_48well_discon_180614_140037.txt'

#cleanData(in_path, out_path)
plot_data(out_path)
