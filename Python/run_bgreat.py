from __future__ import division
import os, re
import numpy as np
import pandas as pd
import  matplotlib.pyplot as plt
import scipy.signal
import bgreat as bg

mydir = os.path.expanduser("~/GitHub/Task2_Traits/")

def checkTemp(df):
    temp_min = min(df['Temp_C'])
    temp_max = max(df['Temp_C'])
    temp_diff = temp_max - temp_min
    if temp_diff > 3:
        print "Temperature difference greater than 3C, check for temperature effects"


def clean_file(path_IN, path_OUT, wells = 48):
    IN = open(path_IN, 'r')
    OUT = open(path_OUT, 'w')
    path_name = path_IN[:-4] + '_wellNames' + path_IN[-4:]
    IN_wells = pd.read_csv(path_name, sep = '\t')
    IN_wells["SampleRep"] = IN_wells['Sample'] + '_' + IN_wells["Rep"].map(str)
    names_dict = dict([(i, a) for i, a in zip(IN_wells.Well, IN_wells.SampleRep)])
    keep_dict = dict([(i, a) for i, a in zip(IN_wells.Well, IN_wells.Keep)])
    for line in IN:
        line_clean = line.strip().split('\t')
        if len(line_clean) == wells + 2:
            if line_clean[0] == 'Time':
                line_clean[1] = 'Temp_C'
                for j, item in enumerate(line_clean):
                    if (item in names_dict) and (keep_dict[item] == True):
                        line_clean[j] = names_dict[item]
                    elif (item in names_dict) and (keep_dict[item] == False):
                        line_clean[j] = 'Remove_' + names_dict[item]
                    else:
                        continue
            print>> OUT, '\t'.join(line_clean)
    OUT.close()

    df = pd.read_csv(path_OUT, sep = '\t')
    cols = [c for c in df.columns if ('empty_' not in c) and ('blank_media_' not in c) and ('blank_water_' not in c) and ('Remove_' not in c)]
    df = df[cols]
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
    df.to_csv(path_OUT, sep = '\t', index = False)

    #checkTemp


def clean_data(figs = False):
    raw_data_path = mydir + 'data/uMax/raw_data/'
    clean_data_path = mydir + 'data/uMax/clean_data/'
    dfs = []
    for x in os.listdir(raw_data_path):
        if x == '.DS_Store':
            continue
        print x
        path_in = raw_data_path + x + '/' + x + '.txt'
        path_out = clean_data_path + x + '.txt'
        clean_file(path_in, path_out, wells = 48)
        df = pd.read_csv(path_out, sep = '\t')
        #df['Time'] = df['Time'].apply(np.floor)
        # 194 time points

        # plot every pop
        t = df.Time.values
        if figs == False:
            continue
        fig_direct = re.split(r'[./]', path_out)
        fig_path = mydir + 'figures/uMax/data_plots/' + fig_direct[-2]
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        for column in df:
            if column == 'Time':
                continue
            od = df[column]
            fig = plt.figure()
            plt.plot(t, od, c = 'black', lw = 2)
            fig.tight_layout(pad = 0.5)
            fig_name = fig_path + '/' + df[column].name + '.png'
            fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
            plt.close()

def merge_data():
    clean_data_path = mydir + 'data/uMax/clean_data/'
    dfs = []
    for x in os.listdir(clean_data_path):
        if x == '.DS_Store':
            continue
        x_split = re.split(r'[._/]', x)
        add_on = x_split[-3] + '_' + x_split[-2]
        path_in = clean_data_path + x
        df = pd.read_csv(path_in, sep = '\t')
        df_column_dict = {}
        for column in df.columns.values:
            if column == 'Time':
                continue
            df_column_dict[column] = column + '_' + add_on
        df = df.rename(columns = df_column_dict)
        dfs.append(df)
    data = reduce(lambda x, y: pd.merge(x, y, on = 'Time', how='outer'), dfs)
    data_OUT = mydir + 'data/uMax/merged_clean_data/data.txt'
    data.to_csv(data_OUT, sep = '\t', index = False)

    data_columns = data.columns.values[1:]
    data_columns_dict = {}
    for data_column in data_columns:
        # [treatment, strain, rep, day, technical rep, plate]
        data_column_split = re.split(r'[-_/]', data_column)
        data_columns_dict[data_column] = [data_column_split[0][1], data_column_split[0][2], \
                                        data_column_split[0][3], data_column_split[1], \
                                        data_column_split[2], data_column_split[3]+'_'+data_column_split[4]]
    metadata = pd.DataFrame.from_dict(data_columns_dict).T
    metadata = metadata.reset_index(drop=False)
    metadata.columns = ['Sample', 'Treatment', 'Strain', 'BiolRep', 'Day', 'TechRep', 'Plate']
    metadata_OUT = mydir + 'data/uMax/merged_clean_data/metadata.txt'
    metadata.to_csv(metadata_OUT, sep = '\t', index = False)

def bgreat():
    data_path = mydir + 'data/uMax/merged_clean_data/data.txt'
    meta_path = mydir + 'data/uMax/merged_clean_data/metadata.txt'
    data = pd.read_csv(data_path, sep = '\t')
    meta = pd.read_csv(meta_path, sep = '\t')
    cols = [c for c in data.columns if (('B' in c) or ('S' in c)) and ('100' in c)]
    data = data[cols]
    meta = meta[meta['Sample'].isin(cols)]

    parent = 'B'
    control = 0
    condition = 0

    meta['strain-regression'] = (meta.Strain!=parent).astype(int)
    meta['condition'] = meta['Treatment']
    meta['interaction'] = meta['strain-regression']*meta.condition
    # try to add multiple
    data = np.log10(data)


    #parent = 'B'
    #meta['strain-regression'] = (meta.Strain!=parent).astype(int)
    #bgreat.setGlobals(_data=data,_meta=meta)
    #bgreat.setGlobals(_parent=parent,_control=control)

    #meta['Interaction'] = ((meta.Treatment!=control) & (meta.Strain!=parent)).astype(int)
    #meta['Condition'] = (meta.Treatment!=control).astype(int)
    #bgreat.setGlobals(_control=control,_meta=meta)

    #mutants = ['S']
    #results = bgreat.testMutants(mutants)
    #print results


    #results = bgreat.testMutantCondition(mutants)

def test_bgreat():
    data_path = mydir + 'data/uMax/merged_clean_data/test_data.csv'
    meta_path = mydir + 'data/uMax/merged_clean_data/test_meta.csv'
    data = pd.read_csv(data_path, sep = ',', index_col=0)
    meta = pd.read_csv(meta_path, sep = ',')

    parent = 'parent'
    control = 'control'
    condition = 'stress'

    meta['strain-regression'] = (meta.strain!=parent).astype(int)
    meta['condition'] = (meta.Condition!=control).astype(int)
    meta['interaction'] = meta['strain-regression']*meta.condition

    bg.setGlobals(_data=data,_meta=meta,_parent=parent,_control=control,_condition=condition)

    mutants = ['mutant']
    results = bg.testMutantCondition(mutants,numPerm=1)

    #bg.setGlobals(_data=data,_meta=meta)

    #parent = 'parent'
    #meta['strain-regression'] = (meta.strain!=parent).astype(int)
    #control = 'control'
    #bg.setGlobals(_parent=parent,_control=control)

    #mutants = ['mutant']
    #results = bg.testMutants(mutants)

    #meta['interaction'] = ((meta.Condition!=control) & (meta.strain!=parent)).astype(int)
    #meta['condition'] = (meta.Condition!=control).astype(int)
    #bg.setGlobals(_control=control,_meta=meta)

    #results = bg.testMutantCondition(mutants)
    #print results

def figure_test(smooth = True):
    data_path = mydir + 'data/uMax/merged_clean_data/data.txt'
    meta_path = mydir + 'data/uMax/merged_clean_data/metadata.txt'
    data = pd.read_csv(data_path, sep = '\t')
    meta = pd.read_csv(meta_path, sep = '\t')

    cols_B0 = [c for c in data.columns if ('B' in c) and ('L2' in c) and ('-100' in c)]
    data_B0 = data[cols_B0]

    cols_S0 = [c for c in data.columns if ('S' in c) and ('L2' in c) and ('-100' in c)]
    data_S0 = data[cols_S0]
    fig = plt.figure()
    for column in data_B0:
        column_zip = zip(data['Time'].values, data_B0[column].values)
        column_zip = [x for x in column_zip if np.isnan(x[1]) == False]
        t = [x[0] for x in column_zip]

        od = [x[1] for x in column_zip]

        plt.plot(t, od, c = 'blue', lw = 2, alpha = 0.7)

    for column in data_S0:
        column_zip = zip(data['Time'].values, data_S0[column].values)
        column_zip = [x for x in column_zip if np.isnan(x[1]) == False]
        t = [x[0] for x in column_zip]
        od = [x[1] for x in column_zip]

        plt.plot(t, od, c = 'black', lw = 2, alpha = 0.7)

    #fig = plt.figure()
    #fig_direct =  re.split(r'[./]',path_OUT)
    fig_name = mydir + 'figures/SSE_grant.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



#clean_data(figs = False)
#merge_data()
#bgreat()
#figure_test()
