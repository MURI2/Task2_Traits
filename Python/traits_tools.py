from __future__ import division
import os
import numpy as np

def get_path():
    return os.path.expanduser("~/GitHub/Task2_Traits")

def get_time(df):
    time_split = [list(map(int, x.split(':')))for x in df['Time'].values]
    time_hours = []
    for x in time_split:
        x_hour = ((x[0] * 60) + x[1]) / 60
        time_hours.append(round(x_hour, 2))
    t =  np.asarray(time_hours)
    return t

def remove_unnecessary_columns(df):
    to_remove = ['blank', 'Time', 'Temp_C']
    return df.loc[:, ~df.columns.str.contains('|'.join(to_remove))]
