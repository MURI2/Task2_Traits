
import pandas as pd
import os, math
mydir = os.path.expanduser("~/GitHub/Task2_Traits")
results = pd.read_csv(mydir + "/data/uMax/merged_clean_data/data.txt", sep = "\t", index_col=0)
#print(L0B3-100_2_170824_195818.iloc[:,:])

print(results.filter(regex = ("L?B?-300_*" )))

#%matplotlib inline
#results.plot(kind = "hist")