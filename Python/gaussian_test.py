#import numpy as np
#from matplotlib import pyplot as plt

import pandas as pd
data = pd.read_csv("~/GitHub/Task2_Traits/data/uMax/merged_clean_data.csv",index_col=0)
meta = pd.read_csv("~/GitHub/Task2_Traits/data/uMax/merged_clean_data.csv")

assert data.shape[1] == meta.shape[0]

import bgreat
bgreat.setGlobals(_data=data,_meta=meta)

parent = 'my-parent'
meta['strain-regression'] = (meta.strain!=parent).astype(int)
bgreat.setGlobals(_parent=parent,_control=control)

mutants = ['m1','m2']
results = testMutants(mutants)

control = 'my-control'
meta['interaction'] = ((meta.Condition!=control) & (meta.strain!=parent)).astype(int)
meta['condition'] = (meta.Condition!=control).astype(int)
bgreat.setGlobals(_control=control,_meta=meta)

results = bgreat.testMutantCondition(mutants)

#from sklearn.gaussian_process import GaussianProcessRegressor
#from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

#np.random.seed(1)

#def f(x):
    #"""The function to predict."""
    #return x * np.sin(x)