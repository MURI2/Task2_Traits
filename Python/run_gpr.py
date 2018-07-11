from __future__ import division
import pandas as pd
import numpy as np
import traits_tools as tt
import bgreat
import time as timemodule
from gp_growth import factory, metric
from gp_growth.data.growth import GrowthData


def get_params(file_name):
    df_path = tt.get_path() + '/data/uMax/clean_data/Task2_48hr_48well_discon_180614_140037.txt'
    df = pd.read_csv(df_path, sep = '\t')
    df.index = tt.get_time(df)
    data = tt.remove_unnecessary_columns(df)

    data = data.loc[:, data.columns.str.contains('|'.join(['L0B2-500', 'L2B5-500']))]
    data = data.loc[:, ~data.columns.str.contains('|'.join(['L0B2-500_1_1', 'L0B2-500_3_2']))]

    strains = [x.split('_')[0] for x in data.columns.values ]
    bio_reps = [x.split('_')[1] for x in data.columns.values ]
    tech_reps = [x.split('_')[2] for x in data.columns.values ]
    condition = ['control' for x in data.columns.values ]
    meta = pd.DataFrame(
        {'strain': strains,
        'Bio': bio_reps,
        'Tech': tech_reps,
        'Condition': condition})
    #parent = 'L0B2-500'
    parent = strains[0]
    control = 'control'
    #mutants = ['L2B5-500']
    meta['strain-regression'] = (meta.strain!=parent).astype(int)
    meta['condition'] = (meta.Condition!=control).astype(int)
    meta['interaction'] = meta['strain-regression']*meta.condition
    bgreat.setGlobals(_data=data,_meta=meta)
    bgreat.setGlobals(_parent=parent,_control=control)
    bgreat.condition = None
    meta.rename(columns={'strain':'Strain'}, inplace=True)
    data = pd.concat([meta, data.transpose()], axis=1)
    # add index as well column for now
    data['well'] = data.index.values
    cols = list(data)
    cols.insert(0, cols.pop(cols.index('well')))
    data = data.ix[:, cols]

    gpFactory = factory.Factory()
    mse = metric.MeanSquareError(factory=gpFactory)
    muMax = metric.MuMax_simple(factory=gpFactory,n=100)
    #lagTime = metric.LagTime(f=gpFactory)
    lagTime = metric.LagTime(f =.5,factory=gpFactory,)

    carryingCapacity = metric.CarryingCapacity_simple(factory=gpFactory)

    output = pd.DataFrame()
    timeind = 8
    time = data.columns[timeind:].astype(float)
    for i in range(output.shape[0],data.shape[0]):
        od = np.log2(data.iloc[i,timeind:].tolist())
        od = od - od[0]
        design = params = data.iloc[i,:timeind]
        params = pd.DataFrame([params.values],columns=params.index,index=[params.name])
        key = pd.DataFrame(design).T
        temp = GrowthData(pd.DataFrame(list(zip(time.tolist(),od.tolist())),columns=['time',params.index[0]]),key)
        edata = temp.getData("gp")

        n = edata.shape[0]
        train = list(range(n))
        test = list(range(0,n,5))
        [train.remove(x) for x in test]
        cv = [train,test]
        edata_train = edata.iloc[train,:]
        edata_test = edata.iloc[test,:]

        # GP MSE
        try:
            start = timemodule.time()
            gp = gpFactory.build(edata_train, optimize=True)
            stop = timemodule.time()
            params['gp_mse'] = mse.compute(edata_test,gp)
            params['gp_muMax_mean'],params['gp_muMax_std'] = muMax.compute(predictive_data=temp.data.values,model=gp)
            params['gp_CarryingCapacity_mean'],params['gp_CarryingCapacity_std'] = carryingCapacity.compute(predictive_data=temp.data.values,model=gp)
            #params['gp_LagTime_mean'],params['gp_LagTime_std'] = lagTime.compute(predictive_data=temp.data.values,model=gp)
            #params['gp_LagTime_mean'] = lagTime.compute(gpFactory.buildInputFixed(time_min=0,time_max=max(edata.time),size=200,convert=False),gp)
            #print(lagTime.compute(predictive_data = test_data, model = gp1))

            params['gp_loglikelihood'] = gp.log_likelihood()
            params['gp_time'] = stop - start

        except Exception as e:
            params['gp_mse'] = [np.nan]
            params['gp_muMax_mean'],params['gp_muMax_std'] = [np.nan,np.nan]
            params['gp_CarryingCapacity_mean'],params['gp_CarryingCapacity_std'] = [np.nan,np.nan]
            #params['gp_LagTime_mean'],params['gp_LagTime_std'] = [np.nan,np.nan]
            params['gp_loglikelihood'] = [np.nan]
            params['gp_time'] = [np.nan]

        if output.shape[0]==0:
            output = params
        else:
            output = output.append(params)


    df_out = tt.get_path() + '/data/uMax/params_gpr/Task2_48hr_48well_discon_180614_140037.txt'
    output.to_csv(df_out, sep = '\t', index = False)




file_name = 'Task2_48hr_48well_discon_180614_140037'
out_df = get_params(file_name)
