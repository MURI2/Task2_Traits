from __future__ import division
import pandas as pd
import numpy as np
import traits_tools as tt
import gpr_tools as gt
import bgreat
import matplotlib.pyplot as plt


df_path = tt.get_path() + '/data/uMax/clean_data/Task2_48hr_48well_discon_180614_140037.txt'
df = pd.read_csv(df_path, sep = '\t')
df.index = tt.get_time(df)
data = tt.remove_unnecessary_columns(df)

data = data.loc[:, data.columns.str.contains('|'.join(['L0B2-500', 'L2B5-500']))]
data = data.loc[:, ~data.columns.str.contains('|'.join(['L0B2-500_1_1', 'L0B2-500_3_2']))]
#print(data)
#del data['L0B2-500_1_1', 'L0B2-500_3_2']


strains = [x.split('_')[0] for x in data.columns.values ]
bio_reps = [x.split('_')[1] for x in data.columns.values ]
tech_reps = [x.split('_')[2] for x in data.columns.values ]
condition = ['control' for x in data.columns.values ]

meta = pd.DataFrame(
    {'strain': strains,
    'Condition': condition})

#print(data.columns.values)
#print(meta.strain.values)

parent = 'L0B2-500'
control = 'control'
mutants = ['L2B5-500']

meta['strain-regression'] = (meta.strain!=parent).astype(int)
meta['condition'] = (meta.Condition!=control).astype(int)
meta['interaction'] = meta['strain-regression']*meta.condition

bgreat.setGlobals(_data=data,_meta=meta)
bgreat.setGlobals(_parent=parent,_control=control)

#results = bgreat.testMutantControl(mutants,numPerm=20,dims=['time','strain-regression'])


bgreat.condition = None
gp = bgreat.buildGP(bgreat.selectStrain('L2B5-500'))

xpred = np.zeros((100,2))
xpred[:50,0] = np.linspace(data.index.min(),data.index.max())
xpred[50:,0] = np.linspace(data.index.min(),data.index.max())

xpred[50:,1] = 1

mu,cov = gp.predict(xpred,full_cov=True)
var = np.diag(cov)
mu = mu[:,0]



fig = plt.figure()
plt.plot(xpred[:50,0],mu[:50],label=parent);
plt.fill_between(xpred[:50,0],mu[:50]-2*np.sqrt(var[:50]),mu[:50]+2*np.sqrt(var[:50]),alpha=.1)

plt.plot(xpred[:50,0],mu[50:],label='L2B5-500')
plt.fill_between(xpred[:50,0],mu[50:]-2*np.sqrt(var[50:]),mu[50:]+2*np.sqrt(var[50:]),alpha=.1)

plt.legend(fontsize=20)
fig_name = tt.get_path() + '/figures/uMax/test_gpr.png'
fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()



op = np.zeros((50,100))
op[np.arange(50),np.arange(50)] = -1
op[np.arange(50),np.arange(50)+50] = 1


oddeltaMu = np.dot(op,mu)
oddeltaCov = np.dot(np.dot(op,cov),op.T)
oddeltaVar = np.diag(oddeltaCov)

print(oddeltaMu)



fig = plt.figure()
#plt.figure(figsize=(10,5))
plt.plot([data.index.min(),data.index.max()],[0,0],'k',lw=3)
plt.plot(xpred[:50,0],oddeltaMu,label='od-delta')
plt.fill_between(xpred[:50,0],oddeltaMu-2*np.sqrt(oddeltaVar),oddeltaMu+2*np.sqrt(oddeltaVar),alpha=.1)
plt.legend(fontsize=20)
fig_name = tt.get_path() + '/figures/uMax/test_od_delta.png'
fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()




# work on getting  umax....
