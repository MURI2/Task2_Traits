from __future__ import division
import pandas as pd
import numpy as np
import scipy.signal
import scipy.stats as stats
from statsmodels.base.model import GenericLikelihoodModel
import  matplotlib.pyplot as plt

import os, re

np.random.seed(123456789)

mydir = os.path.expanduser("~/GitHub/Task2_Traits/")

# Modified Gompertz Equation
def m_gop(t, b0, A, umax, L):
    t = np.asarray(t)
    term = -np.exp(((umax*np.exp(1)) / A) * (L - t) + 1)
    return b0 + (A*np.exp(term))

# function to generate confidence intervals based on Fisher Information criteria
def CI_FIC(results):
    # standard errors = square root of the diagnol of a variance-covariance matrix
    ses = np.sqrt(np.absolute(np.diagonal(results.cov_params())))
    cfs = results.params
    lw = cfs - (1.96*ses)
    up = cfs +(1.96*ses)
    return (lw, up)


class modifiedGompertz(GenericLikelihoodModel):
    def __init__(self, endog, exog, **kwds):
        super(modifiedGompertz, self).__init__(endog, exog, **kwds)
        #print len(exog)

    def nloglikeobs(self, params):
        b0 = params[0]
        A = params[1]
        umax = params[2]
        L = params[3]
        z = params[4]
        # probability density function (pdf) is the same as dnorm
        exog_pred = m_gop(self.endog, b0 = b0, A = A, umax = umax, L = L)
        # need to flatten the exogenous variable
        LL = -stats.norm.logpdf(self.exog.flatten(), loc=exog_pred, scale=np.exp(z))
        return LL

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, method="bfgs", **kwds):

        if start_params is None:
            b0_start = 1
            A_start = 2
            umax_start = 0.5
            L_start = 0.8
            z_start = 0.8

            start_params = np.array([b0_start, A_start, umax_start, L_start, z_start])

        return super(modifiedGompertz, self).fit(start_params=start_params,
                                maxiter=maxiter, method = method, maxfun=maxfun,
                                **kwds)

def cleanData(path_IN, path_OUT, wells = 48):
    #IN = open(path_IN, 'r')
    OUT = open(path_OUT, 'w')
    path_name = path_IN[:-4] + '_wellNames' + path_IN[-4:]
    IN_wells = pd.read_csv(path_name, sep = '\t')
    IN_wells["SampleRep"] = IN_wells['Sample'] + '_' + IN_wells["Bio_Rep"].map(str) + '_' + IN_wells["Tech_rep"].map(str)
    names_dict = dict([(i, a) for i, a in zip(IN_wells.Well, IN_wells.SampleRep)])
    #print(IN)
    with  open(path_IN, 'r') as f:
        contents = f.read()
        print(contents)
        for line in IN:
            print(line)
            line_clean = line.strip().split('\t')
            print(line_clean)
            if len(line_clean) == wells + 2:
                line_clean
                if line_clean[0] == 'Time':
                    print(line_clean)
                    line_clean[1] = 'Temp_C'
                    for j, item in enumerate(line_clean):
                        if item in names_dict:
                            line_clean[j] = names_dict[item]
                print>> OUT, '\t'.join(line_clean)
        OUT.close()


def checkTemp(df):
    temp_min = min(df['Temp_C'])
    temp_max = max(df['Temp_C'])
    temp_diff = temp_max - temp_min
    if temp_diff > 3:
        print("Temperature difference greater than 3C, check for temperature effects")



#def modGompGrowth(IN_file_name, interceptGuess=0.1, delta = 0.05, synergy=True, \
#    temp = True, smooth = True):
#    IN = pd.read_csv(IN_file_name, sep = '\t')
#    time_split = [map(int, x.split(':'))for x in IN['Time'].values]
#    time_hours = []
#    for x in time_split:
#        x_hour = ((x[0] * 60) + x[1]) / 60
#        time_hours.append(round(x_hour, 2))
#    #IN['Time'] = pd.to_datetime(IN['Time'], format='%H:%M:%S')
#    #IN['Minutes'] = IN['Time'].dt.hour * 60 + IN['Time'].dt.minute
#    IN_wells_name = IN_file_name.replace('clean_data', 'raw_data')
#    split_in = re.split(r'[./]', IN_wells_name)#[-2]
#    split_in[-1] = split_in[-2]
#    merged_in = '/'.join(split_in)
#    IN_wells_name = merged_in + '_wellNames.txt'
#    IN_wells = pd.read_csv(IN_wells_name, sep = '\t')
#    IN_wells["SampleRep"] = IN_wells['Sample'] + '_' + IN_wells["Rep"].map(str)
#    cutoff_dict = dict([(i, a) for i, a in zip(IN_wells.SampleRep, IN_wells.x_axis_cutoff)])
#    if len(set(IN_wells.Keep.values)) == 1 and list(set(IN_wells.Keep.values))[0] == False:
#        print("No usable data")
#    else:
#        to_remove = IN_wells.loc[IN_wells['Keep'] == False].SampleRep.values
#        IN_file_name_split =  re.split(r'[./]', IN_file_name)
#        IN_file_name_split[-3] = 'params'
#        IN_file_name_split[-2] = IN_file_name_split[-2].replace('_clean', '_params')
#        path_OUT = '/'.join(IN_file_name_split[:-1]) + '.' + IN_file_name_split[-1]
#        OUT = open(path_OUT, 'w')
#        print>> OUT, 'Sample', 'Rep', 'b0', 'A', 'umax', 'L', 'z', 'umax_lw', 'umax_up', \
#                'umax_lw_FI', 'umax_up_FI'
#        #t = IN['Minutes'].values / 60
#        t =  np.asarray(time_hours)
#        to_keep = []
#        ignore_columns = ['Time', 'Temp_C', 'Minutes', 'blank_water_1', 'blank_water_2', 'blank_media_1']
#        for column in IN:
#            if IN[column].name in ignore_columns or IN[column].name in to_remove:
#                continue
#            #if IN[column].name != 'L1D1-100_2.1':
#            #    continue
#            #growth curve
#            s = IN[column].values
#            s_name = IN[column].name
#            print s_name
#            t_column = [t_i for t_i in t if t_i < cutoff_dict[s_name]]
#            s = s[:len(t_column)]
#            if max(s) - min(s) < delta:
#                print "Observed change in OD is not greater than " +  str(delta) + \
#                      " in well " + IN[column].name
#                continue
#            if smooth == True:
#                taps = [1/11] * 11
#                s_2 = scipy.signal.lfilter(taps, 1.0, s)
#                s_2[0:5] = s[0:5]
#                s_2[5:-5] = s_2[10:]
#                s_2[-6:] = s[-6:]
#            else:
#                s_2 = s
#            s_2_max =  np.argmax(s_2)
#            if len(s_2) < s_2_max + 20:
#                t_trim = t_column
#                s_trim = s_2
#                #print t
#                #print t_trim
#            else:
#                #t_trim = t[:s_2_max + 10]
#                #s_trim = s_2[:s_2_max + 10]
#                t_trim = t_column
#                s_trim = s_2

#            # we're going to loop through the following combination of
#            # parameter values
#            umax_start_list = [0.05,0.1,1]
#            L_start_list = [-5,-0.5,0.1,5,10,20]
#            z_start_list = [-2,-0.5]
#            # and while keeping the following initial values constant
#            b0_start = interceptGuess
#            A_start = max(s_trim)
#            #A_start = 1.5
#            def loop_params(_model, _umax_start_list, _L_start_list, _z_start_list, _b0_start, _A_start):
#                _results = []
#                for _umax_start in _umax_start_list:
#                    for _L_start in _L_start_list:
#                        for _z_start in _z_start_list:
#                            # b0, A, umax, L, z
#                            start_params = [_b0_start, _A_start, _umax_start, _L_start, _z_start]
#                            _result = model.fit(start_params = start_params, method="lbfgs", \
#                                bounds= [(-1,1), (_A_start * 0.66,10), (0,5), (-20,20), (-20, 20)], \
#                                disp = False)
#                            _results.append(_result)
#                AICs = [_result.aic for _result in _results]
#                _best = _results[AICs.index(min(AICs))]
#                return _best
#            while True:
#                model = modifiedGompertz(t_trim, s_trim)
#                best = loop_params(model, umax_start_list, L_start_list, z_start_list, b0_start, A_start)
#                if best.endog[-1] + 1 <= cutoff_dict[s_name]:
#                    print("Re-fitting the model")
#                else:
#                    break

#            best_CI_FIC = CI_FIC(best)
#            best_CI = best.conf_int()
#            best_params = best.params

            #print best.mle_settings
            #boot_mean, boot_std, boot_samples = best.bootstrap(nrep=500, store=True)
#            line = IN[column].name.split('_')[0]
#            rep = IN[column].name.split('_')[1]
#            print>> OUT, line, rep, best_params[0], best_params[1], best_params[2], \
#                    best_params[3], best_params[4], best_CI[2][0], best_CI[2][1], \
#                    best_CI_FIC[0][2], best_CI_FIC[1][2]

#            fig = plt.figure()
#            plt.plot(best.endog, best.exog, c = 'black', lw = 2)
#            y_pred = m_gop(best.endog, best_params[0], best_params[1], best_params[2], best_params[3])
#            plt.plot(best.endog, y_pred, c = 'blue', lw = 2)
#            #fig.tight_layout(pad = 0.5)
#            fig_direct =  re.split(r'[./]',path_OUT)
#            #fig_direct[-4] = 'figures/uMax'
#            #fig_direct[-3] = 'model_fits'
#            fig_path = mydir + 'figures/uMax/model_fits/' + fig_direct[-2]
#            #fig_direct.remove(fig_direct[-1])
#            #fig_direct = '/'.join(fig_direct)
#            if not os.path.exists(fig_path):
#                os.makedirs(fig_path)
#            fig_name = fig_path + '/' + IN[column].name + '.png'
#            fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
#            plt.close()

#        OUT.close()
