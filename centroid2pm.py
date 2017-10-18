#!/usr/bin/env python2
# vim: set ft=python

import numpy as np
import matplotlib.pyplot as plt
from re import search
from matplotlib.backends.backend_pdf import PdfPages
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from itertools import chain

cent = [i.split('\t') for i in open('./centroids_jg_first224.2dg').\
        read().split('\n')[:-1]]
cenh = cent[0]
cent = [[float(j) if not search('[a-zA-Z]', j) else np.nan if j == 'nan'
         else j for j in i] for i in cent[1:]]

pdf = PdfPages('centroid_regression_jg.pdf')

def pm(obj, cntrd):
    mjd = np.mean(cntrd[...,-2:], 1)
    # design matrix
    Ax = np.matrix([np.ones(np.shape(cntrd)[0]), mjd]).T
    Ay = np.matrix([np.ones(np.shape(cntrd)[0]), mjd]).T

    b_x = cntrd[...,0]
    b_y = cntrd[...,1]

    res_x = sm.OLS(b_x, Ax).fit()
    prstd_x, ci_l_x, ci_u_x = wls_prediction_std(res_x)
    res_y = sm.OLS(b_y, Ay).fit()
    prstd_y, ci_l_y, ci_u_y = wls_prediction_std(res_y)

    fig, pan = plt.subplots(1,2)
    pan[0].scatter(mjd, cntrd[...,0], color='deepskyblue')
    pan[0].plot(mjd, res_x.fittedvalues, color='darkorange')
    pan[0].plot(mjd, ci_u_x, '--', color='sandybrown')
    pan[0].plot(mjd, ci_l_x, '--', color='sandybrown')
    pan[0].set_title('x coordinates')
    pan[0].set_xlabel('MJD')
    pan[1].scatter(mjd, cntrd[...,1], color='deepskyblue')
    pan[1].plot(mjd, res_y.fittedvalues, color='darkorange')
    pan[1].plot(mjd, ci_u_y, '--', color='sandybrown')
    pan[1].plot(mjd, ci_l_y, '--', color='sandybrown')
    pan[1].set_title('y coordinates')
    pan[1].set_xlabel('MJD')
    fig.suptitle(obj)
    fig.tight_layout()
    pdf.savefig()
    fig.show()

    return (-(res_x.params[1] * 2.75) * 365, (res_y.params[1] * 2.75) * 365,
            np.sqrt(((res_x.params[1] * 2.75) * 365)**2 +\
                    ((res_y.params[1] * 2.75) * 365)**2),
            res_x.params, res_x.pvalues, res_x.rsquared, res_x.rsquared_adj,
            res_y.params, res_y.pvalues, res_y.rsquared, res_y.rsquared_adj)

def prepCntrd(x):
    cen = np.array(x[4:], dtype='float32').reshape((len(x)-4)/6,6)
    cen = cen[np.sum(np.isnan(cen), 1) < 1]
    return cen

test = pm(cent[12][0], prepCntrd(cent[12]))

for obj in cent:
