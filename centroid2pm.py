#!/usr/bin/env python2
# vim: set ft=python

import numpy as np
import matplotlib.pyplot as plt
from re import search
from matplotlib.backends.backend_pdf import PdfPages
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std

cent = [i.split('\t') for i in open('./centroids_jg_first224.2dg').\
        read().split('\n')[:-1]]
cenh = cent[0]
cent = [[float(j) if not search('[a-zA-Z]', j) else np.nan if j == 'nan'
         else j for j in i] for i in cent[1:]]

ppcent = [i.split('\t') for i in open('./preposhib_centroids_jg_first224.2dg').\
          read().split('\n')[:-1]]

def pmAllEpoch(obj, cntrd, plot=True):
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

    if plot:
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
        plt.close('All')

    return [-(res_x.params[1] * 2.75) * 365, (res_y.params[1] * 2.75) * 365,
            np.sqrt(((res_x.params[1] * 2.75) * 365)**2 +\
                    ((res_y.params[1] * 2.75) * 365)**2)] +\
            [float(i) for i in np.nditer(res_x.params)] +\
            [float(i) for i in np.nditer(res_x.pvalues)] +\
            [float(i) for i in np.nditer(res_x.rsquared)] +\
            [float(i) for i in np.nditer(res_y.params)] +\
            [float(i) for i in np.nditer(res_y.pvalues)] +\
            [float(i) for i in np.nditer(res_y.rsquared)]

def prepCntrd(x):
    cen = np.array(x[4:], dtype='float32').reshape((len(x)-4)/6,6)
    cen = cen[np.sum(np.isnan(cen), 1) < 1]
    return cen

def pmPrePost(OBJ):
    d = [float(i) for i in OBJ[4:]]
    pm_ra = -((d[4] - d[0]) / (d[-1] - d[-2])) * 2.75 * 365
    pm_dec = ((d[5] - d[1]) / (d[-1] - d[-2])) * 2.75 * 365
    pm_tot = np.sqrt(pm_ra**2 + pm_dec**2)
    return [pm_ra, pm_dec, pm_tot]

pdf = PdfPages('centroid_regression_jg.pdf')

with open('jg_first224_pm_all', 'wb') as out:
    out.write('\t'.join(['AllWISE_ID', 'RA', 'DEC', 'all_pmRA', 'all_pmDEC',
                         'all_pmTOT', 'all_x_intrcpt', 'all_x_coeff',
                         'all_x_p_intrcpt', 'all_x_p_coeff', 'all_x_R2',
                         'all_y_intrcpt', 'all_y_coeff', 'all_y_p_intrcpt',
                         'all_y_p_coeff', 'all_y_R2']) + '\n')
    for OBJ in cent:
        out.write('\t'.join([str(i) for i in OBJ[:3] +\
                             pmAllEpoch(OBJ[0], prepCntrd(OBJ), False)]) + '\n')

with open('jg_first224_pm_preposthib', 'wb') as out:
    out.write('\t'.join(['AllWISE_ID', 'RA', 'DEC', 'pp_pmRA', 'pp_pmDEC', 
                         'pp_pmTOT']) + '\n')
    for OBJ in ppcent:
        out.write('\t'.join([str(i) for i in OBJ[:3] + pmPrePost(OBJ)])+\
                  '\n')

pdf.close()
