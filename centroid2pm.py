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





cen = np.array(cent[0][4:], dtype='float32').reshape((len(cent[0])-4)/6,6)
cen = cen[np.sum(np.isnan(cen), 1) < 1]

# design matrix
Ax = np.matrix([np.ones(np.shape(cen)[0]), np.mean(cen[...,-2:], 1)]).T
Ay = np.matrix([np.ones(np.shape(cen)[0]), np.mean(cen[...,-2:], 1)]).T

b_x = cen[...,0]
b_y = cen[...,1]

res_x = sm.OLS(b_x, Ax).fit()
prstd_x, ci_l_x, ci_u_x = wls_prediction_std(res_x)
res_y = sm.OLS(b_y, Ay).fit()
prstd_y, ci_l_y, ci_u_y = wls_prediction_std(res_y)

fig, pan = plt.subplots(1,2)
pan[0].scatter(np.mean(cen[...,-2:], 1), cen[...,0], color='deepskyblue')
pan[0].plot(np.mean(cen[...,-2:], 1), res_x.fittedvalues, color='darkorange')
pan[1].scatter(np.mean(cen[...,-2:], 1), cen[...,1], color='deepskyblue')
pan[1].plot(np.mean(cen[...,-2:], 1), res_y.fittedvalues, color='darkorange')
fig.tight_layout()
fig.show()

pm = (-(res_x.params[1] * 2.75) * 365, (res_y.params[1] * 2.75) * 365, 
     np.sqrt(((res_x.params[1] * 2.75) * 365)**2 +\
             ((res_y.params[1] * 2.75) * 365)**2))
################################################################################

################################################################################
## same but different, more stats                                             ##
################################################################################

obj = ['J110125.95+540052.8',  '165.3581626', '54.0146782']

records = getCoadds(template, findex, float(obj[1]), float(obj[2]))

w2 = [(i[1], i[0][0]) for i in records if i[0][1] == 2]

Ab = np.vstack([coadd2Mat(i[0], 4) for i in w2])
A = Ab[...,:3]
A[...,2] = A[...,2]/10000
A = sm.add_constant(A)
b = Ab[...,4]

res = sm.OLS(b, A).fit()
res.summary()

################################################################################
## Estimate PM for WISEA J110125.95+540052.8 using centroid estimates         ##
################################################################################

cent = [i.split('\t') for i in open('centroids.2dg').read().split('\n')[:-1]]
cenh = cent[0]
Cen = [float(i) if not search('[a-zA-Z]', i) else i for i in cent[1]]
cen = np.array(Cen[4:]).reshape(6,6)

pxfit = np.polyfit(cen[...,0], cen[...,1], 1, full=True)
wcsfit = np.polyfit(cen[...,2], cen[...,3], 1, full=True)

y1px = [pxfit[0][1] + pxfit[0][0] * i for i in cen[...,0]]

fig, pan = plt.subplots(1,1)
pan.scatter(cen[...,0], cen[...,1])
pan.plot(cen[...,0], y1px)
fig.show()


avgmjd = np.mean(cen[...,4:5], 1)

pdf = PdfPages('regression')
#estimate y per unit mjd
pxfity = np.polyfit(avgmjd, cen[...,1], 1, full=True)
y1px = [pxfity[0][1] + pxfity[0][0] * i for i in avgmjd]
fig, pan = plt.subplots(1,2)
pan[0].scatter(avgmjd, cen[...,1])
pan[0].plot(avgmjd, y1px)
pan[0].tick_params(axis='both', which='major', labelsize=10)
pan[0].set_title('y coordinates')
pan[0].set_ylabel('Pixels')
pan[0].set_xlabel('MJD')
pmDEC = (pxfity[0][0] * 2.75) * 365 #px=3: 0.10398721796262224
#estimate x per unit mjd
pxfitx = np.polyfit(avgmjd, cen[...,0], 1, full=True)
y2px = [pxfitx[0][1] + pxfitx[0][0] * i for i in avgmjd]
pan[1].scatter(avgmjd, cen[...,0])
pan[1].plot(avgmjd, y2px)
pan[1].tick_params(axis='both', which='major', labelsize=10)
pan[1].set_title('x coordinates')
pan[1].set_xlabel('MJD')
fig.tight_layout()
pdf.savefig()
pdf.close()
fig.show()
pmRA = (pxfitx[0][0] * 2.75) * 365 #px=3: 0.90799224214238516

PM = np.sqrt(pxfity[0][0]**2 + pxfitx[0][0]**2) * 2.75 * 365 #px=3: 0.91392737856481898







