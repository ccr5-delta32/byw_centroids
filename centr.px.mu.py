#!/usr/bin/env python2
# vim: set ft=python

from photutils import centroid_com, centroid_1dg, centroid_2dg
from astropy import wcs
from astropy.io import fits
import astropy.visualization as vis
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from os.path import isfile, isdir
from math import ceil, cos
from re import search

## common part of the URL
url = 'https://faun.rc.fas.harvard.edu/ameisner/unwise/tr_neo2/'

## A function to find coadds for a given RA/DEC
def getCoadds(template, indxtbl, r, d, span=1.56):
    ## template is any coadd FITS object of the set to be queried
    ## This should work if all coadds evaluated have the same pixelscaling as
    ## template.
    result = []
    boxsize = 2
    box = ((1 / cos((d / 360) * (2*np.pi)) * span) * boxsize,
           span*boxsize*1.1)
    indx_subs =\
      indxtbl[1].data[np.bitwise_and(np.bitwise_and(
                                          indxtbl[1].data['RA'] > r - box[0],
                                          indxtbl[1].data['RA'] < r + box[0]),
                                        np.bitwise_and(
                                          indxtbl[1].data['DEC'] < d + box[1],
                                          indxtbl[1].data['DEC'] > d - box[1]))]
    for record in indx_subs:
        template[0].header['CRVAL1'] = record[3]
        template[0].header['CRVAL2'] = record[4]
        w = wcs.WCS(template[0].header)
        chk = w.wcs_world2pix(np.array([[r,d]]),0)
        test = (True in np.isnan(chk), np.amax(chk) >= 2048,
                np.amin(chk) < 0)
        if sum(test) != 0:
            continue
        result.append((record,
                       '/'.join(['e'+((3-len(str(record[2])))*'0')+\
                                 str(record[2]),record[0][:3], record[0],
                                 'unwise-'+record[0]+'-w'+str(record[1])+\
                                 '-img-u.fits'])))
    return result
 
## Get the RA and DEC of candidate objects from the AllWISE catalog. Requires
## the file 'objects',  containing all the target object ID's (simply 1 per line)
## to be present in the current path. The RA and DEC will be retrieved from the
## AllWISE catalog on VizieR and a new file will be written with these included
## (obj_coord). If this file is already present in the current path this will
## not be done again.
if not isfile('obj_coord'):
    from astroquery.vizier import Vizier
    objects = open('objects', 'r').read().split('\n')[:-1]
    obj_coord = []
    for drawf in objects:
        v = Vizier(columns=['AllWISE', 'RAJ2000', 'DEJ2000'],
                   column_filters={'AllWISE':drawf})
        obj_coord.append([i for i in v.query_object(drawf,
                            catalog='AllWISE')[0][0]])
    with open('obj_coord', 'wb') as out:
        for i in obj_coord:
            out.write('\t'.join([str(j) for j in i]) + '\n')

## open a list of objects with their RA and DEC (see above)
obj_coord = [i.split('\t') for i in open('obj_coord','r').\
               read().split('\n')[:-1]]

## Open a list containing, for each object, the number of pixels below, to the
## left, above, and to the right of the pixel covering the RA and DEC of the
## object retrieved from VizieR should be included in all epochs. I manually
## curated this list based on previous results to improve the convergence of the
## 2-dimensional Gaussian centroid method. By specifying a select number of
## pixels around the object noisy pixels and other objects close-by have been
## removed.
## If this file does not exist you can specify how many pixels will to
## included on each side, for all objects.
if isfile('obj_px'):
    pix = [i.split('\t') for i in open('./obj_px').read().split('\n')[1:-1]]
else:
    npx = 3 ## change this to the desired number of pixels
    pix = [ [i[0]] + [npx] * 4 for i in obj_coord ]

## Open the coadd index table either from the server or, if available, from the
## current path
if isfile('./tr_neo2_index.fits'):
    findex = fits.open('./tr_neo2_index.fits')
else:
    findex = fits.open(url+'tr_neo2_index.fits', lazy_load_hdus=False)

## A list of lists to which all the intermediate results are appended
result_hdr = ['ID', 'RA', 'DEC', 'path', 'coadd_id', 'epoch',
              'x_px_2dg', 'y_px_2dg', 'x_wcs_2dg', 'y_wcs_2dg', 'MJDMIN',
              'MJDMAX']
results = [ [] for i in result_hdr]
## another list for intermediate results from pre/post hibernation
result_hdr1 = ['ID', 'RA', 'DEC', 'coadd_id', 'prehib_x_px', 'prehib_y_px',
               'prehib_x_wcs', 'prehib_y_wcs', 'posthib_x_px', 'posthib_y_px',
               'posthib_x_wcs', 'posthib_y_wcs', 'MJDpre', 'MJDpost']
results1 = [ [] for i in result_hdr1]

## Open a pdf for writing images of the objects with the estimated centroids
pdf = PdfPages('centroids.pdf')

## open any coadd below to use while finding coads for specific positions
template = fits.open(url + 'e000/074/0741m743/unwise-0741m743-w2-img-u.fits')

## Modified Julian date of WISE hibernation start
hibmjd = 55609.8333333333

## Prepare some subsections of the data and get the MJD while we are at it
def processCoadd(w2, ipix):
    fovs = []
    MJD = []
    for i, path in enumerate(w2):
        if isfile(path[0]):
            coadd = fits.open(path[0])
        else:
            coadd = fits.open(url + path[0])
        w = wcs.WCS(coadd[0].header)
        objpx = [int(round(j, 0)) for j in \
                 w.wcs_world2pix(np.array([[float(obj[1]),
                                            float(obj[2])]]), 0)[0]]

        fov = coadd[0].data[objpx[1]-int(ipix[1]):objpx[1]+int(ipix[3])+1,
                            objpx[0]-int(ipix[2]):objpx[0]+int(ipix[4])+1]
        ## Add them to the plots
        farr[i].imshow(fov, origin='lower', interpolation='nearest',
                       cmap='Greys')
        farr[i].set_title(' | '.join([obj[0], path[0].split('/')[0],
                                     path[0].split('/')[-1]]), fontsize=4)
        fovs.append(fov)
        MJD.append( (coadd[0].header['MJDMIN'], coadd[0].header['MJDMAX']) )
    return (np.array(fovs), np.array(MJD), objpx, w)

def addArrays(arrays):
    if len(arrays) == 1:
        return arrays[0]
    first = arrays[0]
    for i in range(1, len(arrays)):
        first = first + arrays[i]
    return first

## Iterate over all objects and do stuff with them
for I, obj in enumerate(obj_coord):
    print str(I) + ': ' + obj[0]
    marker= ('+', 'x') #symbols used in the plots for centroids
    ms, mew = 10,1.5
    records = getCoadds(template, findex, float(obj[1]), float(obj[2]))
    if len(records[0]) == 0:
        f.write('\t'.join(obj) + '\n')
        continue
    w2 = [(i[1], i[0][0]) for i in records if i[0][1] == 2]
    fig, farr = plt.subplots(int(ceil(float(len(w2)) / 3 )), 3, sharex='col',
                             sharey='row' )
    farr = farr.reshape(len(farr)*3)
    [i.set_xlim([-1, sum([int(pix[I][2]), int(pix[I][4]), 1])]) for i in farr]
    [i.set_ylim([-1, sum([int(pix[I][1]), int(pix[I][3])], 1)]) for i in farr]

    # Get coadds and keep a small subsection, and get the MJD of the coadds
    fovandmjd = processCoadd(w2, pix[I])
    objpx = fovandmjd[2]
    w = fovandmjd[3]
    
    # estimate centroids using all epochs
    for J, fov in enumerate(fovandmjd[0]):
        croids = centroid_2dg(fov) # 2D Gaussian fitted to the data
        farr[J].plot(croids[0], croids[1], color='#00cccc', marker=marker[1],
                ms=ms, mew=mew)
        wcs_2dg =\
        w.wcs_pix2world(np.array([[objpx[0]-int(pix[I][2])+croids[0],
                                   objpx[1]-int(pix[I][1])+croids[1]]]), 0)
        results[0].append(obj[0])
        results[1].append(obj[1])
        results[2].append(obj[2])
        results[3].append(w2[J][1])
        results[4].append('/'.join(w2[J][0].split('/')[:-1]))
        results[5].append(w2[J][0].split('/')[0])
        results[6].append(croids[0]+objpx[0]+int(pix[I][2]))
        results[7].append(croids[1]+objpx[1]+int(pix[I][1]))
        results[8].append(wcs_2dg[0][0])
        results[9].append(wcs_2dg[0][1])
        results[10].append(fovandmjd[1][J][0])
        results[11].append(fovandmjd[1][J][1])

    plt.suptitle(obj[0] + ' | RA=' + obj[1] + ', DEC=' + obj[2])
    pdf.savefig()

    # estimate centroids in data pre and post hibernation summed, respectively
    fig, farr = plt.subplots(1,1)
    farr.set_xlim([-1, sum([int(pix[I][2]), int(pix[I][4]), 1])])
    farr.set_ylim([-1, sum([int(pix[I][1]), int(pix[I][3]), 1])])

    ## Add all pixel intensities from pre-hibernation epochs and add to the plot
    prefov = addArrays(fovandmjd[0][np.amax(fovandmjd[1], 1) <= hibmjd])
    farr.imshow(prefov, origin='lower', interpolation='nearest', cmap='Greys')
    ## estimate centroid
    pre_croids = centroid_2dg(prefov)
    ## plot centroid
    farr.plot(pre_croids[0], pre_croids[1], color='#FF8C00', marker=marker[0],
                 ms=ms, mew=mew)
    plt.suptitle(obj[0] + ' | added pixel intensities pre-hibernation')
    pdf.savefig()

    ## Add all pixel intensities from post-hibernation epochs and add to the plot
    fig, farr = plt.subplots(1,1)
    farr.set_xlim([-1, sum([int(pix[I][2]), int(pix[I][4]), 1])])
    farr.set_ylim([-1, sum([int(pix[I][1]), int(pix[I][3]), 1])])

    postfov = addArrays(fovandmjd[0][np.amax(fovandmjd[1], 1) > hibmjd])
    farr.imshow(postfov, origin='lower', interpolation='nearest', cmap='Greys')
    ## estimate centroid
    post_croids = centroid_2dg(postfov)
    ## plot centroid
    farr.plot(post_croids[0], post_croids[1], color='#FF8C00', marker=marker[0],
                 ms=ms, mew=mew)

    results1[0].append(obj[0])
    results1[1].append(obj[1])
    results1[2].append(obj[2])
    results1[3].append(w2[I][1])
    results1[4].append(pre_croid[0]+objpx[0]+int(pix[I][2]))
    results1[5].append(pre_croid[1]+objpx[1]+int(pix[I][1]))
    prewcs = w.wcs_pix2world(np.array([[objpx[0] - int(pix[I][2]) +\
                                        pre_croids[0],
                                        objpx[1] - int(pix[I][1]) +\
                                        pre_croids[1]]]), 0)
    results1[6].append(prewcs[0][0])
    results1[7].append(prewcs[0][1])
    results1[8].append(post_croid[0]+objpx[0]+int(pix[I][2]))
    results1[9].append(post_croid[1]+objpx[1]+int(pix[I][1]))
    postwcs = w.wcs_pix2world(np.array([[objpx[0] - int(pix[I][2]) +\
                                        post_croids[0],
                                        objpx[1] - int(pix[I][1]) +\
                                        post_croids[1]]]), 0)
    results1[10].append(postwcs[0][0])
    results1[11].append(postwcs[0][1])
    results1[12].append(np.mean(fovandmjd[1][np.amax(fovandmjd[1], 1)<=hibmjd]))
    results1[13].append(np.mean(fovandmjd[1][np.amax(fovandmjd[1], 1)>hibmjd]))

    plt.suptitle(obj[0] + ' | added pixel intensities post-hibernation')
    pdf.savefig()
    plt.close('all')
pdf.close()

## write temporary results to file
#with open('arb_px_run1_tmp', 'wb') as out:
#    for i in range(len(results[0])):
#        out.write('\t'.join([str(j[i]) for j in results]) + '\n')

## Produce final table now because beforehand not all epochs included are known
unipoch = sorted(set(results[5])) # all unique epochs included in the results
## LOL!
import itertools
final_hdr = ['ID', 'RA', 'DEC', 'coadd_id'] +\
            list(itertools.chain.from_iterable(
                 list(itertools.chain.from_iterable(
                      zip([[i + ['_px_x', '_px_y', '_wcs_x' ,'_wcs_y',
                                 '_MJDMIN', '_MJDMAX'][j] for j in
                           range(6)] for i in unipoch])))))
## An array for the sorted-out data
final = np.empty([len(final_hdr), len(set(results[3]))],
                    dtype=[ (j, 'O') if i in [0, 3] else (j, 'float64') for i, j
                           in enumerate(final_hdr) ])
final[:] = np.NaN

## Process the intermediate results and populate the final table for all epochs
for I, coadd in enumerate(set(results[3])):
    res = [np.nan for i in range(len(final_hdr))]
    dresult = [[ k[i] for k in results ] for i in range(len(results[0])) if
               results[3][i] == coadd ]
    res[0] = dresult[0][0]
    res[1:3] = [float(i) for i in dresult[0][1:3]]
    res[3] = dresult[0][3]
    for d in dresult:
        hdr_indx = [i for i in range(len(final_hdr)) if\
                    search(d[5],final_hdr[i])]
        res[hdr_indx[0]:hdr_indx[-1]+1] = d[-6:]
    final[0][I] = tuple(res)

## Write final table for all epochs
with open('centroids.2dg', 'wb') as out:
    out.write('\t'.join(final_hdr) + '\n')
    for row in np.sort(final,order='ID')[0]:
        out.write('\t'.join([str(i) for i in row]) + '\n')

## write results for pre/post hibernation
with open('preposhib_centroids.2dg', 'wb') as out:
    out.write('\t')
    for i in range(len(results1[0])):
        out.write('\t'.join([str(j[i]) for j in results1]) + '\n')





