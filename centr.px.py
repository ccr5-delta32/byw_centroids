#!/usr/bin/env python2
# vim: set ft=python

from photutils import centroid_com, centroid_1dg, centroid_2dg
from astropy import wcs
from astropy.io import fits
import astropy.visualization as vis
import matplotlib.pyplot as plt
import numpy
from matplotlib.backends.backend_pdf import PdfPages
from os.path import isfile, isdir
from math import ceil, cos
from re import search

## A new and improved method to find coadds for a given RA/DEC
## inspired by the http://unwise.me code in github
def getCoadds(template, indxtbl, r, d, span=1.56):
    ## template is any coadd FITS object of the set to be queried
    ## This should work if all coadds evaluated have the same pixelscaling as
    ## template.
    result = []
    # TODO wrapping at 0, 360, 90, -90
    boxsize = 2 
    box = ((1 / cos((d / 360) * (2*numpy.pi)) * span) * boxsize,
           span*boxsize*1.1)
    indx_subs =\
      indxtbl[1].data[numpy.bitwise_and(numpy.bitwise_and(
                                          indxtbl[1].data['RA'] > r - box[0],
                                          indxtbl[1].data['RA'] < r + box[0]),
                                        numpy.bitwise_and(
                                          indxtbl[1].data['DEC'] < d + box[1],
                                          indxtbl[1].data['DEC'] > d - box[1]))]
    for record in indx_subs:
        template[0].header['CRVAL1'] = record[3]
        template[0].header['CRVAL2'] = record[4]
        w = wcs.WCS(template[0].header)
        chk = w.wcs_world2pix(numpy.array([[r,d]]),0)
        test = (True in numpy.isnan(chk), numpy.amax(chk) >= 2048,
                numpy.amin(chk) < 0)
        if sum(test) != 0:
            continue
        result.append((record,
                       '/'.join(['e'+((3-len(str(record[2])))*'0')+\
                                 str(record[2]),record[0][:3], record[0],
                                 'unwise-'+record[0]+'-w'+str(record[1])+\
                                 '-img-u.fits'])))
    return result
 
## Get the RA and DEC of candidate objects from the ALlWISE catalog. Requires
## the file 'objects',  containing all the target objects (simply 1 per line)
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

## Open a list containing, for each object, the number of pixels below, to the
## left, above, and to the right of the pixel covering the RA and DEC of the
## object retrieved from VizieR should be included in all epochs. I manually
## curated this list based on previous results to improve the convergence of the
## 2-dimensional Gaussian centroid method. By specifying a select number of
## pixels around the object noisy pixels and other objects close-by have been
## removed.
## If this file does not exist you can specify how many pixels will to
## included on each side, for all objects (not tested yet!).
if isfile('obj_px'):
    pix = [i.split('\t') for i in open('./obj_px').read().split('\n')[1:-1]]
else:
    npx = 7 ## change this to the desired number of pixels
    pix = [ [i] + [npx] * 4 for i in objects ]

## A list of lists to which all the intermediate results are appended
result_hdr = ['ID', 'RA', 'DEC', 'path', 'coadd_id', 'epoch',
              'x_px_com', 'y_px_com', 'x_wcs_com', 'y_wcs_com', 'x_px_1dg',
              'y_px_1dg', 'x_wcs_1dg', 'y_wcs_1dg', 'x_px_2dg', 'y_px_2dg',
              'x_wcs_2dg', 'y_wcs_2dg', 'MJDMIN','MJDMAX']
results = [ [] for i in result_hdr]

## Open a pdf for writing images of the objects with the estimated centroids
pdf = PdfPages('centroids.pdf')

## open a list of objects with their RA and DEC (see above)
obj_coord = [i.split('\t') for i in open('obj_coord','r').\
               read().split('\n')[:-1]]

## common part of the URL
url = 'https://faun.rc.fas.harvard.edu/ameisner/unwise/tr_neo2/'

## Open the coadd index table either from the server or, if available, from the
## current path
if isfile('./tr_neo2_index.fits'):
    findex = fits.open('./tr_neo2_index.fits')
else:
    findex = fits.open(url+'tr_neo2_index.fits', lazy_load_hdus=False)

## In case no coadd can be found for a given object it will be saved in this
f = open('check_these_because_no_records', 'wb')

## open any coadd below to use while finding coads for specific positions
template = fits.open('./e000/074/0741m743/unwise-0741m743-w2-img-u.fits')

## Iterate over all objects and do stuff with them
for I, obj in enumerate(obj_coord):
    print str(I) + ': ' + obj[0]
    atrimb = [1, -4000, 4000] # like wiseview: ['linear',clip_low,'trim_bright']
    marker= '+' # symbol used in the plots to indicate the computed centroid
    ms, mew = 20, 2
    records = getCoadds(template, findex, float(obj[1]), float(obj[2]))
    if len(records[0]) == 0:
        f.write('\t'.join(obj) + '\n')
        continue
    w2 = [(j[1], str(j[0][2])) for i, j in enumerate(records) if\
          j[0][1]==2]
    fig, farr = plt.subplots(int(ceil(float(len(w2)) / 3 )), 3, sharex='col',
                             sharey='row' )
    farr = farr.reshape(len(farr)*3)
    [i.set_xlim([-1, sum([int(pix[I][2]), int(pix[I][4]), 1])]) for i in farr]
    [i.set_ylim([0, sum([int(pix[I][1]), int(pix[I][3])])]) for i in farr]
    for i, path in enumerate(w2):
        if isfile(path[0]):
            coadd = fits.open(path[0])
        else:
            coadd = fits.open(url + path[0])
        w = wcs.WCS(coadd[0].header)
        objpx = [int(round(j, 0)) for j in \
                 w.wcs_world2pix(numpy.array([[float(obj[1]),
                                               float(obj[2])]]), 0)[0]]
        fov = coadd[0].data[objpx[1]-int(pix[I][1]):objpx[1]+int(pix[I][3])+1,
                            objpx[0]-int(pix[I][2]):objpx[0]+int(pix[I][4])+1]
        fov = numpy.clip(fov, atrimb[1], atrimb[2])
        farr[i].imshow(fov, origin='lower', interpolation='nearest',
                       cmap='Greys', norm=vis.mpl_normalize.\
                       ImageNormalize(stretch = vis.AsinhStretch(atrimb[0])))
        farr[i].set_title(' | '.join([obj[0], path[0].split('/')[0],
                                     path[0].split('/')[-1]]), fontsize=4)
        croids = centroid_com(fov) # center of mass from 2D image moments
        croids = numpy.append(croids, centroid_1dg(fov)) # 1D Gaussians fit to the marginal x and y
        croids = numpy.append(croids, centroid_2dg(fov)) # 2D Gaussian fitted to the 2D distribution
        farr[i].plot(croids[0], croids[1], color='#00aa00', marker=marker,
                ms=ms, mew=mew)
        farr[i].plot(croids[2], croids[3], color='#0000ff', marker=marker,
                ms=ms, mew=mew)
        farr[i].plot(croids[4], croids[5], color='#00cccc', marker=marker,
                ms=ms, mew=mew)
        wcs_com =\
        w.wcs_pix2world(numpy.array([[objpx[0]-int(pix[I][2])+croids[0],
                                      objpx[1]-int(pix[I][1])+croids[1]]]), 0)
        wcs_1dg =\
        w.wcs_pix2world(numpy.array([[objpx[0]-int(pix[I][2])+croids[2],
                                      objpx[1]-int(pix[I][1])+croids[3]]]), 0)
        wcs_2dg =\
        w.wcs_pix2world(numpy.array([[objpx[0]-int(pix[I][2])+croids[4],
                                      objpx[1]-int(pix[I][1])+croids[5]]]), 0)
        results[0].append(obj[0])
        results[1].append(obj[1])
        results[2].append(obj[2])
        results[3].append(records[0][0][0])
        results[4].append('/'.join(path[0].split('/')[:-1]))
        results[5].append(path[0].split('/')[0])
        results[6].append(croids[0]+objpx[0]+int(pix[I][2]))
        results[7].append(croids[1]+objpx[1]+int(pix[I][1]))
        results[8].append(wcs_com[0][0])
        results[9].append(wcs_com[0][1])
        results[10].append(croids[2]+objpx[0]+int(pix[I][2]))
        results[11].append(croids[3]+objpx[1]+int(pix[I][1]))
        results[12].append(wcs_1dg[0][0])
        results[13].append(wcs_1dg[0][1])
        results[14].append(croids[4]+objpx[0]+int(pix[I][2]))
        results[15].append(croids[5]+objpx[1]+int(pix[I][1]))
        results[16].append(wcs_2dg[0][0])
        results[17].append(wcs_2dg[0][1])
        results[18].append(coadd[0].header[-1])
        results[19].append(coadd[0].header[-2])
    plt.suptitle(obj[0] + ' | RA=' + obj[1] + ', DEC=' + obj[2])
    pdf.savefig()
    plt.close('all')
pdf.close()
f.close()

## write temporary results to file
with open('arb_px_run1_tmp', 'wb') as out:
    for i in range(len(results[0])):
        out.write('\t'.join([str(j[i]) for j in results]) + '\n')

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
final = numpy.empty([len(final_hdr), len(obj_coord)],
                    dtype=[ (j, 'O') if i in [0, 3] else (j, 'float64') for i, j
                           in enumerate(final_hdr) ])
final[:] = numpy.NaN

## Process the intermediate results and populate the final table
## Inspection of the pdf shows that the centroid_2dg is most robust so only this
## is included
for I, obj in enumerate(obj_coord):
    res = [numpy.nan for i in range(len(final_hdr))]
    res[0] = obj[0]
    res[1:3] = [float(i) for i in obj[1:]]
    dresult = [[ k[i] for k in results ] for i in range(len(results[0])) if
               results[0][i] == obj[0] ]
    if len(dresult) == 0: 
        res[3] = None
        final[0][I] = tuple(res)
        continue
    res[3] = dresult[0][3]
    for d in dresult:
        hdr_indx = [i for i in range(len(final_hdr)) if\
                    search(d[5],final_hdr[i])]
        res[hdr_indx[0]:hdr_indx[-1]+1] = d[-6:]
    final[0][I] = tuple(res)

## Write final table
with open('centroids.2dg', 'wb') as out:
    out.write('\t'.join(final_hdr) + '\n')
    for row in final[0]:
        out.write('\t'.join([str(i) for i in row]) + '\n')
