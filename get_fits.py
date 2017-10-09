#!/usr/bin/env python2
# vim: set ft=python:

''' This will get all the relevant data and save a local copy. It requires a
    text file named 'objects' containing the list of candidate brown objs. '''

from astropy.io import fits
from subprocess import call, check_output
from os.path import isfile, isdir
from re import search
from astropy import wcs
import numpy
from math import cos


url = 'https://faun.rc.fas.harvard.edu/ameisner/unwise/tr_neo2/'

## Get the RA and DEC of candidate objects from the AllWISE catalog,
## requires a text file with target object AllWISE designators, 1 per line.
if not isfile('obj_coord'):
    import astroquery
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
else:
    obj_coord = [i.split('\t') for i in open('obj_coord','r').\
                   read().split('\n')[:-1]]

## Get the indextable for all the coadds and save a local copy if this has not
## been done yet
if not isfile('tr_neo2_index.fits'):
    findex = fits.open(url+'tr_neo2_index.fits', lazy_load_hdus=False)
    findex.writeto('tr_neo2_index.fits')
else:
    findex = fits.open('tr_neo2_index.fits', lazy_load_hdus=False)

## A new and improved method to find coadds for a given RA/DEC
## inspired by the http://unwise.me code in github
def getCoadds(template, indxtbl, r, d, span=1.56):
    ## template is any coadd FITS object of the set to be queried
    ## This should work if all coadds evaluated have the same pixelscaling as
    ## template.
    result = []
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

## Save a local copy of all the w2 coadds required for the target objects to be
## able to work also when the server is down
## We need a FITS object of any 1 coadd to use when looking for specific ones
if isfile('e000/074/0741m743/unwise-0741m743-w2-img-u.fits'):
    template = fits.open('e000/074/0741m743/unwise-0741m743-w2-img-u.fits')
else:
    template = fits.open(url +
                         'e000/074/0741m743/unwise-0741m743-w2-img-u.fits')

for obj in obj_coord:
    records = getCoadds(template, findex, float(obj[1]), float(obj[2]))
    w2 = [(i[1], i[0][0]) for i in records if i[0][1] == 2]
    for path in w2:
        if isfile(path[0]):
            continue
        coadd = fits.open(url + path[0])
        pth = path[0].split('/')
        for i, j in enumerate(pth[:-1]):
            if not isdir('/'.join(pth[:i+1])):
                call(['mkdir', '/'.join(pth[:i+1])])
        coadd.writeto(path[0])
