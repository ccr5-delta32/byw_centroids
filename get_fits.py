#!/usr/bin/env python2
# vim: set ft=python:

''' This will get all the relevant data and save a local copy. It requires a
    text file named 'objects' containing the list of candidate brown objs. '''

from astropy.io import fits
from subprocess import call, check_output
from os.path import isfile, isdir
from re import search

url = 'https://faun.rc.fas.harvard.edu/ameisner/unwise/tr_neo2/'

## Get the RA and DEC of candidate objects from the AllWISE catalog,
## requires a text file with target object AllWISE designators, 1 per line.
if not isfile('obj_coord'):
    import astroquery
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

indxtbl = findex[1].data

## Inefficient but effective function to find all rows in the indextable
## relevant to a specific position in the sky. Returns those rows as well as the
## paths to the coadds to be combined with the base url
## TODO: fins a better way to do this because it does not work for every
## possible coordinates
def findCoadds(table, RA, DEC):
    result = []
    for i in table:
        if i[3] > (RA - 0.78) and i[3] < (RA + 0.78) and i[4] > (DEC - 0.78)\
           and i[4] < (DEC + 0.78):
            result.append(i)
    paths = [ '/'.join(['e'+((3-len(str(i[2])))*'0')+str(i[2]), i[0][:3], i[0],
              'unwise-'+i[0]+'-w'+str(i[1])+'-img-u.fits']) for i in result ]
    return (result, paths)

## Save a local copy of all the w2 coadds required for the target objects to be
## able to work also when the server is down
for obj in obj_coord:
    records = findCoadds(indxtbl, float(obj[1]), float(obj[2]))
    w2 = [(records[1][i], str(j[2])) for i, j in enumerate(records[0]) if\
          j[1]==2]
    for path in w2:
        coadd = fits.open(url + path[0])
        pth = path[0].split('/')
        for i, j in enumerate(pth[:-1]):
            if not isdir('/'.join(pth[:i+1])):
                call(['mkdir', '/'.join(pth[:i+1])])
        coadd.writeto(path[0])

## Get all the epochs covered by the target objects for determining the header
## of the final table a priori
lsd = [i for i in check_output(['ls', '-d', '*/'], shell=True).split('\n')[:-1]\
       if search('e\d\d\d', i) ]
with open('epochs_included' ,'wb') as out:
    out.write('\n'.join(lsd) + '\n')

