## Centroids of select objects for WISE-only position measurements.

## what it does
centr.px.py estimates the subpixel centroids of objects of interest from the
AllWISE catalog using the
[photutils](https://photutils.readthedocs.io/en/stable/#) package for python.
It provides the 3 functions:
1. centroid_com, which estimates the center of mass from 2d image moments
2. centroid_1dg, which estimates the centroids from 1d Gaussians fit to the
   marginalx any distributions.
3. centroid_2dg, which estimates the centroid from a 2d Gaussian fitted to the
   2d distribution.

The script applies all of these and the results are plotted over the actual
image and saved in a pdf ('centroids.pdf'). The markers are colored green, blue,
cyan to indicate the respective estimation method.

In get_fits.py there is some code amongst others that will save a local copy of
all the relevant coadds, which can save time.

## how it does it
In principle only a list of candidate objects is needed that for now needs be
called 'objects'. RA and DEC of these objects is retrieved from the AllWISE
catalog in VizieR, and saved to a new file called 'obj_coord' (unless this
has already been done).

At first I selected a fixed but small number of pixels around each object and
estimated the subpixel centroid in that. However, nearby objects and noisy
pixels turned out to have large effects on the estimates. Eventually I specified
for each object how many pixels to include below, to the left, above, and to the
right of the pixel covering the object coordinates after visual inspection of
the images in the resulting pdf file (optional file 'obj_px'). Although this
is not particularly elegant it is effective and seems to greatly improve the
accuracy of the estimations. Nevertheless, some images are so noisy, or have
such high background intensities that currently the subpixel centroid estimate
is unreliable; this can be specific to a single epoch only. The images in the
pdf should help to identify them.

The final result is a table (centroids.2dg) with in the first 4 columns:
object_id, object_RA, object_DEC, and coadd_id, followed by the 6 columns per
epoch: px_x, px_y (subpixel centroids in pixels), wcs_x, wcs_y (subpixel
centroids in world coordinates), MJDMAX, MJDMIN (MJD values from the coadd
FITS header). Where coadd_id is None I did not find a coadd containing the
object coordinates yet using my method. Only the centroid_2dg results are
included bcause these were most robust as became apparent from eyeballing the
images in the pdf.
