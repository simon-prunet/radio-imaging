#!/usr/bin/env python

"""\
Code used to handle RASCIL images
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

from astropy.io import fits
import numpy as np
import os.path

from ska_sdp_datamodels.image import Image
from rascil.processing_components.image.operations import import_image_from_fits
from rascil.processing_components.simulation.testing_support import replicate_image
from ska_sdp_func_python.imaging import create_image_from_visibility
from ska_sdp_datamodels.visibility.vis_model import Visibility

def to_rascil_format(fitsfile,postfix='_ext',overwrite=True):
    """
    to_rascil_format Takes a regular 2D image and completes its header to make it a 4D rascil FITS image

    :fitsfile: filename of 2d image
    :postfix: extension of new file
    :overwrite: overwrite existing file
    """

    hdulist = fits.open(fitsfile)
    header = hdulist[0].header
    data = hdulist[0].data
    if (header['NAXIS']==2):
        print('This is a regular 2D image. Will be transformed into a RASCIL 4D image')
    elif(header['NAXIS']==4):
        print('This already looks like a RASCIL 4D image ? Exiting.')
        return
    else:
        print('NAXIS value is neither 2 nor 4. Unknown format. Exiting')
        return

    # Add 2 axes: 1 for Stokes, 1 for frequencies
    header['NAXIS']=4
    header['WCSAXES']=4
    
    header['NAXIS3']=1
    header['NAXIS4']=1

    header['CRPIX3']=1.0
    header['CRPIX4']=1.0
    header['CRVAL3']=1.0
    header['CRVAL4']=100000000.0 # Hz

    header['CDELT3']=1.0
    header['CDELT4']=100000.0 # Hz

    header['CTYPE3']='STOKES'
    header['CTYPE4']='FREQ'

    header['CUNIT4']='Hz'

    header['RADESYS']='ICRS' # Needed ?

    header['CTYPE1']='RA---SIN'
    header['CTYPE2']='DEC--SIN'

    # Reshape data in 4D. Beware that quick axes come first in FITS, and last in python (C-order)
    ny,nx = data.shape
    data = data.reshape([1,1,ny,nx])

    hdu = fits.PrimaryHDU(data,header)
    prefix = fitsfile.split('.fits')[0]
    outfile = prefix+postfix+'.fits'
    hdu.writeto(outfile,overwrite=overwrite)
    return


# Almost cut and paste from create_test_image (rascil library)
def create_image_from_fits(
    fitsfile,
    cellsize=None,
    frequency=None,
    channel_bandwidth=None,
    phasecentre=None,
    polarisation_frame=None,
) -> Image:
    """create_image_from_fits creates a RASCIL image from a given fits file

    :param fitsfile: Input image FITS file
    :param cellsize: angular resolution of pixel pixel size in dg
    :param frequency: Frequency (array) in Hz
    :param channel_bandwidth: Channel bandwidth (array) in Hz
    :param phasecentre: Phase centre of image (SkyCoord)
    :param polarisation_frame: Polarisation frame
    :return: Image
    """
    

    if not os.path.exists(fitsfile):
        print ("Input FITS file %s does not exist, exiting."%fitsfile)
        return None

    im = import_image_from_fits(fitsfile)

    if frequency is None:
        frequency = [1e8]
    if polarisation_frame is None:
        polarisation_frame = im.image_acc.polarisation_frame
    im = replicate_image(im, frequency=frequency, polarisation_frame=polarisation_frame)

    wcs = im.image_acc.wcs.deepcopy()

    if cellsize is not None:
        wcs.wcs.cdelt[0] = -180.0 * cellsize / np.pi
        wcs.wcs.cdelt[1] = +180.0 * cellsize / np.pi
    if frequency is not None:
        wcs.wcs.crval[3] = frequency[0]
    if channel_bandwidth is not None:
        wcs.wcs.cdelt[3] = channel_bandwidth[0]
    else:
        if len(frequency) > 1:
            wcs.wcs.cdelt[3] = frequency[1] - frequency[0]
        else:
            wcs.wcs.cdelt[3] = 0.001 * frequency[0]
    wcs.wcs.radesys = "ICRS"
    wcs.wcs.equinox = 2000.00

    if phasecentre is None:
        phasecentre = im.image_acc.phasecentre
    else:
        wcs.wcs.crval[0] = phasecentre.ra.deg
        wcs.wcs.crval[1] = phasecentre.dec.deg
        # WCS is 1 relative
        wcs.wcs.crpix[0] = im["pixels"].data.shape[3] // 2 + 1
        wcs.wcs.crpix[1] = im["pixels"].data.shape[2] // 2 + 1

    return Image.constructor(
        im["pixels"].data, wcs=wcs, polarisation_frame=polarisation_frame
    )


def create_empty_image(vis, npixel, cellsize):
    """create_empty_image creates an empty RASCIL image from a given set of visibilities and parameters

    :param vis: visibilities
    :param npixel: number of pixels on a dimension, assumed to be square
    :param cellsize: angular resolution of pixel
    :return: Empty rascil Image
    """
    return create_image_from_visibility(vis, npixel=npixel, cellsize=cellsize, polarisation_frame=vis.visibility_acc.polarisation_frame)


def add_to_image(image, nparr):
    """add_to_image adds a numpy array to some given rascil image

    :param image: RASCIL image
    :param nparr: array to add
    :return: rascil Image with added numpy array values
    """
    image.pixels.data[0,0,:,:] += nparr

    return image


