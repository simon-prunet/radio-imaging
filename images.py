from astropy.io import fits
import numpy as np

from ska_sdp_datamodels.image import Image

def to_rascil_format(fitsfile,postfix='_ext',overwrite=True):
    '''
    Takes a regular 2D image and complete its header
    to make it a 4D rascil FITS image
    '''

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
    """Create a useful test image

    This is the test image M31 widely used in ALMA and other simulations. It is actually part of an Halpha region in
    M31.

    :param fitsfile: Input image FITS file
    :param cellsize: pixel size in dg
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
        wcs.wcs.cdelt[0] = -180.0 * cellsize / numpy.pi
        wcs.wcs.cdelt[1] = +180.0 * cellsize / numpy.pi
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

    return create_image_from_array(
        im["pixels"].data, wcs=wcs, polarisation_frame=polarisation_frame
    )



