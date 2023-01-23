from astropy.io import fits
import numpy as np

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
    outfile = prefix+postfix.'.fits'
    hdu.writeto(outfile,overwrite=overwrite)



