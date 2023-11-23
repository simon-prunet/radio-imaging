from astropy.io import fits
from astropy import wcs
import sys

ifilename = sys.argv[1]
center_x = int(sys.argv[2])
center_y = int(sys.argv[3])
halfdim = int(sys.argv[4])
hfilename = sys.argv[5]
ofilename = sys.argv[6]

f = fits.open(ifilename)
h = fits.open(hfilename)
w = wcs.WCS(f[0].header)

newf = fits.PrimaryHDU()
newf.data = [[f[0].data[center_x - halfdim:center_x + halfdim, center_y - halfdim:center_y + halfdim]]]
newf.header = h[0].header

hdulist = fits.HDUList([newf])
hdulist.writeto(ofilename, overwrite=True)
hdulist.close()