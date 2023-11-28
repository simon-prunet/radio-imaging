from astropy.io import fits
from astropy import wcs
import sys

ifilename = sys.argv[1]
feather_size = int(sys.argv[2])
feather_end = int(sys.argv[3])
ofilename = sys.argv[4]

f = fits.open(ifilename)

data = f[0].data
while len(data.shape) > 2:
    data = data[0]

for y, row in enumerate(data):
    for x, pixel in enumerate(row):
        y_fixed = min(y, len(data) - y) - feather_end
        x_fixed = min(x, len(row) - x) - feather_end

        feather = 0

        if y_fixed >= 0 and x_fixed >= 0:
            tx = x_fixed / feather_size
            ctx = 3 * tx * tx - 2 * tx * tx * tx
            ty = y_fixed / feather_size
            cty = 3 * ty * ty - 2 * ty * ty * ty

            feather_x = 1 if x_fixed > feather_size else ctx
            feather_y = 1 if y_fixed > feather_size else cty

            feather = min(feather_x, feather_y)

        data[x, y] *= feather

newf = fits.PrimaryHDU()
newf.data =[[data]]
newf.header = f[0].header

hdulist = fits.HDUList([newf])
hdulist.writeto(ofilename, overwrite=True)
hdulist.close()