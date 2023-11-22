import astropy
from astropy.io import fits

import sys

import numpy as np
from matplotlib import pyplot as plt

import math
from mpl_toolkits.axes_grid1 import make_axes_locatable

def readFitsTo2Dnparr(filename):
	res = fits.open(filename)[0].data
	while len(res.shape) > 2:
		res = res[0]
	return res

def signaltonoise(sky, inputImage):
	error_img = sky-inputImage
	return 20*math.log10(np.linalg.norm(sky)/np.linalg.norm(error_img)), np.absolute(error_img)

path = sys.argv[1]

sky = readFitsTo2Dnparr("sky.fits")
dirty = readFitsTo2Dnparr("dirty.fits")
first_deconv = readFitsTo2Dnparr(path+"/deconv_iteration_0_channel_0.fits")
last_deconv = readFitsTo2Dnparr(path+"/final_deconv.fits")
last_residual = readFitsTo2Dnparr(path+"/dirty_iteration_5_channel_0.fits")

sd1, err_first = signaltonoise(sky, first_deconv)
sd2, err_last = signaltonoise(sky, last_deconv)
sddirty = signaltonoise(sky, dirty)

print(sddirty)

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True

f, (ax1, ax2, ax3) = plt.subplots(1, 3)
ax1.set_title('Deconv, major cycle 1, SNR=' + str(f'{sd1:.2f}'), fontweight='bold', fontsize='19')
ax2.set_title('Deconv, major cycle 5, SNR=' + str(f'{sd2:.2f}'), fontweight='bold', fontsize='19')
ax3.set_title('Long baseline residual, major cycle 5', fontweight='bold', fontsize='19')

ax1.set_xticks([])
ax1.set_yticks([])

ax2.set_xticks([])
ax2.set_yticks([])

ax3.set_xticks([])
ax3.set_yticks([])

im1 = ax1.imshow(err_first, interpolation='nearest', cmap='turbo', origin='lower')
im2 = ax2.imshow(err_last, interpolation='nearest', cmap='turbo', origin='lower')
im3 = ax3.imshow(last_residual, interpolation='nearest', cmap='turbo', origin='lower')

divider1 = make_axes_locatable(ax1)
divider2 = make_axes_locatable(ax2)
divider3 = make_axes_locatable(ax3)
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cax2 = divider2.append_axes("right", size="5%", pad=0.05)
cax3 = divider3.append_axes("right", size="5%", pad=0.05)

cb1 = f.colorbar(im1, cax=cax1, orientation='vertical')
cb2 = f.colorbar(im2, cax=cax2, orientation='vertical')
cb3 = f.colorbar(im3, cax=cax3, orientation='vertical')

cb1.ax.tick_params(labelsize=15)
cb2.ax.tick_params(labelsize=15)
cb3.ax.tick_params(labelsize=15)

plt.show()