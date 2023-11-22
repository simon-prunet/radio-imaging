import astropy
from astropy.io import fits

import sys

import numpy as np
from matplotlib import pyplot as plt

import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
from matplotlib import cm

def readFitsTo2Dnparr(filename):
	res = fits.open(filename)[0].data
	while len(res.shape) > 2:
		res = res[0]
	return res

def signaltonoise(sky, inputImage):
	error_img = sky-inputImage
	return 20*math.log10(np.linalg.norm(sky)/np.linalg.norm(error_img)), np.absolute(error_img)

sky = readFitsTo2Dnparr("sky.fits")

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True

axes = [None, None, None]

f, (axes[0], axes[1], axes[2]) = plt.subplots(1, 3)

titles = ["Full", "Low", "Multistep"]

for i in range(0, 3):
	path = sys.argv[i+1]
	deconv = readFitsTo2Dnparr(path+"/final_deconv.fits")
	sd1, err = signaltonoise(sky, deconv)

	axes[i].set_xticks([])
	axes[i].set_yticks([])

	axes[i].set_title(titles[i] + " SNR=" + str(f'{sd1:.2f}'), fontweight='bold', fontsize='19')
	im = axes[i].imshow(err, interpolation='nearest', origin='lower', cmap='turbo', norm=colors.LogNorm(vmin=0.00001, vmax=1))

	divider = make_axes_locatable(axes[i])
	cax = divider.append_axes("right", size="5%", pad=0.05)

	cb = f.colorbar(im, cax=cax, orientation='vertical')
	cb.ax.tick_params(labelsize=15)


plt.show()