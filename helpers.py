import astropy
from astropy.io import fits

import matplotlib.pyplot as plt
import numpy

import csv

def write_to_csv(data, filename):
	with open(filename, 'w', newline='') as file:
	    writer = csv.writer(file)
	    writer.writerow(data)

def write_nparr_to_fits(data, filename):
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(filename, overwrite=True)
    hdulist.close()

def compute_snr(gt, recon):
	return 20 * numpy.log10(numpy.linalg.norm(gt) / numpy.linalg.norm(gt-recon))

def compute_rmse(gt, recon):
	return numpy.sqrt(numpy.mean((gt-recon) ** 2))

def compute_maxabserr(gt, recon):
	return numpy.max(numpy.abs(gt - recon))

def compute_errstd(gt, recon):
	return numpy.std(gt - recon)

def exp_growth(x, low, high, steepness = 2):
	return low * (high / low) ** (x ** steepness)

def readFits(filename):
	dat = fits.open(filename)[0].data
	while len(dat.shape) > 2:
	    dat = dat[0]

	return dat

def plot1D(x, y, xlabel, ylabel):
	plt.plot(x, y)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.show()

def plotGDP(gt, dirty, psf, cmap):
	fig, axes = plt.subplots(1, 3)
	axes[0].set_title("True Sky")
	axes[0].imshow(gt, cmap=cmap, origin='lower')

	axes[1].set_title("Dirty Image")
	axes[1].imshow(dirty, cmap=cmap, origin='lower')

	axes[2].set_title("PSF")
	axes[2].imshow(psf, cmap=cmap, origin='lower')

	plt.show()


def read_csv(filename):
	data = []
	with open(filename, newline='') as file:
	    reader = csv.reader(file, delimiter=',')
	    for row in reader:
	    	data += row

	return data