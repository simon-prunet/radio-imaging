import astropy
from astropy.io import fits

import matplotlib.pyplot as plt
import matplotlib
import numpy

import csv

from scipy import signal

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

def compute_ssim(gt, recon):
	gt_mean = numpy.mean(gt)
	recon_mean = numpy.mean(recon)
	#gt_variance = numpy.variance(gt)
	#recon_variance = numpy.variance(recon)
	covariance_mat = numpy.cov([gt.flatten(), recon.flatten()])

	gt_variance = covariance_mat[0, 0]
	recon_variance = covariance_mat[1, 1]
	covariance = covariance_mat[0, 1]

	dynamic_range = numpy.max(gt) - numpy.min(gt) * 2

	c1 = (0.01 * dynamic_range) ** 2
	c2 = (0.03 * dynamic_range) ** 2
	c3 = c2 / 2

	l = (2 * gt_mean * recon_mean + c1) / (gt_mean ** 2 + recon_mean ** 2 + c1)
	c = (2 * numpy.sqrt(gt_variance) * numpy.sqrt(recon_variance) + c2) / (gt_variance + recon_variance + c2)
	s = (covariance + c3) / (numpy.sqrt(gt_variance) * numpy.sqrt(recon_variance) + c3)

	return (l ** 1) * (c ** 1) * (s ** 3)

	#return ((2 * gt_mean * recon_mean + c1) * (2 * covariance + c2)) / ((gt_mean ** 2 + recon_mean ** 2 + c1) * (gt_variance + recon_variance + c2))


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

def plotSNRvsSSIM(lambdas, path, snr_idx, ssim_idx, gt, cmap):
	snr_fn = path + "lambda_" + str(lambdas[snr_idx]) + ".fits"
	ssim_fn = path + "lambda_" + str(lambdas[ssim_idx]) + ".fits"

	snrfile = readFits(snr_fn)
	ssimfile = readFits(ssim_fn)

	err_snr = numpy.abs(gt - snrfile)
	err_ssim = abs(gt - ssimfile)

	fig, axes = plt.subplots(1, 2)
	axes[0].set_title("Best SNR Image $\lambda = " + "{:.2f}".format(lambdas[snr_idx]) + "$")
	im = axes[0].imshow(snrfile, cmap=cmap, origin='lower')

	axes[1].set_title("Best SSIM Image $\lambda = " + "{:.2f}".format(lambdas[ssim_idx]) + "$")
	im = axes[1].imshow(ssimfile, cmap=cmap, origin='lower')

	plt.colorbar(im, ax=axes.ravel().tolist(), shrink = 0.7) 

	plt.show()

	fig, axes = plt.subplots(1, 2)

	axes[0].set_title("SNR Absolute Error Image")
	im = axes[0].imshow(err_snr, cmap=cmap, origin='lower')

	axes[1].set_title("SSIM Absolute Error Image")
	im = axes[1].imshow(err_ssim, cmap=cmap, origin='lower')
	
	plt.colorbar(im, ax=axes.ravel().tolist(), shrink = 0.7) 

	plt.show()


def read_csv(filename):
	data = []
	with open(filename, newline='') as file:
	    reader = csv.reader(file, delimiter=',')
	    for row in reader:
	    	data += row

	return data

def computeErrorForSavedResults(gt, lambdas, path, ofilename, errorMetric):
	errors = [0] * len(lambdas)
	for i, val in enumerate(lambdas):
		curr_filename = path + "lambda_" + str(val) + ".fits"
		curr_file = readFits(curr_filename)
		errors[i] = errorMetric(gt, curr_file)

	write_to_csv(errors, path + ofilename)

def circularConv(gt, psf):
	psf_fft = numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.ifftshift(psf)))
	gt_fft = numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.ifftshift(gt)))

	return (numpy.fft.ifftshift(numpy.fft.fft2(numpy.fft.fftshift(psf_fft * gt_fft)))).real

def linearConv(gt, psf):
	return signal.fftconvolve(gt, psf, mode='same')