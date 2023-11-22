import numpy
import sys
import os

import astropy
from astropy.io import fits

import csv

def write_to_csv(data, filename):
	with open(filename, 'w', newline='') as file:
	    writer = csv.writer(file)
	    writer.writerows(map(lambda x: [x], data))

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

gt_filename = sys.argv[1]
dirty_filename = sys.argv[2]
psf_filename = sys.argv[3]

lambda_low = float(sys.argv[4])
lambda_high = float(sys.argv[5])
num_samples = int(sys.argv[6])

output_filename = sys.argv[7]

testname = sys.argv[8]

lambdas = numpy.arange(0, 1, 1 / float(num_samples))
for i, v in enumerate(lambdas):
	lambdas[i] = exp_growth(v, lambda_low, lambda_high)

snrs = [0] * len(lambdas)
rmses = [0] * len(lambdas)
abserr = [0] * len(lambdas)
errstd = [0] * len(lambdas)

dirty = fits.open(dirty_filename)[0].data
if len(dirty.shape) > 2:
    while len(dirty.shape) > 2:
        dirty = dirty[0]

    write_nparr_to_fits(dirty, dirty_filename)

max_val = numpy.max(numpy.abs(dirty))

gt = fits.open(gt_filename)[0].data
while len(gt.shape) > 2:
    gt = gt[0]

for i, curr_lambda in enumerate(lambdas):
	curr_output_name = testname + "lambda_" + str(curr_lambda) + ".fits"
	#normalization of lambda, otherwise ideal value of lambda becomes more dependent on image amplitudes
	curr_lambda /= max_val if max_val > 1e-10 else 9999999999
	os.system("julia make_fullres.jl " + str(curr_lambda) + " " + psf_filename + " " + dirty_filename + " " + curr_output_name)
	recon = fits.open(curr_output_name)[0].data
	snrs[i] = compute_snr(gt, recon)
	rmses[i] = compute_rmse(gt, recon)
	abserr[i] = compute_maxabserr(gt, recon)
	errstd[i] = compute_errstd(gt, recon)

write_to_csv([snrs, rmses, errstd, abserr], output_filename)