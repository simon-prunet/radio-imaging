#!/usr/bin/env python

"""\
notebook specific code used for regularization related tests
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import numpy
import matplotlib

from scipy import signal
import matplotlib.pyplot as plt

from radioimaging.util import util
from radioimaging.evaluation import evaluation
import split
from radioimaging.deconvolution import deconvolve

def computeErrorForSavedResults(gt, lambdas, path, ofilename, errorMetric):
    errors = [0] * len(lambdas)
    for i, val in enumerate(lambdas):
        curr_filename = path + "lambda_" + str(val) + ".fits"
        curr_file = util.fromfits(curr_filename)
        errors[i] = errorMetric(gt, curr_file)

    write_to_csv(errors, path + ofilename)

def plotSNRvsSSIM(lambdas, path, snr_idx, gt, cmap, same_scale = False):
    snr_fn = path + "lambda_" + str(lambdas[snr_idx]) + ".fits"
    snrfile = util.fromfits(snr_fn)
    err_snr = numpy.abs(gt - snrfile)
    snr_title = "Best SNR Image $\lambda = " + "{:.2f}".format(lambdas[snr_idx]) + "$"

    util.plotNImages([snrfile, err_snr], [snr_title, "SNR absolute error"], cmap, same_scale = same_scale)


def lambdatests_allvis(lambdas, path, dirty, psf, gt, wavelet_type_idx, niter, runtests, lp_filter=None, stats_type = "psnr"):
    if lp_filter is not None:
        gt = signal.fftconvolve(gt, lp_filter, mode='same')

    if runtests:
        snrs = [0] * len(lambdas)

        tmp_dirty_fn = "tmp_dirty.fits"
        tmp_psf_fn = "tmp_psf.fits"

        util.tofits(dirty, tmp_dirty_fn)
        util.tofits(psf, tmp_psf_fn)

        for i, curr_lambda in enumerate(lambdas):
            curr_output_name = path + "lambda_" + str(curr_lambda) + ".fits"
            os.system("julia ../julia_rascil_scripts/make_fullres.jl " + str(curr_lambda) + " " + tmp_psf_fn + " " \
                      + tmp_dirty_fn + " " + str(wavelet_type_idx) + " " + str(niter) + " " + curr_output_name)
            recon = util.fromfits(curr_output_name)

            if lp_filter is not None:
                recon = signal.fftconvolve(recon, lp_filter, mode='same')

            if stats_type == "psnr":
                snrs[i] = evaluation.compute_snr(gt, recon)
            elif stats_type == "ispace_resid_norm":
                recon = signal.fftconvolve(recon, psf, mode='same')
                identity_filt = split.bandpass(dirty, 0, dirty.shape[0] * 2)
                dirty = signal.fftconvolve(dirty, identity_filt, mode='same')

                snrs[i] = numpy.linalg.norm(recon - dirty)
            else:
                snrs[i] = evaluation.compute_snr(gt, recon)
            
        util.write_to_csv(snrs, path + "snr.dat")

    snrs = util.read_csv(path + "snr.dat")
    snrs = [float(x) for x in snrs]

    return snrs

def lambdatests_highvis(lambdas, path, dirty, psf, gt, lowres, wavelet_type_idx, niter, cut_center, cut_hw, runtests):
    if runtests:
        snrs = [0] * len(lambdas)

        tmp_dirty_fn = "tmp_dirty.fits"
        tmp_psf_fn = "tmp_psf.fits"
        tmp_low_fn = "tmp_low.fits"

        util.tofits(dirty, tmp_dirty_fn)
        util.tofits(psf, tmp_psf_fn)
        util.tofits(lowres, tmp_low_fn)

        vis_noise = numpy.mean(deconvolve.compute_windowed_var(dirty, 5))
        recon_noise = vis_noise / 1000

        for i, curr_lambda in enumerate(lambdas):
            curr_output_name = path + "lambda_" + str(curr_lambda) + ".fits"
            os.system("julia ../julia_rascil_scripts/make_multistep.jl " + str(curr_lambda) + " " + tmp_psf_fn + \
                      " " + tmp_dirty_fn + " " + tmp_low_fn + " " + str(wavelet_type_idx) + " " + str(niter) + " " + \
                      str(recon_noise) + " " + str(vis_noise) + " " + str(cut_center) + " " + str(cut_hw) + " " + curr_output_name)
            recon = fits.open(curr_output_name)[0].data
            snrs[i] = evaluation.compute_snr(gt, recon)
        
        util.write_to_csv(snrs, path + "snr.dat")

    snrs = util.read_csv(path + "snr.dat")
    snrs = [float(x) for x in snrs]

    return snrs


def deconvolve(dirty, psf, lowres, niter, wavelet_type_idx, curr_maj_iter, initial_lambda, lambda_mul, cut_center=20, cut_halfwidth=5, variance_window=5, recon_variance_factor=1, lowin=None, vis_variance=None, recon_variance=None):
    curr_lambda = initial_lambda * (lambda_mul ** curr_maj_iter)
    curr_lambda *= numpy.linalg.norm(dirty)

    tmp_psf_name = "tmp_psf.fits"
    tmp_res_name = "tmp_residual.fits"
    tmp_lowin_name = "tmp_lowin.fits"
    tmp_output_name = "tmp_output.fits"

    util.tofits(psf, tmp_psf_name)
    util.tofits(dirty, tmp_res_name)

    #low resolution step if no constraint
    if lowin is None:
        os.system("julia ../julia_rascil_scripts/make_lowres.jl " + str(curr_lambda) + " " + tmp_psf_name + " " + tmp_res_name + " " + str(wavelet_type_idx) + " " + str(niter) + " " + tmp_output_name)
    else:
        #auto calculate sigma and eta if none provided
        if vis_variance is None:
            vis_variance = numpy.mean(deconvolve.compute_windowed_var(dirty, variance_window))
            recon_variance = vis_variance / recon_variance_factor

        util.tofits(lowres, tmp_lowin_name)

        os.system("julia ../julia_rascil_scripts/make_multistep.jl " + str(curr_lambda) + " " + tmp_psf_name + " " + tmp_res_name + " " + tmp_lowin_name + " " + str(wavelet_type_idx) + " " + str(niter) + " " + \
                str(recon_variance) + " " + str(vis_variance) + " " + str(cut_center) + " " + str(cut_halfwidth) + " " + tmp_output_name)

    deconvolved = util.fromfits(tmp_output_name)

    return deconvolved