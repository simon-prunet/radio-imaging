#!/usr/bin/env python

"""\
Deconvolution code that takes care of most of the boilerplate code, typically calls some external code from either RASCIL, Julia, or
the python fista implementation
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import numpy
import os

from ska_sdp_func_python.imaging import create_image_from_visibility

from radioimaging.util import util

def compute_windowed_var(image, window):
    """
    compute_windowed_var computes the sliding windowed variance of a given image

    :image: input image
    :window: window size
    :return: image containing the per-pixel estimated variance
    """

    estimated_variance = numpy.zeros(image.shape)

    #convolving initial signal with a variance estimation kernel with a size 2xwindow_hsize+1
    for y in range(image.shape[1]):
        for x in range(image.shape[0]):
            start_x = max(x - window, 0)
            end_x = min(x + window + 1, image.shape[0])
            start_y = max(y - window, 0)
            end_y = min(y + window + 1, image.shape[1])

            estimated_variance[x, y] = numpy.var(image[start_x:end_x, start_y:end_y])
            
    return estimated_variance

def deconvolve_single(dirty, psf, niter, wavelet_type_idx, curr_maj_iter, initial_lambda, lambda_mul, script_root=""):
    """
    deconvolve_single deconvolves a full resolution image using the fista l1 deconvolution implemented in julia

    :dirty: dirty image
    :psf: point spread function
    :niter: number of fista iterations
    :wavelet_type_idx: wavelet index, 1 for daubechies and 2 for iuwt
    :curr_maj_iter: current major iteration, used to determine actual regularization parameter
    :initial_lambda: lambda used for the first major iteration
    :lambda_mul: multiplier to apply to lambda every major cycle
    :script_root: root directory of julia script
    :return: deconvolved full resolution image
    """
    res = numpy.array(dirty)
    np_psf = numpy.array(psf)

    curr_lambda = initial_lambda * (lambda_mul ** curr_maj_iter)
    curr_lambda *= numpy.linalg.norm(dirty)

    tmp_psf_name = "tmp_psf.fits"
    tmp_res_name = "tmp_residual.fits"
    tmp_output_name = "tmp_output.fits"

    util.tofits(psf, tmp_psf_name)
    util.tofits(dirty, tmp_res_name)

    os.system("julia --threads 32 " + script_root + "/make_fullres.jl " + str(curr_lambda) + " " + tmp_psf_name + " " + tmp_res_name + " " + str(wavelet_type_idx) + " " + str(niter) + " " + tmp_output_name)

    deconvolved = util.fromfits(tmp_output_name)

    return deconvolved


def deconvolve_multipartition_single(dirty, psf, niter, wavelet_type_idx, curr_maj_iter, initial_lambda, script_root=""):
    """
    deconvolve_multipartition_single deconvolves a full resolution image using the fista l1 deconvolution implemented in julia.
    This version differs from the original as it uses lambda max to regularize instead of lambda

    :dirty: dirty image
    :psf: point spread function
    :niter: number of fista iterations
    :wavelet_type_idx: wavelet index, 1 for daubechies and 2 for iuwt
    :curr_maj_iter: current major iteration, used to determine actual regularization parameter
    :initial_lambda: initial multiplier applied to lambda_max, increased as major cycles progresses
    :script_root: root directory of julia script
    :return: deconvolved full resolution image
    """
    res = numpy.array(dirty)
    np_psf = numpy.array(psf)

    tmp_psf_name = "tmp_psf.fits"
    tmp_res_name = "tmp_residual.fits"
    tmp_output_name = "tmp_output.fits"

    util.tofits(psf, tmp_psf_name)
    util.tofits(dirty, tmp_res_name)

    os.system("julia --threads 32  " + script_root + "/make_fullres_mp.jl " + str(initial_lambda) + " " + tmp_psf_name + " " + tmp_res_name + " " + str(wavelet_type_idx) + " " + str(niter) + " " + tmp_output_name)

    deconvolved = util.fromfits(tmp_output_name)

    return deconvolved

def deconvolve(step, dirty, psf, prev_estimates, niter, wavelet_type_idx, curr_maj_iter, initial_lambda, lambda_mul, cut_center, cut_halfwidth, variance_window, recon_variance_factor, script_root=""):
    """
    deconvolve deconvolves a low or high resolution image using the fista l1 deconvolution implemented in julia for the parallel
    interleaved reconstruction method

    :step: 0 for low, 1 for high
    :dirty: partial resolution dirty image
    :psf: partial resolution point spread function
    :prev_estimates: list of previous reconstructed images for each resolution 
    :niter: number of fista iterations
    :wavelet_type_idx: wavelet index, 1 for daubechies and 2 for iuwt
    :curr_maj_iter: current major iteration, used to determine actual regularization parameter
    :initial_lambda: initial multiplier applied to lambda, increased as major cycles progresses
    :lambda_mul: multiplier to apply to lambda every major cycle
    :cut_center: center of circle bisecting transition area, in pixel units
    :cut_halfwidth: halfwidth of transition area in pixel units
    :variance_window: window size for computing variance of dirty
    :recon_variance_factor: reconstruction variance is set to the dirty variance multiplied by this factor
    :script_root: root directory of julia script
    :return: deconvolved full resolution image from second major cycle onwards, otherwise a deconvolved partial resolution image
    """
    res = numpy.array(dirty)
    np_psf = numpy.array(psf)

    curr_lambda = initial_lambda

    if step == 0:
        curr_lambda = initial_lambda * (lambda_mul ** (curr_maj_iter))
    else:
        if curr_maj_iter > 1:
            curr_lambda = initial_lambda * (lambda_mul ** (curr_maj_iter - 1))

    tmp_psf_name = "tmp_psf_" + str(step) + ".fits"
    tmp_res_name = "tmp_residual_" + str(step) + ".fits"
    tmp_constraint_name = "tmp_constraint_" + str(step) + ".fits"
    tmp_output_name = "tmp_output_" + str(step) + ".fits"

    constraint = prev_estimates[1] - prev_estimates[0] if step == 0 else prev_estimates[0] - prev_estimates[1]

    vis_variance = 1#numpy.mean(compute_windowed_var(dirty, variance_window))
    constraint_variance = 1#numpy.mean(compute_windowed_var(constraint, variance_window))

    low_variance = vis_variance if step == 0 else constraint_variance
    high_variance = constraint_variance if step == 0 else vis_variance

    util.tofits(psf, tmp_psf_name)
    util.tofits(dirty, tmp_res_name)

    curr_lambda *= (numpy.linalg.norm(dirty) + numpy.linalg.norm(constraint))

    util.tofits(constraint, tmp_constraint_name)

    os.system("julia --threads 32  " + script_root + "/make_multistep_interleaved.jl " + str(curr_lambda) + " " + tmp_psf_name + " " + tmp_res_name + " " + tmp_constraint_name + " " + str(wavelet_type_idx) + " " + str(niter) + " " + \
            str(low_variance) + " " + str(high_variance) + " " + str(cut_center) + " " + str(cut_halfwidth)  + " " + str(step) + " " + str(curr_maj_iter) + " " + tmp_output_name)

    deconvolved = util.fromfits(tmp_output_name)

    return deconvolved



#interleaved deconvolution for multiple partitions
def deconvolve_multipartition(partition, dirty, psf, prev_estimates, niter, wavelet_type_idx, curr_maj_iter, initial_lambda, lambda_mul, ells, delta, variance_window, dirty_var, deconv_partitions, script_root=""):
    """
    deconvolve_multipartition deconvolves a partial resolution image using the fista l1 deconvolution implemented in julia for the parallel
    interleaved reconstruction method

    :partition: an index corresponding to the resolution
    :dirty: partial resolution dirty image
    :psf: partial resolution point spread function
    :prev_estimates: list of previous reconstructed images for each resolution 
    :niter: number of fista iterations
    :wavelet_type_idx: wavelet index, 1 for daubechies and 2 for iuwt
    :curr_maj_iter: current major iteration, used to determine actual regularization parameter
    :initial_lambda: initial multiplier applied to lambda_max, increased as major cycles progresses
    :lambda_mul: multiplier to apply to lambda_max every major cycle
    :ells: centers of circles bisecting transition areas, in pixel units
    :delta: halfwidth of transition areas in pixel units. This is for now assumed to be constant across all transition areas
    :variance_window: window size for computing variance of reconstructed images of each partition
    :dirty_var: variance of dirty image of current partition
    :deconv_partitions: list specifying whether to deconvolve a specific partition, 1 for yes, 0 for no. This is used for cases when
    the larger resolutions have very few pixels in fourier space, but are also close to fully sampled. In this case, the dirty can be
    used directly as the deconvolved image.
    :script_root: root directory of julia script
    :return: deconvolved full resolution image from second major cycle onwards, otherwise a deconvolved partial resolution image
    """
    res = numpy.array(dirty)
    np_psf = numpy.array(psf)

    tmp_psf_name = "tmp_psf_" + str(partition) + ".fits"
    tmp_res_name = "tmp_residual_" + str(partition) + ".fits"
    
    tmp_output_name = "tmp_output_" + str(partition) + ".fits"

    n_partitions = len(prev_estimates) if prev_estimates is not None else 1
    sigma2s = [0] * n_partitions

    constraint_param = ""
    sigma2s_param = ""
    deconv_param = ""

    if curr_maj_iter > 0:
        for i, est_image in enumerate(prev_estimates):
            tmp_constraint_name = "tmp_constraint_" + str(partition) + "_" + str(i) + ".fits"
            constraint_param += tmp_constraint_name + " "

            if i == partition:
                util.tofits(est_image, tmp_constraint_name)
                sigma2s[i] = dirty_var

                continue

            constraint_image = est_image - prev_estimates[partition]
            util.tofits(constraint_image, tmp_constraint_name)
            sigma2s[i] = numpy.mean(compute_windowed_var(constraint_image, variance_window))

    util.tofits(np_psf, tmp_psf_name)
    util.tofits(res, tmp_res_name)

    ells_param = ""

    for ell in ells:
        ells_param += str(ell) + " "

    for sigma2 in sigma2s:
        sigma2s_param += str(sigma2) + " "

    #partition is +1 due to julia being 1 indexed
    command = "julia --threads 32  " + script_root + "/make_multipartition.jl " + tmp_psf_name + " " + tmp_res_name + " " + str(initial_lambda) + " " + str(niter) + " " + str(n_partitions) + " " + str(partition + 1) \
         + " " + str(curr_maj_iter)  + " " + str(delta) + " " + tmp_output_name + " " + str(deconv_partitions[partition]) + " " + constraint_param + ells_param + sigma2s_param

    os.system(command)

    deconvolved = util.fromfits(tmp_output_name)

    return deconvolved


def deconvolve_multistep(dirty, psf, constraint, niter, wavelet_type_idx, curr_maj_iter, initial_lambda, lambda_mul, cut_center, cut_halfwidth, variance_window, recon_variance_factor, script_root=""):
    """
    deconvolve_multistep performs the second step deconvolution of the multistep reconstruction, using a reconstructed low-resolution
    image as an additional data fidelity constraint

    :dirty: high resolution dirty image
    :psf: high resolution point spread function
    :constraint: previously reconstructed low resolution image minus all already reconstructed low-resolution information
    :niter: number of fista iterations
    :wavelet_type_idx: wavelet index, 1 for daubechies and 2 for iuwt
    :curr_maj_iter: current major iteration, used to determine actual regularization parameter
    :initial_lambda: initial multiplier applied to lambda, increased as major cycles progresses
    :lambda_mul: multiplier to apply to lambda every major cycle
    :cut_center: center of circle bisecting transition area, in pixel units
    :cut_halfwidth: halfwidth of transition area in pixel units
    :variance_window: window size for computing variance of dirty
    :recon_variance_factor: reconstruction variance is set to the dirty variance multiplied by this factor
    :script_root: root directory of julia script
    :return: deconvolved full resolution image
    """

    res = numpy.array(dirty)
    np_psf = numpy.array(psf)

    curr_lambda = initial_lambda * (lambda_mul ** (curr_maj_iter))

    tmp_psf_name = "tmp_psf.fits"
    tmp_res_name = "tmp_residual.fits"
    tmp_constraint_name = "tmp_constraint.fits"
    tmp_output_name = "tmp_output.fits"

    vis_variance = numpy.mean(compute_windowed_var(dirty, variance_window))
    constraint_variance = numpy.mean(compute_windowed_var(constraint, variance_window))

    low_variance = constraint_variance
    high_variance = vis_variance

    util.tofits(psf, tmp_psf_name)
    util.tofits(dirty, tmp_res_name)

    curr_lambda *= (numpy.linalg.norm(dirty) + numpy.linalg.norm(constraint))

    util.tofits(constraint, tmp_constraint_name)

    os.system("julia --threads 32  " + script_root + "/make_multistep_interleaved.jl " + str(curr_lambda) + " " + tmp_psf_name + " " + tmp_res_name + " " + tmp_constraint_name + " " + str(wavelet_type_idx) + " " + str(niter) + " " + \
            str(low_variance) + " " + str(high_variance) + " " + str(cut_center) + " " + str(cut_halfwidth)  + " 1 " + str(curr_maj_iter + 1) + " " + tmp_output_name)

    deconvolved = util.fromfits(tmp_output_name)

    return deconvolved