#!/usr/bin/env python

"""\
Utility code used by various different modules and notebooks
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import numpy
import gc
from scipy import stats

from ska_sdp_func_python.visibility import subtract_visibility
from ska_sdp_func_python.imaging import invert_ng, predict_ng
from ska_sdp_func_python.grid_data import griddata_visibility_reweight
from ska_sdp_func_python.visibility import convert_visibility_to_stokesI
from ska_sdp_datamodels.visibility.vis_model import Visibility

from radioimaging.visibility import ingest
from radioimaging.images import images

def compute_snr(gt, recon):
    """
    compute_snr computes the signal to noise ratio of some reconstructed image against a ground truth

    :gt: ground truth
    :recon: reconstructed image
    :return: snr
    """
    difnorm = numpy.linalg.norm(gt-recon)
    if difnorm == 0:
        return 0
        
    return 10 * numpy.log10(numpy.linalg.norm(gt) / difnorm)

def compute_rmse(gt, recon):
    """
    compute_rmse computes the root mean squared error of some reconstructed image against a ground truth

    :gt: ground truth
    :recon: reconstructed image
    :return: rmse
    """
    return numpy.sqrt(numpy.mean((gt-recon) ** 2))


def compute_jackknifed_residual_bychannel(sky_estimate, ms_name, npixel, cellsize, weighting, robustness, weight_grid, channel_start, channel_end, data_descriptors, algorithm='ng'):
    """
    compute_jackknifed_residual_bychannel computes an ideal residual image by randomly inverting half of the visibilities, this is 
    used to compare against reconstructions of real datasets where the ground truths are not available.
    This function performs this channel by channel to account for measurement sets that are too large to fit into memory at once

    :sky_estimate: empty image template with the same dimensionality as the residual
    :ms_name: measurement set containing the visibilities
    :npixel: residual image dimensions, assumed to be square
    :cellsize: angular resolution per pixel
    :weighting: weighting strategy
    :robustness: used for briggs weighting
    :weight_grid: precomputed weight grid used for briggs and uniform weighting schemes
    :channel_start: first channel to read
    :channel_end: last channel to read
    :data_descriptors: spectral windows to read
    :algorithm: gridding algorithm
    """

    final_residual = None

    channels = range(channel_start, channel_end + 1)
    rng = numpy.random.default_rng(42)

    for dd in data_descriptors:
        for curr_channel in channels:
            [measured_vis], _ = ingest.create_visibility_from_ms(ms_name, start_chan=curr_channel, end_chan=curr_channel, selected_dds=[dd])
            measured_vis = convert_visibility_to_stokesI(measured_vis)

            for (i, j, k, l), vis in numpy.ndenumerate(measured_vis.vis):
                if rng.random() > 0.5:
                    measured_vis.vis.data[i, j, k, l] *= -1

            measured_vis = griddata_visibility_reweight(measured_vis, weight_grid[0], weighting=weighting, robustness=robustness, sumwt=weight_grid[1])

            estimated_vis = measured_vis.copy(deep=True)
            estimated_vis = predict_ng(estimated_vis, sky_estimate, context=algorithm)
            residual_vis = subtract_visibility(measured_vis, estimated_vis)

            if final_residual is None:
                final_residual = images.create_empty_image(measured_vis, npixel, cellsize)

            channel_residual, sumwt = invert_ng(residual_vis, final_residual, context=algorithm)

            final_residual = images.add_to_image(final_residual, channel_residual.pixels.data[0,0,:,:])

            gc.collect()

    gc.collect()

    return final_residual


def t_test_image(recon_resid, target, window_size, test_type="wasserstein"):
    """
    t_test_image performs a per-pixel statistical test between a reconstruction and target residual.
    Uses a sliding window to obtain the distributions as only one realization is assumed 

    :recon_resid: reconstruction residual
    :target: target residual
    :window_size: sliding window size (assumed to be square)
    :test_type: type of test, only wasserstein, mwu, and student are currently supported
    :return: image of either distances or pvalues, depending on the test type used
    """
    assert(recon_resid.shape == target.shape)
    assert(len(recon_resid.shape) == 2)

    output_pvals = numpy.zeros(recon_resid.shape)

    for y in range(output_pvals.shape[1]):
        for x in range(output_pvals.shape[0]):
            start_x = max(x - window_size, 0)
            end_x = min(x + window_size + 1, output_pvals.shape[0])
            start_y = max(y - window_size, 0)
            end_y = min(y + window_size + 1, output_pvals.shape[1])

            recon_window = recon_resid[start_x:end_x, start_y:end_y]
            target_window = target[start_x:end_x, start_y:end_y]

            if test_type == "wasserstein":
                ttest_res = stats.wasserstein_distance(recon_window.flatten(), target_window.flatten())
                output_pvals[x, y] = ttest_res 
            elif test_type == "mwu":
                ttest_res = stats.mannwhitneyu(recon_window.flatten(), target_window.flatten())
                output_pvals[x, y] = ttest_res.pvalue
            elif test_type == "student":
                ttest_res = stats.ttest_ind(recon_window.flatten(), target_window.flatten(), equal_var=False)
                output_pvals[x, y] = ttest_res.pvalue
            else:
                print("test not supported")

    return output_pvals, numpy.linalg.norm(output_pvals)

def compute_maxabserr(gt, recon):
    """
    compute_maxabserr computes the maximum absolute error of some reconstructed image against a ground truth

    :gt: ground truth
    :recon: reconstructed image
    :return: max absolute error
    """
    return numpy.max(numpy.abs(gt - recon))

def compute_ssim(gt, recon):
    """
    compute_ssim computes the ssim of some reconstructed image against a ground truth

    :gt: ground truth
    :recon: reconstructed image
    :return: ssim
    """
    gt_mean = numpy.mean(gt)
    recon_mean = numpy.mean(recon)
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