import images
import helpers

from ska_sdp_func_python.visibility import subtract_visibility
from ska_sdp_func_python.imaging import invert_ng, predict_ng, create_image_from_visibility
from ska_sdp_datamodels.gridded_visibility import create_griddata_from_image
from ska_sdp_func_python.grid_data import grid_visibility_weight_to_griddata,griddata_visibility_reweight

from visibilities import compute_residual, compute_psf, compute_weights, create_empty_image

import numpy
import os

#this is currently a stub, to be implemented
#interleaved deconvolution, currently assumes 2 partitions
def deconvolve(step, dirty, psf, prev_estimates, niter, wavelet_type_idx, curr_maj_iter, initial_lambda, lambda_mul, cut_center, cut_halfwidth, variance_window, recon_variance_factor):
    res = numpy.array(dirty)
    np_psf = numpy.array(psf)

    curr_lambda = initial_lambda * (lambda_mul ** curr_maj_iter)
    curr_lambda *= numpy.linalg.norm(dirty)

    tmp_psf_name = "tmp_psf.fits"
    tmp_res_name = "tmp_residual.fits"
    tmp_constraint_name = "tmp_constraint.fits"
    tmp_output_name = "tmp_output.fits"

    vis_variance = numpy.mean(helpers.compute_windowed_var(dirty, variance_window))

    low_variance = vis_variance if step == 0 else vis_variance / recon_variance_factor
    high_variance = vis_variance / recon_variance_factor if step == 0 else vis_variance

    helpers.write_nparr_to_fits(psf, tmp_psf_name)
    helpers.write_nparr_to_fits(dirty, tmp_res_name)

    constraint = prev_estimates[1] - prev_estimates[0] if step == 0 else prev_estimates[0] - prev_estimates[1] 
    #helpers.plotNImages([constraint], ["constraint step " + str(step) + " maj cycle " + str(curr_maj_iter)], cmap="turbo")

    helpers.write_nparr_to_fits(constraint, tmp_constraint_name)

    os.system("julia ../julia_rascil_scripts/make_multistep_interleaved.jl " + str(curr_lambda) + " " + tmp_psf_name + " " + tmp_res_name + " " + tmp_constraint_name + " " + str(wavelet_type_idx) + " " + str(niter) + " " + \
            str(low_variance) + " " + str(high_variance) + " " + str(cut_center) + " " + str(cut_halfwidth)  + " " + str(step) + " " + str(curr_maj_iter) + " " + tmp_output_name)

    deconvolved = helpers.readFits(tmp_output_name)

    return deconvolved


#stub for now, implemented later
def add_to_image(image, nparr):
    image.pixels.data[0,0,:,:] += nparr

    return image