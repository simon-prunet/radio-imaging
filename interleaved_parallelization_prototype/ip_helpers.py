#helper file for interleaved prototype. This may contain some duplicated code from the various python scripts in ../notebooks_and_pyscripts. This is done primarily to have a localized standalone prototype

import astropy
from astropy.io import fits

from ska_sdp_func_python.visibility import subtract_visibility
from ska_sdp_func_python.imaging import invert_ng, predict_ng, create_image_from_visibility, advise_wide_field
from ska_sdp_datamodels.gridded_visibility import create_griddata_from_image
from ska_sdp_func_python.grid_data import grid_visibility_weight_to_griddata,griddata_visibility_reweight
import numpy
import os

def tofits(data, filename):
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(filename, overwrite=True)
    hdulist.close()

def fromfits(filename):
    dat = fits.open(filename)[0].data
    while len(dat.shape) > 2:
        dat = dat[0]

    return dat

def compute_windowed_var(image, window):
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

def compute_residual(sky_estimate, vis, npixel, cellsize):
    #degridding, get visibilities of sky estimate
    vest = vis.copy(deep=True)
    vest = predict_ng(vest, sky_estimate, context='ng')

    #subtraction to obtain residual vis
    vres = subtract_visibility(vis, vest)

    #obtain dirty image and psf
    model = create_image_from_visibility(vres,cellsize=cellsize,npixel=npixel)
    dirty, sumwt = invert_ng(vres, model, context='ng')

    return dirty

def compute_psf(vis, npixel, cellsize):
    model = create_image_from_visibility(vis,cellsize=cellsize,npixel=npixel)
    psf, sumwt = invert_ng(vis, model, context='ng', dopsf=True)

    return psf

def compute_weights(vis, npixel, cellsize, weighting, robustness=0.0):
    if (weighting != "natural"):
        model = create_image_from_visibility(vis, cellsize=cellsize, npixel=npixel)
        grid_weights = create_griddata_from_image(model, polarisation_frame=model.image_acc.polarisation_frame)
        grid_weights = grid_visibility_weight_to_griddata(vis, grid_weights)
        vis = griddata_visibility_reweight(vis, grid_weights[0], weighting=weighting, robustness=robustness, sumwt=grid_weights[1])
    else:
        vis = griddata_visibility_reweight(vis, None, weighting=weighting)

    return vis

def create_empty_image(vis, npixel, cellsize):
    return create_image_from_visibility(vis, npixel=npixel, cellsize=cellsize)

#interleaved deconvolution, currently assumes 2 partitions
def deconvolve(step, dirty, psf, prev_estimates, niter, wavelet_type_idx, curr_maj_iter, initial_lambda, lambda_mul, cut_center, cut_halfwidth, variance_window, recon_variance_factor):
    res = numpy.array(dirty)
    np_psf = numpy.array(psf)

    curr_lambda = initial_lambda * (lambda_mul ** curr_maj_iter)
    curr_lambda *= numpy.linalg.norm(dirty)

    tmp_psf_name = "tmp_psf_" + str(step) + ".fits"
    tmp_res_name = "tmp_residual_" + str(step) + ".fits"
    tmp_constraint_name = "tmp_constraint_" + str(step) + ".fits"
    tmp_output_name = "tmp_output_" + str(step) + ".fits"

    vis_variance = numpy.mean(compute_windowed_var(dirty, variance_window))

    low_variance = vis_variance if step == 0 else vis_variance / recon_variance_factor
    high_variance = vis_variance / recon_variance_factor if step == 0 else vis_variance

    tofits(psf, tmp_psf_name)
    tofits(dirty, tmp_res_name)

    constraint = prev_estimates[1] - prev_estimates[0] if step == 0 else prev_estimates[0] - prev_estimates[1] 

    tofits(constraint, tmp_constraint_name)

    os.system("julia julia/make_multistep_interleaved.jl " + str(curr_lambda) + " " + tmp_psf_name + " " + tmp_res_name + " " + tmp_constraint_name + " " + str(wavelet_type_idx) + " " + str(niter) + " " + \
            str(low_variance) + " " + str(high_variance) + " " + str(cut_center) + " " + str(cut_halfwidth)  + " " + str(step) + " " + str(curr_maj_iter) + " " + tmp_output_name)

    deconvolved = fromfits(tmp_output_name)

    return deconvolved


#stub for now, implemented later
def add_to_image(image, nparr):
    image.pixels.data[0,0,:,:] += nparr

    return image