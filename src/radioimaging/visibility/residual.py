#!/usr/bin/env python

"""\
This file contains code that produces both residual images and psfs from either a given set of RASCIL visibilities, or from a measurement set filename.
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import numpy
import time
import gc

from ska_sdp_func_python.visibility import subtract_visibility
from ska_sdp_func_python.imaging import invert_ng, predict_ng, create_image_from_visibility, advise_wide_field
from ska_sdp_datamodels.gridded_visibility import create_griddata_from_image
from ska_sdp_func_python.grid_data import griddata_visibility_reweight
from ska_sdp_func_python.visibility import convert_visibility_to_stokesI
from ska_sdp_func_python.grid_data import grid_visibility_weight_to_griddata, griddata_visibility_reweight

from radioimaging.visibility import ingest
from radioimaging.images import images

def compute_residual(sky_estimate, vis, npixel, cellsize, include_weight=False):
    """
    compute_residual computes the residual between a sky estimate and a set of initial visibility measurements. It uses the improved w-stacking algorithm
    as the operator

    :sky_esimate: estimate of the sky
    :vis: measured visibilities, in the RASCIL format
    :npixel: image dimensions, the same is used for both x and y axis
    :cellsize: angular resolution of each pixel in radians
    :return: residual image as a RASCIL Image
    """


    #degridding, get visibilities of sky estimate
    vest = vis.copy(deep=True)
    vest = predict_ng(vest, sky_estimate, context='ng')

    #subtraction to obtain residual vis
    vres = subtract_visibility(vis, vest)

    #obtain dirty image and psf
    model = create_image_from_visibility(vres,cellsize=cellsize,npixel=npixel, polarisation_frame=vres.visibility_acc.polarisation_frame)

    dirty, sumwt = invert_ng(vres, model, context='ng')

    if include_weight:
        return dirty, sumwt
    else:
        return dirty



def compute_residual_from_ms(sky_estimate, ms_name, channel_start, channel_end, data_descriptors, npixel, cellsize, weighting, robustness=0, weight_grid=None, algorithm='ng', bda=False):
    """
    compute_residual_from_ms computes the residual between a sky estimate and a set of initial visibility measurements which is stored in a measurement set. 
    This is done on a channel-by-channel basis so that arbitrarily large measurement sets can be processed in memory so as long as each channel can fit in memory
    It uses the nifty-gridder implementation, thus the operators are limited to what is provided in this repository

    :sky_esimate: estimate of the sky
    :ms_name: filename of measurement set
    :channel_start: first channel to read of each data descriptor
    :channel_end: last channel to read of each data descriptor
    :data_descriptors: data descriptors to read from
    :npixel: image dimensions, the same is used for both x and y axis
    :cellsize: angular resolution of each pixel in radians
    :weighting: how visibilities should be weighted (eg. briggs, natural, uniform)
    :robustness: parameter for briggs weighting
    :weight_grid: weight grid used for uniform and briggs weighting
    :algorithm: gridding algorithm to use, needs to be one supported by nifty
    :return: residual image as a RASCIL Image and timing information for profiling
    """

    final_residual = None

    read_from_disk_timings = []
    convert_polarization_timings = []
    weight_timings = []
    malloc_timings = []
    predict_timings = []
    subtract_timings = []
    invert_timings = []
    add_timings = []

    channels = range(channel_start, channel_end + 1)

    weight = 0

    for i, dd in enumerate(data_descriptors):
        curr_weight_grid = weight_grid[i]
        
        for curr_channel in channels:
            read_start = time.time()
            [measured_vis], _ = ingest.create_visibility_from_ms(ms_name, start_chan=curr_channel, end_chan=curr_channel, selected_dds=[dd], use_weight_spec=bda)
            polarization_start = time.time()
            
            measured_vis = convert_visibility_to_stokesI(measured_vis)
            
            #measured_vis.frequency.data[0] = weight_grid[0].griddata_acc.griddata_wcs.sub([4]).wcs.crval[0]

            weight_start = time.time()
            measured_vis = griddata_visibility_reweight(measured_vis, curr_weight_grid[0], weighting=weighting, robustness=robustness, sumwt=curr_weight_grid[1])
            allocate_start = time.time()
            estimated_vis = measured_vis.copy(deep=True)
            predict_start = time.time()
            estimated_vis = predict_ng(estimated_vis, sky_estimate, context=algorithm)
            subtract_start = time.time()
            residual_vis = subtract_visibility(measured_vis, estimated_vis)

            if final_residual is None:
                final_residual = images.create_empty_image(measured_vis, npixel, cellsize)

            curr_residual_model = images.create_empty_image(measured_vis, npixel, cellsize)

            invert_start = time.time()
            channel_residual, sumwt = invert_ng(residual_vis, curr_residual_model, context=algorithm)

            prev_weight = weight
            weight += sumwt[0,0]

            add_start = time.time()
            if prev_weight == 0:
                final_residual = images.add_to_image(final_residual, channel_residual.pixels.data[0,0,:,:])
            else:
                final_residual.pixels.data[0,0,:,:] *= prev_weight / weight
                final_residual = images.add_to_image(final_residual, channel_residual.pixels.data[0,0,:,:] * (sumwt[0,0] / weight))

            channel_residual = None
            channel_end = time.time()

            gc.collect()

            read_from_disk_timings.append(polarization_start - read_start)
            convert_polarization_timings.append(weight_start - polarization_start)
            weight_timings.append(allocate_start - weight_start)
            malloc_timings.append(predict_start - allocate_start)
            predict_timings.append(subtract_start - predict_start)
            subtract_timings.append(invert_start - subtract_start)
            invert_timings.append(add_start - invert_start)
            add_timings.append(channel_end - add_start)

    

    read_from_disk_total = sum(read_from_disk_timings)
    convert_polarization_total = sum(convert_polarization_timings)
    weight_total = sum(weight_timings)
    malloc_total = sum(malloc_timings)
    predict_total = sum(predict_timings)
    subtract_total = sum(subtract_timings)
    invert_total = sum(invert_timings)
    add_total = sum(add_timings)

    gc.collect()

    return final_residual, [read_from_disk_total, convert_polarization_total, weight_total, malloc_total, predict_total, subtract_total, invert_total, add_total]


def compute_psf(vis, npixel, cellsize, include_weight_and_model=False):
    """
    compute_psf calculates a psf when given a set of visibilities
    the visibilities are assumed to have already been weighted

    :vis: visibilities given in the RASCIL format
    :npixel: image dimensions, the same is used for both x and y axis
    :cellsize: angular resolution of each pixel in radians
    :return: psf as a RASCIL Image 
    """

    model = create_image_from_visibility(vis,cellsize=cellsize,npixel=npixel, polarisation_frame=vis.visibility_acc.polarisation_frame)
    psf, sumwt = invert_ng(vis, model, context='ng', dopsf=True)

    if include_weight_and_model:
        return psf, model, sumwt
    else:
        return psf


def compute_psf_by_channel(ms_name, channel_start, channel_end, data_descriptors, npixel, cellsize, weighting, robustness=0, weight_grid=None, algorithm='ng', bda=False):
    """
    compute_psf_by_channel calculates a psf from visibilities stored in a measurement set. This is done on a channel-by-channel basis so that 
    arbitrarily large measurement sets can be processed in memory so as long as each channel can fit in memory
    It uses the nifty-gridder implementation, thus the operators are limited to what is provided in this repository

    :ms_name: filename of measurement set
    :channel_start: first channel to read of each data descriptor
    :channel_end: last channel to read of each data descriptor
    :data_descriptors: data descriptors to read from
    :npixel: image dimensions, the same is used for both x and y axis
    :cellsize: angular resolution of each pixel in radians
    :weighting: how visibilities should be weighted (eg. briggs, natural, uniform)
    :robustness: parameter for briggs weighting
    :weight_grid: weight grid used for uniform and briggs weighting
    :algorithm: gridding algorithm to use, needs to be one supported by nifty
    :return: psf as a RASCIL Image and timing information for profiling
    """

    model = None
    psf = None

    read_timings = []
    polarization_timings = []
    weight_timings = []
    invert_timings = []
    add_timings = []

    channels = range(channel_start, channel_end + 1)

    weight = 0

    for i, dd in enumerate(data_descriptors):
        curr_weight_grid = weight_grid[i]
        for curr_channel in channels:
            read_start = time.time()
            [vis], _ = ingest.create_visibility_from_ms(ms_name, start_chan=curr_channel, end_chan=curr_channel, selected_dds=[dd], use_weight_spec=bda)
            polar_start = time.time()
            vis = convert_visibility_to_stokesI(vis)
            weight_start = time.time()

            #vis.frequency.data[0] = curr_weight_grid[0].griddata_acc.griddata_wcs.sub([4]).wcs.crval[0]

            vis = griddata_visibility_reweight(vis, curr_weight_grid[0], weighting=weighting, robustness=robustness, sumwt=curr_weight_grid[1])

            model = create_image_from_visibility(vis, cellsize=cellsize, npixel=npixel, polarisation_frame=vis.visibility_acc.polarisation_frame)

            invert_start = time.time()
            curr_psf, sumwt = invert_ng(vis, model, context=algorithm, dopsf=True)
            prev_weight = weight
            weight += sumwt[0,0]

            add_start = time.time()
            if psf is None:
                psf = curr_psf
            else:
                psf.pixels.data[0,0,:,:] *= prev_weight / weight
                psf = images.add_to_image(psf, curr_psf.pixels.data[0,0,:,:] * (sumwt[0,0] / weight))
            channel_end = time.time()

            gc.collect()

            read_timings.append(polar_start - read_start)
            polarization_timings.append(weight_start - polar_start)
            weight_timings.append(invert_start - weight_start)
            invert_timings.append(add_start - invert_start)
            add_timings.append(channel_end - add_start)


    read_total = sum(read_timings)
    polarization_total = sum(polarization_timings)
    weight_total = sum(weight_timings)
    invert_total = sum(invert_timings)
    add_total = sum(add_timings)

    gc.collect()

    return psf, model, [read_total, polarization_total, weight_total, invert_total, add_total], weight

def visibilities_from_image(vt,fitsfile,scale_factor=1.0,return_cellsize=True, return_image=False, override_cellsize=False, ocellsize=1):
    """
    visibilities_from_image performs degridding to obtain visibility values from a fits file. Is essentially a wrapper function
    for degridding

    :vt: empty visibilities with correct uvw coordinates
    :fitsfile: filename of fits file
    :scale_factor: spatial scaling to manually modify object size
    :return_cellsize: whether or not to return the pixel angular resolution in radians
    :return_image: whether or not to return the image that is degridded
    :cellsize: angular resolution of each pixel in radians
    :override_cellsize: flag to determine if the cellsize should be overwritten
    :ocellsize: overwritten cell size angular resolution
    :return: degridded visibilities, as well as cell size and image depending on whether the flags were set
    """
    advice = advise_wide_field(vt, guard_band_image=3.0, delA=0.1, facets=1, 
        oversampling_synthesised_beam=4.0)
    cellsize = ocellsize if override_cellsize else scale_factor*advice['cellsize']

    im = images.create_image_from_fits(fitsfile,frequency=vt.frequency.data,cellsize=cellsize,phasecentre=vt.phasecentre)
    ivt = predict_ng(vt,im,context='ng')
    if return_cellsize and return_image:
        return(ivt,cellsize,im)
    elif (return_cellsize and not return_image):
        return(ivt,cellsize)
    elif (return_image and not return_cellsize):
        return(ivt,im)
    else:
        return(ivt)

def dirty_psf_from_visibilities(vt,cellsize,npix=512,weighting="uniform",robustness=0.0, override_cellsize=True):
    """
    dirty_psf_from_visibilities obtains a dirty image and a psf from some given set of visibilities. This is a wrapper function
    that calls gridding twice, once for the psf and once for the dirty

    :vt: visibilities
    :cellsize: angular resolution of pixels in radians
    :npix: number of pixels per dimension of the images. Assumed to be square
    :weighting: weighting scheme for the visibilities
    :robustness: parameter used for briggs weighting
    :override_cellsize: if this flag is not set, the ideal cellsize is automatically calculated and used rather than the passed
    cellsize
    :return: dirty image and psf
    """

    # First create empty rascil Image instance from visibilities
    model = create_image_from_visibility(vt,cellsize=cellsize,npixel=npix, override_cellsize=override_cellsize, polarisation_frame=vt.visibility_acc.polarisation_frame)
    print ("Model image plate scale (arcsec) is %e"%numpy.abs((model.image_acc.wcs.wcs.cdelt[0]*3600)))
    # Reweight visibilities if not natural weighting
    if (weighting != "natural"):
        grid_weights = create_griddata_from_image(model,polarisation_frame=model.image_acc.polarisation_frame)
        grid_weights = grid_visibility_weight_to_griddata(vt,grid_weights)
        vt=griddata_visibility_reweight(vt, grid_weights[0], weighting=weighting, 
                                        robustness=robustness, sumwt=grid_weights[1])

    dirty, sumwt = invert_ng(vt, model, context='ng')
    psf, sumwt   = invert_ng(vt, model, context='ng', dopsf=True)

    return (dirty,psf)