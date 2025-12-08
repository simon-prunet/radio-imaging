#!/usr/bin/env python

"""\
This file contains code that computes the weights for some given set of visibilities
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import numpy
import os
import time
import gc

from ska_sdp_func_python.imaging import create_image_from_visibility
from ska_sdp_datamodels.gridded_visibility import create_griddata_from_image
from ska_sdp_func_python.grid_data import grid_visibility_weight_to_griddata, griddata_visibility_reweight, griddata_merge_weights
from ska_sdp_func_python.visibility import convert_visibility_to_stokesI

from radioimaging.visibility import ingest

def compute_weights(vis, npixel, cellsize, weighting, robustness=0.0):
    """
    compute_weights computes the weights 

    :vis: measured unweighted visibilities, in the RASCIL format
    :npixel: image dimensions, the same is used for both x and y axis
    :cellsize: angular resolution of each pixel in radians
    :weighting: how visibilities should be weighted (eg. briggs, natural, uniform)
    :robustness: parameter for briggs weighting
    :return: weighted visibilities
    """

    if (weighting != "natural"):
        model = create_image_from_visibility(vis, cellsize=cellsize, npixel=npixel, polarisation_frame=vis.visibility_acc.polarisation_frame)
        grid_weights = create_griddata_from_image(model, polarisation_frame=model.image_acc.polarisation_frame)
        grid_weights = grid_visibility_weight_to_griddata(vis, grid_weights)
        vis = griddata_visibility_reweight(vis, grid_weights[0], weighting=weighting, robustness=robustness, sumwt=grid_weights[1])
    else:
        vis = griddata_visibility_reweight(vis, None, weighting=weighting)

    return vis


def compute_weights_griddata_from_ms(ms_name, channel_start, channel_end, data_descriptors, npixel, cellsize, bda=False):
    """
    compute_weights_griddata_from_ms computes the weights grid from some given measurement set
    This is typically used in conjunction with the compute residual or psf functions which ingest the measurement set by channel as well 

    :ms_name: filename of measurement set
    :channel_start: first channel to read of each data descriptor
    :channel_end: last channel to read of each data descriptor
    :data_descriptors: data descriptors to read from
    :npixel: image dimensions, the same is used for both x and y axis
    :cellsize: angular resolution of each pixel in radians
    :return: weights griddata and timing information for profiling
    """

    total_grid = None
    model = None

    read_timings = []
    polarization_timings = []
    grid_timings = []

    channels = range(channel_start, channel_end + 1)

    total_vis = 0

    for dd in data_descriptors:
        for curr_channel in channels:
            read_start = time.time()
            [vis], num_vis = ingest.create_visibility_from_ms(ms_name, start_chan=curr_channel, end_chan=curr_channel, selected_dds=[dd], use_weight_spec=bda)
            total_vis += num_vis
            pol_start = time.time()
            vis = convert_visibility_to_stokesI(vis)

            #if model is None:
            model = create_image_from_visibility(vis, cellsize=cellsize, npixel=npixel, polarisation_frame=vis.visibility_acc.polarisation_frame)

            grid_start = time.time()
            curr_grid_weights = create_griddata_from_image(model, polarisation_frame=model.image_acc.polarisation_frame)
            
            curr_grid_weights = grid_visibility_weight_to_griddata(vis, curr_grid_weights)
            
            if total_grid is None:
                total_grid = curr_grid_weights
            else:
                total_grid = griddata_merge_weights([total_grid, curr_grid_weights])

            channel_end = time.time()

            gc.collect()

            read_timings.append(pol_start - read_start)
            polarization_timings.append(grid_start - pol_start)
            grid_timings.append(channel_end - grid_start)

    read_total = sum(read_timings)
    polarization_total = sum(polarization_timings)
    weight_total = sum(grid_timings)

    gc.collect()

    return total_grid, [read_total, polarization_total, weight_total], total_vis