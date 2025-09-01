#!/usr/bin/env python

"""\
Script that performs the multistep reconstruction across two partitions.
Example execution: python ip_multistep.py <config_filename>
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import json
from rascil.processing_components import create_visibility_from_ms
import numpy
import time
import sys
from pathlib import Path

from radioimaging.visibility import weights, residual
from radioimaging.util import util
from radioimaging.images import images
from radioimaging.deconvolution import deconvolve

wavelet_type_dict = {"daubechies" : 0, "iuwt" : 1}
config_filename = sys.argv[1]


with open(config_filename) as f: 
    data = f.read() 

config = json.loads(data)

low_ms_name = config["lowres-dataset"] 
high_ms_name = config["highres-dataset"]
npixels = config["npixels"]
cellsize = config["cellsize"]
weighting = config["weighting"]
wavelet_idx = wavelet_type_dict[config["wavelet_dict"]]
robustness = config["robustness"]
channel_start = int(config["channel_start"])
channel_end = int(config["channel_end"])
init_lambda_low = config["init_lambda_low"]
lambda_mul_low = config["lambda_mul_low"]
init_lambda_high = config["init_lambda_high"]
lambda_mul_high = config["lambda_mul_high"]
data_descriptors = config["data_descriptors"]
output_dir = config["output_dir"] + "_multistep/"

Path(output_dir).mkdir(parents=True, exist_ok=True)

weight_grid, weight_timings, num_vis = weights.compute_weights_griddata_from_ms(low_ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize)
psf, estimate_low, psf_timings, weight = residual.compute_psf_by_channel(low_ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize, weighting, robustness=robustness, weight_grid=weight_grid)

for i in range(config["nmajcyc"]):
    resid, resid_timings = residual.compute_residual_from_ms(estimate_low, low_ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize, weighting, robustness=robustness, weight_grid=weight_grid)
    util.tofits(resid.pixels.data[0,0,:,:], output_dir + "low_residual_" + str(i) + ".fits")
    deconvolved = deconvolve.deconvolve_single(resid.pixels.data[0,0,:,:], psf.pixels.data[0,0,:,:], config["nfistaiter"], wavelet_idx, i, init_lambda_low, lambda_mul_low, script_root="julia")
    util.tofits(deconvolved, output_dir + "low_deconv_" + str(i) + ".fits")
    estimate_low = images.add_to_image(estimate_low, deconvolved)
    util.tofits(estimate_low.pixels.data[0,0,:,:], output_dir + "low_recon_" + str(i) + ".fits")

weight_grid, weight_timings, num_vis = weights.compute_weights_griddata_from_ms(high_ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize)
psf, estimate_full, psf_timings, weight = residual.compute_psf_by_channel(high_ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize, weighting, robustness=robustness, weight_grid=weight_grid)

for i in range(config["nmajcyc"]):
    resid, resid_timings = residual.compute_residual_from_ms(estimate_full, high_ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize, weighting, robustness=robustness, weight_grid=weight_grid)
    util.tofits(resid.pixels.data[0,0,:,:], output_dir + "high_residual_" + str(i) + ".fits")
    deconvolved = deconvolve.deconvolve_multistep(resid.pixels.data[0,0,:,:], psf.pixels.data[0,0,:,:], estimate_low.pixels.data[0,0,:,:], config["nfistaiter"], wavelet_idx, i, init_lambda_high, lambda_mul_high, \
            config["sep_center"], config["sep_hw"], config["visvar_window"], config["reconvar_factor"], script_root="julia")
    util.tofits(deconvolved, output_dir + "full_deconv_" + str(i) + ".fits")
    estimate_full = images.add_to_image(estimate_full, deconvolved)
    util.tofits(estimate_full.pixels.data[0,0,:,:], output_dir + "full_recon_" + str(i) + ".fits")
