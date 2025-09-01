#!/usr/bin/env python

"""\
Script that performs the l1 reconstruction on some dataset.
Example execution: python ip_l1.py <config_filename>
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import json
import numpy
import time
import sys
from pathlib import Path

from radioimaging.util import util
from radioimaging.visibility import residual, weights
from radioimaging.images import images
from radioimaging.deconvolution import deconvolve

wavelet_type_dict = {"daubechies" : 0, "iuwt" : 1}
config_filename = sys.argv[1]

with open(config_filename) as f: 
    data = f.read()

config = json.loads(data)

ms_name = config["full-dataset"]
npixels = config["npixels"]
cellsize = config["cellsize"]
weighting = config["weighting"]
robustness = config["robustness"]
channel_start = int(config["channel_start"])
channel_end = int(config["channel_end"])
wavelet_idx = wavelet_type_dict[config["wavelet_dict"]]
data_descriptors = config["data_descriptors"]
output_dir = config["output_dir"] + "_l1/"

Path(output_dir).mkdir(parents=True, exist_ok=True)


timings_file = output_dir + "mc_timings"
breakdown_file = output_dir + "mc_timings_breakdown"

recon_start = time.time()

weight_grid, weight_timings, num_vis = weights.compute_weights_griddata_from_ms(ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize)
print("num vis: " + str(num_vis))

psf, estimate, psf_timings, weight = residual.compute_psf_by_channel(ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize, weighting, robustness=robustness, weight_grid=weight_grid)

util.write_to_csv([num_vis], breakdown_file)
util.write_to_csv(weight_timings, breakdown_file)
util.write_to_csv(psf_timings, breakdown_file)

util.tofits(psf.pixels.data[0,0,:,:], output_dir + "psf.fits")

init_lambda = config["init_lambda_full"]
lambda_mul = config["lambda_mul_full"]

nmaj = config["nmajcyc"] + 1

mc_start = time.time()
util.write_to_csv([mc_start - recon_start], timings_file)

for i in range(nmaj):
    curr_mc_start = time.time()
    resid, resid_timings = residual.compute_residual_from_ms(estimate, ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize, weighting, robustness=robustness, weight_grid=weight_grid)

    util.tofits(resid.pixels.data[0,0,:,:], output_dir + "residual_" + str(i) + ".fits")

    deconvolve_start = time.time()
    deconvolved = deconvolve.deconvolve_single(resid.pixels.data[0,0,:,:], psf.pixels.data[0,0,:,:], config["nfistaiter"], wavelet_idx, i, init_lambda, lambda_mul, script_root="julia")

    util.tofits(deconvolved, output_dir + "deconvolved_" + str(i) + ".fits")

    estimate = images.add_to_image(estimate, deconvolved)

    curr_mc_end = time.time()

    resid_timings.append(curr_mc_end - deconvolve_start)

    util.write_to_csv([curr_mc_end - curr_mc_start], timings_file)
    util.write_to_csv(resid_timings, breakdown_file)

util.tofits(estimate.pixels.data[0, 0, ...], output_dir + "final_deconvolved.fits")