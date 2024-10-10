import json
import ip_helpers as iph
import numpy
import time
import sys
from pathlib import Path

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
output_dir = config["output_dir"] + "_serial/"

Path(output_dir).mkdir(parents=True, exist_ok=True)


timings_file = output_dir + "mc_timings"
breakdown_file = output_dir + "mc_timings_breakdown"

recon_start = time.time()

weight_grid, weight_timings, num_vis = iph.compute_weights_griddata_by_channel(ms_name, npixels, cellsize, channel_start, channel_end, data_descriptors)
print("num vis: " + str(num_vis))

psf, estimate, psf_timings = iph.compute_psf_by_channel(ms_name, npixels, cellsize, weight_grid, weighting, robustness, channel_start, channel_end, data_descriptors)

iph.write_to_csv([num_vis], breakdown_file)
iph.write_to_csv(weight_timings, breakdown_file)
iph.write_to_csv(psf_timings, breakdown_file)

iph.tofits(psf.pixels.data[0,0,:,:], output_dir + "psf.fits")

init_lambda = config["init_lambda_full"]
lambda_mul = config["lambda_mul_full"]

nmaj = config["nmajcyc"] * 2

mc_start = time.time()
iph.write_to_csv([mc_start - recon_start], timings_file)

for i in range(nmaj):
    curr_mc_start = time.time()
    residual, resid_timings = iph.compute_residual_bychannel(estimate, ms_name, npixels, cellsize, weighting, robustness, weight_grid, channel_start, channel_end, data_descriptors)

    iph.tofits(residual.pixels.data[0,0,:,:], output_dir + "residual_" + str(i) + ".fits")

    deconvolve_start = time.time()
    deconvolved = iph.deconvolve_single(residual.pixels.data[0,0,:,:], psf.pixels.data[0,0,:,:], config["nfistaiter"], wavelet_idx, i, init_lambda, lambda_mul)

    iph.tofits(deconvolved, output_dir + "deconvolved_" + str(i) + ".fits")

    estimate = iph.add_to_image(estimate, deconvolved)

    curr_mc_end = time.time()

    resid_timings.append(curr_mc_end - deconvolve_start)

    iph.write_to_csv([curr_mc_end - curr_mc_start], timings_file)
    iph.write_to_csv(resid_timings, breakdown_file)

iph.tofits(estimate.pixels.data[0, 0, ...], output_dir + "final_deconvolved.fits")