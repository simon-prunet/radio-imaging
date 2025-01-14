import json
import ip_helpers as iph
from rascil.processing_components import create_visibility_from_ms
import numpy
import time
import sys
from pathlib import Path

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

weight_grid, weight_timings, num_vis = iph.compute_weights_griddata_by_channel(low_ms_name, npixels, cellsize, channel_start, channel_end, data_descriptors)
psf, estimate_low, psf_timings = iph.compute_psf_by_channel(low_ms_name, npixels, cellsize, weight_grid, weighting, robustness, channel_start, channel_end, data_descriptors)

for i in range(config["nmajcyc"]):
    residual, resid_timings = iph.compute_residual_bychannel(estimate_low, low_ms_name, npixels, cellsize, weighting, robustness, weight_grid, channel_start, channel_end, data_descriptors)
    iph.tofits(residual.pixels.data[0,0,:,:], output_dir + "low_residual_" + str(i) + ".fits")
    deconvolved = iph.deconvolve_single(residual.pixels.data[0,0,:,:], psf.pixels.data[0,0,:,:], config["nfistaiter"], wavelet_idx, i, init_lambda_low, lambda_mul_low)
    iph.tofits(deconvolved, output_dir + "low_deconv_" + str(i) + ".fits")
    estimate_low = iph.add_to_image(estimate_low, deconvolved)
    iph.tofits(estimate_low.pixels.data[0,0,:,:], output_dir + "low_recon_" + str(i) + ".fits")

weight_grid, weight_timings, num_vis = iph.compute_weights_griddata_by_channel(high_ms_name, npixels, cellsize, channel_start, channel_end, data_descriptors)
psf, estimate_full, psf_timings = iph.compute_psf_by_channel(high_ms_name, npixels, cellsize, weight_grid, weighting, robustness, channel_start, channel_end, data_descriptors)

for i in range(config["nmajcyc"]):
    residual, resid_timings = iph.compute_residual_bychannel(estimate_full, high_ms_name, npixels, cellsize, weighting, robustness, weight_grid, channel_start, channel_end, data_descriptors)
    iph.tofits(residual.pixels.data[0,0,:,:], output_dir + "high_residual_" + str(i) + ".fits")
    deconvolved = iph.deconvolve_multistep(residual.pixels.data[0,0,:,:], psf.pixels.data[0,0,:,:], estimate_low.pixels.data[0,0,:,:], config["nfistaiter"], wavelet_idx, i, init_lambda_high, lambda_mul_high, \
            config["sep_center"], config["sep_hw"], config["visvar_window"], config["reconvar_factor"])
    iph.tofits(deconvolved, output_dir + "full_deconv_" + str(i) + ".fits")
    estimate_full = iph.add_to_image(estimate_full, deconvolved)
    iph.tofits(estimate_full.pixels.data[0,0,:,:], output_dir + "full_recon_" + str(i) + ".fits")
