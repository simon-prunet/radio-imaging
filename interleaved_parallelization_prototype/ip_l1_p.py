#!/usr/bin/env python

"""\
Script that performs parallel spatial frequency l1 reconstruction on some dataset for two partitions.
Example execution: mpiexec -n 3 python ip_l1_p.py <config_filename>
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

from mpi4py import MPI
import json
import numpy
import time
import sys
from pathlib import Path

from radioimaging.visibility import weights, residual, ingest
from radioimaging.util import util
from radioimaging.images import images, filters
from radioimaging.deconvolution import deconvolve

wavelet_type_dict = {"daubechies" : 0, "iuwt" : 1}
config_filename = sys.argv[1]

def combine(reconstructions, npixels, average=True):
    final_recon = numpy.zeros((npixels, npixels))

    for recon in reconstructions:
        while len(recon.shape) > 2:
            recon = recon[0]

        final_recon += recon

    if average:
        final_recon /= len(reconstructions)

    return final_recon

#master node, responsible for reading config, sending this to individual nodes, and gathering and combining reconstructed images
def master():
    recon_start = time.time()
    with open(config_filename) as f: 
        data = f.read() 

    config = json.loads(data)
    config = comm.bcast(config, root=0)

    output_dir = config["output_dir"] + "_l1_p/"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    timings_file = output_dir + "mc_timings"

    load_tokens = comm.gather(None, root=0)

    startup_end = time.time()
    util.write_to_csv([startup_end - recon_start], timings_file)

    curr_recon_vl = None
    curr_recon_vh = None

    for i in range(config["nmajcyc"]):
        mc_start = time.time()
        #for now doing this once every major cycle, but in practice only needed at the end of the reconstruction
        reconstructions = comm.gather(None, root=0)
        if curr_recon_vl is None:
            curr_recon_vl = reconstructions[1]
            curr_recon_vh = reconstructions[2]
        else:
            curr_recon_vl += reconstructions[1]
            curr_recon_vh += reconstructions[2]

        curr_recon = combine([curr_recon_vh, curr_recon_vl], config["npixels"], average=(i > 0))

        util.tofits(curr_recon, output_dir + "deconv_combined_" + str(i) + ".fits")
        util.tofits(reconstructions[1], output_dir + "deconv_vl_" + str(i) + ".fits")
        util.tofits(reconstructions[2], output_dir + "deconv_vh_" + str(i) + ".fits")

        mc_end = time.time()

        util.write_to_csv([mc_end - mc_start], timings_file)

#reconstruction nodes, responsible for reconstructing the individual full resolution images. For now assumes only 2 partitions
def recon(step):
    config = comm.bcast(None, root=0)

    ms_name = config["lowres-dataset"] if step == 0 else config["highres-dataset"]
    npixels = config["npixels"]
    cellsize = config["cellsize"]
    weighting = config["weighting"]
    wavelet_idx = wavelet_type_dict[config["wavelet_dict"]]
    robustness = config["robustness"]
    channel_start = int(config["channel_start"])
    channel_end = int(config["channel_end"])
    init_lambda = config["init_lambda_low"] if step == 0 else config["init_lambda_high"]
    lambda_mul = config["lambda_mul_low"] if step == 0 else config["lambda_mul_high"]
    data_descriptors = config["data_descriptors"]
    output_dir = config["output_dir"] + "_parallel/"

    breakdown_file = output_dir + "mc_timings_breakdown_" + str(step)

    weight_grid, weight_timings, num_vis = weights.compute_weights_griddata_from_ms(ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize)
    print(str(step) + " num vis: " + str(num_vis))
    psf, estimate, psf_timings, weight = residual.compute_psf_by_channel(ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize, weighting, robustness=robustness, weight_grid=weight_grid)
    util.tofits(psf.pixels.data[0,0,:,:], output_dir + "psf_" + str(step) + ".fits")
    other_estimate = ingest.create_image_from_ms(ms_name, npixels, cellsize)

    barrier_start = time.time()
    comm.gather(step, root=0)
    barrier_end = time.time()

    util.write_to_csv([num_vis], breakdown_file)
    util.write_to_csv(weight_timings, breakdown_file)
    util.write_to_csv(psf_timings, breakdown_file)
    util.write_to_csv([barrier_end - barrier_start], breakdown_file)

    for i in range(config["nmajcyc"]):
        send_start = time.time()
        #We use constraints after the first major cycle, which are injected into our objective function. These constraints are sent to and obtained from the other reconstruction node
        if i > 0:
            if step == 0:
                comm.Send(estimate.pixels.data, dest=2)
                comm.Recv(other_estimate.pixels.data, source=2)
            elif step == 1:
                comm.Recv(other_estimate.pixels.data, source=1)
                comm.Send(estimate.pixels.data, dest=1)
        send_end = time.time()

        resid, resid_timings = residual.compute_residual_from_ms(estimate, ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize, weighting, robustness=robustness, weight_grid=weight_grid)

        util.tofits(resid.pixels.data[0,0,:,:], output_dir + "residual_" + str(step) + "_" + str(i) + ".fits")

        deconv_start = time.time()
        constraint = [estimate.pixels.data[0,0,:,:], other_estimate.pixels.data[0,0,:,:]] if step == 0 else [other_estimate.pixels.data[0,0,:,:], estimate.pixels.data[0,0,:,:]]
        deconvolved = deconvolve.deconvolve(step, resid.pixels.data[0,0,:,:], psf.pixels.data[0,0,:,:], constraint, config["nfistaiter"], wavelet_idx, i, init_lambda, lambda_mul, \
            config["sep_center"], config["sep_hw"], config["visvar_window"], config["reconvar_factor"], script_root="julia")
        deconv_end = time.time()

        estimate = images.add_to_image(estimate, deconvolved)

        sendrecon_start = time.time()
        comm.gather(deconvolved, root=0)
        sendrecon_end = time.time()

        resid_timings.insert(0, send_end - send_start)
        resid_timings.append(deconv_end - deconv_start)
        resid_timings.append(sendrecon_end - sendrecon_start)

        util.write_to_csv(resid_timings, breakdown_file)


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    master()
elif rank == 1 or rank == 2:
    recon(rank - 1)