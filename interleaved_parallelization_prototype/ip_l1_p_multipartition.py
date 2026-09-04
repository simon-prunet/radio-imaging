#!/usr/bin/env python

"""\
Script that performs parallel spatial frequency l1 reconstruction on some dataset for multiple partitions. The number of processes allocated needs to be the number
of partitions +1 to account for the master node.
Example execution for 3 partitions: mpiexec -n 4 python ip_l1_p_multipartition.py <config_filename>
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"


from mpi4py import MPI
import json
from rascil.processing_components import create_visibility_from_ms
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

    output_dir = config["output_dir"] + "pl1/"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    timings_file = output_dir + "mc_timings"

    wgts = comm.allgather(None)

    comm.Barrier()

    startup_end = time.time()
    util.write_to_csv([startup_end - recon_start], timings_file)

    for i in range(config["nmajcycl1"]):
        mc_start = time.time()
        all_estimates = comm.allgather(None)

        recons_combined = combine(all_estimates[1:], config["npixels"], average=(i > 0))

        util.tofits(recons_combined, output_dir + "recon_combined_" + str(i) + ".fits")

        for j, curr_recon in enumerate(all_estimates[1:]):
            util.tofits(curr_recon, output_dir + "recon_" + str(j) + "_" + str(i) + ".fits")

        mc_end = time.time()

        util.write_to_csv([mc_end - mc_start], timings_file)

#reconstruction nodes, responsible for reconstructing the individual full resolution images. For now assumes only 2 partitions
def recon(partition):
    with open(config_filename) as f: 
        data = f.read() 

    config = json.loads(data)

    all_datasets = config["datasets"]
    ms_name = all_datasets[partition]

    npixels = config["npixels"]
    cellsize = config["cellsize"]
    weighting = config["weighting"]
    wavelet_idx = wavelet_type_dict[config["wavelet_dict"]]
    robustness = config["robustness"]
    channel_start = int(config["channel_start"])
    channel_end = int(config["channel_end"])
    data_descriptors = list(range(int(config["data_descriptor_start"]), int(config["data_descriptor_end"]) + 1))
    bda = config["bda"] if config["bda"] is not None else False

    ells = config["ells"]
    ells.sort()

    delta = config["delta"]

    variance_window = config["visvar_window"]

    init_lambdas = config["init_lambdas"]
    init_lambda = init_lambdas[partition]
    step = (1 - init_lambda) / config["nmajcycl1"]

    lambdas = [init_lambda + i * step for i in range(config["nmajcycl1"])]

    n_fista_iter = config["nfistaiter"]
    
    output_dir = config["output_dir"] + "pl1/"

    k = config["lambda_growth_steepness"]

    deconv_partitions = config["deconv_partitions"]

    breakdown_file = output_dir + "mc_timings_breakdown_" + str(partition)

    weight_grid, weight_timings, num_vis = weights.compute_weights_griddata_from_ms(ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize, bda=bda)
    psf, estimate, psf_timings, weight = residual.compute_psf_by_channel(ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize, weighting, robustness=robustness, weight_grid=weight_grid, bda=bda)

    #leaving it here because it might be useful later
    wgts = numpy.zeros(len(init_lambdas))
    wgts = comm.allgather(weight)
    total_weight = numpy.sum(wgts[1:])

    psfsum = numpy.sum(psf)

    util.tofits(psf.pixels.data[0,0,:,:], output_dir + "psf_" + str(partition) + ".fits")

    #synchronization barrier primarily used to start the major cycles at the same time. This isn't strictly needed, but is handy for taking stats, as we can determine the amount
    #of waiting time before the first major cycle
    barrier_start = time.time()
    comm.Barrier()
    barrier_end = time.time()

    util.write_to_csv([num_vis], breakdown_file)
    util.write_to_csv(weight_timings, breakdown_file)
    util.write_to_csv(psf_timings, breakdown_file)
    util.write_to_csv([barrier_end - barrier_start], breakdown_file)

    sendrecon_start = sendrecon_end = 0
    first_res_var = 0

    pmc = config["nmajcycl1"]

    for i in range(config["nmajcycl1"]):
        send_start = time.time()
        #We use constraints after the first major cycle, which are injected into our objective function. These constraints are sent to and obtained from the other reconstruction node
        prev_estimates = None
        if i > 0:
            sendrecon_start = time.time()
            prev_estimates = comm.allgather(estimate.pixels.data[0,0,:,:])
            prev_estimates = prev_estimates[1:]
            sendrecon_end = time.time()


        send_end = time.time()

        resid, resid_timings = residual.compute_residual_from_ms(estimate, ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize, weighting, robustness=robustness, weight_grid=weight_grid, bda=bda)

        if i == 0:
            first_res_var = numpy.mean(deconvolve.compute_windowed_var(resid.pixels.data[0,0,:,:], variance_window))

        util.tofits(resid.pixels.data[0,0,:,:], output_dir + "residual_" + str(partition) + "_" + str(i) + ".fits")

        t = float(i) / (float(pmc) - 1)
        curr_lambda = init_lambda + (1 - init_lambda) * ((numpy.exp(k * t) - 1) / (numpy.exp(k) - 1))

        deconv_start = time.time()
        deconvolved = deconvolve.deconvolve_multipartition(partition, resid.pixels.data[0,0,:,:], psf.pixels.data[0,0,:,:], prev_estimates, n_fista_iter, i, curr_lambda, \
            ells, delta, variance_window, first_res_var, deconv_partitions)
        deconv_end = time.time()

        util.tofits(deconvolved, output_dir + "deconv_" + str(partition) + "_" + str(i) + ".fits")

        estimate = images.add_to_image(estimate, deconvolved)

        resid_timings.insert(0, send_end - send_start)
        resid_timings.append(deconv_end - deconv_start)
        resid_timings.append(sendrecon_end - sendrecon_start)

        util.write_to_csv(resid_timings, breakdown_file)

    #this is mainly to send to master
    prev_estimates = comm.allgather(estimate.pixels.data[0,0,:,:])


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    master()
else:
    recon(rank - 1)