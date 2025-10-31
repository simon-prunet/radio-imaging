#!/usr/bin/env python

"""\
Script that performs parallel spatial frequency ms-clean reconstruction on some dataset for multiple partitions. The number of processes allocated needs to be the number
of partitions +1 to account for the master node.
Example execution for 3 partitions: mpiexec -n 4 python ip_msclean_p_multipartition.py <config_filename>
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"


from mpi4py import MPI
import json
import numpy
import time
import sys
from pathlib import Path
from ska_sdp_func_python.image.cleaners import msclean
import gc

from radioimaging.visibility import weights, residual, ingest
from radioimaging.util import util
from radioimaging.images import images, filters

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

    output_dir = config["output_dir"] + "pmsc/"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    timings_file = output_dir + "mc_timings"

    psfs = comm.allgather(None)
    wgts = comm.allgather(None)


    comm.Barrier()

    startup_end = time.time()
    util.write_to_csv([startup_end - recon_start], timings_file)

    for i in range(config["nmajcycmsc"]):
        mc_start = time.time()
        comm.Barrier()

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
    robustness = config["robustness"]
    channel_start = int(config["channel_start"])
    channel_end = int(config["channel_end"])
    data_descriptors = range(int(config["data_descriptor_start"]), int(config["data_descriptor_end"]) + 1)
    bda = config["bda"] if config["bda"] is not None else False

    output_dir = config["output_dir"] + "pmsc/"
    
    ells = config["ells"]
    ells.sort()
    delta = config["delta"]

    var_window = config["visvar_window"]

    msc_niter = int(config["msclean_iter"])
    thresh = config["clean_thresh"]
    #scales = config["clean_scales"]
    #scales.sort()
    scales = config["fullres_clean_scales"]
    fracthresh = config["clean_fracthresh"]
    gain = config["clean_gain"]

    deconv_partitions = config["deconv_partitions"]

    breakdown_file = output_dir + "mc_timings_breakdown_" + str(partition)

    weight_grid, weight_timings, num_vis = weights.compute_weights_griddata_from_ms(ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize, bda=bda)
    psf, estimate, psf_timings, weight = residual.compute_psf_by_channel(ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize, weighting, robustness=robustness, weight_grid=weight_grid, bda=bda)
    util.tofits(psf.pixels.data[0,0,:,:], output_dir + "psf_" + str(partition) + ".fits")

    sigma2s = [1] * (len(ells) + 1)

    frs, _, _ = filters.create_filters_mstep(psf.pixels.data.shape[-1] // 2, [delta]*len(ells), ells, sigma2s)
    frns = [numpy.array(fr) for fr in frs]
    fs2ds = [filters.freq1d_to_radial2d(f1d, psf.pixels.data.shape[-1])[1] for f1d in frns]

    wgts = numpy.zeros(len(sigma2s))

    #share psfs, needed for all major-cycles after the first
    psf_send_start = time.time()
    psfs = comm.allgather(psf.pixels.data[0,0,:,:])
    psfs = psfs[1:]
    wgts = comm.allgather(weight)
    wgts = wgts[1:]
    psf_send_end = time.time()

    total_weight = sum(wgts)

    wgts = [weight / total_weight for weight in wgts]

    joint_psf = util.convolve2d(psfs[0] * wgts[0], fs2ds[0])
    for i, currpsf in enumerate(psfs[1:]):
        joint_psf += util.convolve2d(currpsf * wgts[i+1], fs2ds[i+1])

    util.tofits(joint_psf, output_dir + "joint_psf" + str(partition) + ".fits")

    barrier_start = time.time()
    comm.Barrier()
    barrier_end = time.time()

    util.write_to_csv([num_vis], breakdown_file)
    util.write_to_csv(weight_timings, breakdown_file)
    util.write_to_csv(psf_timings, breakdown_file)
    util.write_to_csv([psf_send_end - psf_send_start], breakdown_file)
    util.write_to_csv([barrier_end - barrier_start], breakdown_file)

    first_mc_scales = config["first_mc_clean_scales"][partition]

    skip_deconv = False#deconv_partitions[partition] == 0

    for i in range(config["nmajcycmsc"]):
        gc.collect()
        #We use constraints after the first major cycle, which are injected into our objective function. These constraints are sent to and obtained from the other reconstruction node
        prev_estimates = None
        if i > 0:
            #technically the gather is blocking, this piece of code is just so that idling time and tranfer time counts are separated
            idle_start = time.time()
            comm.Barrier()
            idle_end = time.time()

            sendrecon_start = time.time()
            prev_estimates = comm.allgather(estimate.pixels.data[0,0,:,:])
            prev_estimates = prev_estimates[1:]
            sendrecon_end = time.time()

        resid, resid_timings = residual.compute_residual_from_ms(estimate, ms_name, channel_start, channel_end, data_descriptors, npixels, cellsize, weighting, robustness=robustness, weight_grid=weight_grid, bda=bda)

        util.tofits(resid.pixels.data[0,0,:,:], output_dir + "residual_" + str(partition) + "_" + str(i) + ".fits")

        deconv_start = time.time()

        curr_residual = curr_psf = None
        deconvolved = None

        if skip_deconv:
            if i > 0:
                for j, sigma2 in enumerate(sigma2s):
                    if j != partition:
                        if deconvolved is None:
                            deconvolved = util.convolve2d(prev_estimates[j], fs2ds[j])
                        else:
                            deconvolved += util.convolve2d(prev_estimates[j], fs2ds[j])
            else:
                deconvolved = resid.pixels.data[0,0,:,:] / numpy.sum(psf["pixels"].data[0, 0, :, :])
        else:
            if i > 0:
                curr_residual = util.convolve2d(wgts[partition] * resid["pixels"].data[0, 0, :, :], fs2ds[partition])
                for j, sigma2 in enumerate(sigma2s):
                    if j != partition:
                        constraint = prev_estimates[j] - prev_estimates[partition]
                        constraint = util.convolve2d(constraint, wgts[j] * psfs[j])
                        constraint = util.convolve2d(constraint, fs2ds[j])

                        curr_residual += constraint

                curr_psf = joint_psf
            else:
                curr_residual = resid.pixels.data[0,0,:,:]
                curr_psf = psf["pixels"].data[0, 0, :, :]

            util.tofits(curr_residual, output_dir + "curr_residual_" + str(partition) + "_" + str(i) + ".fits")

            deconvolved, _ = msclean(curr_residual, curr_psf, None, None, gain, thresh, msc_niter, first_mc_scales if i == 0 else scales, fracthresh)

        deconv_end = time.time()

        util.tofits(deconvolved, output_dir + "deconv_" + str(partition) + "_" + str(i) + ".fits")

        estimate = images.add_to_image(estimate, deconvolved)

        resid_timings.append(deconv_end - deconv_start)
        resid_timings.append(sendrecon_end - sendrecon_start if i > 0 else 0)
        resid_timings.append(idle_end - idle_start if i > 0 else 0)

        util.write_to_csv(resid_timings, breakdown_file)

    comm.Barrier()
    prev_estimates = comm.allgather(estimate.pixels.data[0,0,:,:])

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    master()
else:
    recon(rank - 1)