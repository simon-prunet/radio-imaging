from mpi4py import MPI
import json
import ip_helpers as iph
import numpy
import time
import sys
from pathlib import Path
from ska_sdp_func_python.image.cleaners import msclean
import gc

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

    output_dir = config["output_dir"] + "_pmsc/"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    timings_file = output_dir + "mc_timings"

    psfs = comm.allgather(None)
    weights = comm.allgather(None)


    comm.Barrier()

    startup_end = time.time()
    iph.write_to_csv([startup_end - recon_start], timings_file)

    for i in range(config["nmajcycmsc"]):
        mc_start = time.time()
        comm.Barrier()

        all_estimates = comm.allgather(None)

        recons_combined = combine(all_estimates[1:], config["npixels"], average=(i > 0))

        iph.tofits(recons_combined, output_dir + "recon_combined_" + str(i) + ".fits")

        for j, curr_recon in enumerate(all_estimates[1:]):
            iph.tofits(curr_recon, output_dir + "recon_" + str(j) + "_" + str(i) + ".fits")

        mc_end = time.time()

        iph.write_to_csv([mc_end - mc_start], timings_file)

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
    data_descriptors = config["data_descriptors"]

    output_dir = config["output_dir"] + "_pmsc/"
    
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

    breakdown_file = output_dir + "mc_timings_breakdown_" + str(partition)

    weight_grid, weight_timings, num_vis = iph.compute_weights_griddata_by_channel(ms_name, npixels, cellsize, channel_start, channel_end, data_descriptors)
    psf, estimate, psf_timings, weight = iph.compute_psf_by_channel(ms_name, npixels, cellsize, weight_grid, weighting, robustness, channel_start, channel_end, data_descriptors)
    iph.tofits(numpy.fft.ifftshift(numpy.real(numpy.fft.fft2(psf.pixels.data[0,0,:,:]))), output_dir + "psf_" + str(partition) + ".fits")

    sigma2s = [1] * (len(ells) + 1)

    filters, ffilters = iph.create_nfilters(npixels, delta, ells, sigma2s)

    weights = numpy.zeros(len(sigma2s))

    #share psfs, needed for all major-cycles after the first
    psf_send_start = time.time()
    psfs = comm.allgather(psf.pixels.data[0,0,:,:])
    psfs = psfs[1:]
    weights = comm.allgather(weight)
    weights = weights[1:]
    psf_send_end = time.time()

    total_weight = sum(weights)

    weights = [weight / total_weight for weight in weights]

    joint_psf = iph.convolve2d(psfs[0] * weights[0], filters[0])
    for i, currpsf in enumerate(psfs[1:]):
        joint_psf += iph.convolve2d(currpsf * weights[i+1], filters[i+1])

    iph.tofits(joint_psf, output_dir + "joint_psf" + str(partition) + ".fits")

    barrier_start = time.time()
    comm.Barrier()
    barrier_end = time.time()

    iph.write_to_csv([num_vis], breakdown_file)
    iph.write_to_csv(weight_timings, breakdown_file)
    iph.write_to_csv(psf_timings, breakdown_file)
    iph.write_to_csv([psf_send_end - psf_send_start], breakdown_file)
    iph.write_to_csv([barrier_end - barrier_start], breakdown_file)

    smallest_scale = iph.find_clsize(filters[partition])
    # smallest_scale_idx = 0
    # while smallest_scale > scales[smallest_scale_idx] and smallest_scale_idx < (len(scales) - 1):
    #     smallest_scale_idx += 1
    # first_mc_scales = scales[smallest_scale_idx:]

    first_mc_scales = config["first_mc_clean_scales"][partition]

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

        residual, resid_timings = iph.compute_residual_bychannel(estimate, ms_name, npixels, cellsize, weighting, robustness, weight_grid, channel_start, channel_end, data_descriptors)

        iph.tofits(residual.pixels.data[0,0,:,:], output_dir + "residual_" + str(partition) + "_" + str(i) + ".fits")

        deconv_start = time.time()

        curr_residual = curr_psf = None

        if i > 0:
            curr_residual = iph.convolve2d(weights[partition] * residual["pixels"].data[0, 0, :, :], filters[partition])
            for j, sigma2 in enumerate(sigma2s):
                if j != partition:
                    constraint = prev_estimates[j] - prev_estimates[partition]
                    constraint = iph.convolve2d(constraint, weights[j] * psfs[j])
                    constraint = iph.convolve2d(constraint, filters[j])

                    curr_residual += constraint

            curr_psf = joint_psf
        else:
            curr_residual = residual.pixels.data[0,0,:,:]
            curr_psf = psf["pixels"].data[0, 0, :, :]

        iph.tofits(curr_residual, output_dir + "curr_residual_" + str(partition) + "_" + str(i) + ".fits")

        deconvolved, _ = msclean(curr_residual, curr_psf, None, None, gain, thresh, msc_niter, first_mc_scales if i == 0 else scales, fracthresh)

        deconv_end = time.time()

        iph.tofits(deconvolved, output_dir + "deconv_" + str(partition) + "_" + str(i) + ".fits")

        estimate = iph.add_to_image(estimate, deconvolved)

        resid_timings.append(deconv_end - deconv_start)
        resid_timings.append(sendrecon_end - sendrecon_start if i > 0 else 0)
        resid_timings.append(idle_end - idle_start if i > 0 else 0)

        iph.write_to_csv(resid_timings, breakdown_file)

    comm.Barrier()
    prev_estimates = comm.allgather(estimate.pixels.data[0,0,:,:])

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    master()
else:
    recon(rank - 1)