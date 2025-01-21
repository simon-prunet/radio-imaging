from mpi4py import MPI
import json
import ip_helpers as iph
from rascil.processing_components import create_visibility_from_ms
import numpy
import time
import sys
from pathlib import Path

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

    output_dir = config["output_dir"] + "_parallel/"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    timings_file = output_dir + "mc_timings"

    load_tokens = comm.gather(None, root=0)

    startup_end = time.time()
    iph.write_to_csv([startup_end - recon_start], timings_file)

    for i in range(config["nmajcyc"]):
        mc_start = time.time()
        #for now doing this once every major cycle, but in practice only needed at the end of the reconstruction
        reconstructions = comm.gather(None, root=0)

        recons_combined = combine(reconstructions, config["npixels"], average=(i > 0))

        iph.tofits(recons_combined, output_dir + "recon_combined_" + str(i) + ".fits")

        for j, curr_recon in enumerate(reconstructions):
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
    wavelet_idx = wavelet_type_dict[config["wavelet_dict"]]
    robustness = config["robustness"]
    channel_start = int(config["channel_start"])
    channel_end = int(config["channel_end"])
    data_descriptors = config["data_descriptors"]

    ells = config["ells"]
    ells.sort()

    delta = config["delta"]

    variance_window = config["visvar_window"]

    init_lambdas = config["init_lambdas"]
    init_lambda = init_lambdas[partition]

    lambda_muls = config["lambda_muls"]
    lambda_mul = lambda_muls[partition]

    n_fista_iter = config["nfistaiter"]
    
    output_dir = config["output_dir"] + "_parallel_multipartition/"

    breakdown_file = output_dir + "mc_timings_breakdown_" + str(partition)

    weight_grid, weight_timings, num_vis = iph.compute_weights_griddata_by_channel(ms_name, npixels, cellsize, channel_start, channel_end, data_descriptors)
    psf, estimate, psf_timings, weight = iph.compute_psf_by_channel(ms_name, npixels, cellsize, weight_grid, weighting, robustness, channel_start, channel_end, data_descriptors)
    iph.tofits(psf.pixels.data[0,0,:,:], output_dir + "psf_" + str(partition) + ".fits")
    other_estimate = iph.create_image_from_ms(ms_name, npixels, cellsize)

    #synchronization barrier primarily used to start the major cycles at the same time. This isn't strictly needed, but is handy for taking stats, as we can determine the amount
    #of waiting time before the first major cycle
    barrier_start = time.time()
    comm.Barrier()
    barrier_end = time.time()

    iph.write_to_csv([num_vis], breakdown_file)
    iph.write_to_csv(weight_timings, breakdown_file)
    iph.write_to_csv(psf_timings, breakdown_file)
    iph.write_to_csv([barrier_end - barrier_start], breakdown_file)

    for i in range(config["nmajcyc"]):
        send_start = time.time()
        #We use constraints after the first major cycle, which are injected into our objective function. These constraints are sent to and obtained from the other reconstruction node
        prev_estimates = None
        if i > 0:
            prev_estimates = comm.allgather(estimate.pixels.data[0,0,:,:])
        send_end = time.time()

        residual, resid_timings = iph.compute_residual_bychannel(estimate, ms_name, npixels, cellsize, weighting, robustness, weight_grid, channel_start, channel_end, data_descriptors)

        iph.tofits(residual.pixels.data[0,0,:,:], output_dir + "residual_" + str(partition) + "_" + str(i) + ".fits")

        deconv_start = time.time()
        deconvolved = iph.deconvolve_multipartition(partition, residual.pixels.data[0,0,:,:], psf.pixels.data[0,0,:,:], prev_estimates, n_fista_iter, wavelet_idx, i, init_lambda, lambda_mul, \
            ells, delta, variance_window)
        deconv_end = time.time()

        estimate = iph.add_to_image(estimate, deconvolved)

        sendrecon_start = time.time()
        comm.gather(deconvolved, root=0)
        sendrecon_end = time.time()

        resid_timings.insert(0, send_end - send_start)
        resid_timings.append(deconv_end - deconv_start)
        resid_timings.append(sendrecon_end - sendrecon_start)

        iph.write_to_csv(resid_timings, breakdown_file)


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    master()
else:
    recon(rank - 1)