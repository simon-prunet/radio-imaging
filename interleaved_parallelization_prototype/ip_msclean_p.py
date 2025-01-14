from mpi4py import MPI
import json
import ip_helpers as iph
import numpy
import time
import sys
from pathlib import Path
from ska_sdp_func_python.image.cleaners import msclean

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

    output_dir = config["output_dir"] + "_msclean_p/"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    timings_file = output_dir + "mc_timings"

    startup_end = time.time()
    iph.write_to_csv([startup_end - recon_start], timings_file)

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

        iph.tofits(curr_recon, output_dir + "deconv_combined_" + str(i) + ".fits")
        iph.tofits(reconstructions[1], output_dir + "deconv_vl_" + str(i) + ".fits")
        iph.tofits(reconstructions[2], output_dir + "deconv_vh_" + str(i) + ".fits")

        mc_end = time.time()

        iph.write_to_csv([mc_end - mc_start], timings_file)

#reconstruction nodes, responsible for reconstructing the individual full resolution images. For now assumes only 2 partitions
def recon(step):
    config = comm.bcast(None, root=0)

    ms_name = config["lowres-dataset"] if step == 0 else config["highres-dataset"]
    npixels = config["npixels"]
    cellsize = config["cellsize"]
    weighting = config["weighting"]
    robustness = config["robustness"]
    channel_start = int(config["channel_start"])
    channel_end = int(config["channel_end"])
    wavelet_idx = wavelet_type_dict[config["wavelet_dict"]]
    data_descriptors = config["data_descriptors"]
    output_dir = config["output_dir"] + "_msclean_p/"
    msc_niter = int(config["msclean_iter"])
    delta = config["sep_center"]
    ell = config["sep_hw"]
    var_window = config["visvar_window"]
    thresh = config["clean_thresh"]
    scales = config["clean_scales"]
    small_scales = [0, 1, 2, 4, 6, 10]

    sens = None
    gain = 0.1
    fracthresh = 1e-3

    breakdown_file = output_dir + "mc_timings_breakdown_" + str(step)

    weight_grid, weight_timings, num_vis = iph.compute_weights_griddata_by_channel(ms_name, npixels, cellsize, channel_start, channel_end, data_descriptors)
    psf, estimate, psf_timings, weight = iph.compute_psf_by_channel(ms_name, npixels, cellsize, weight_grid, weighting, robustness, channel_start, channel_end, data_descriptors)
    iph.tofits(numpy.fft.ifftshift(numpy.real(numpy.fft.fft2(psf.pixels.data[0,0,:,:]))), output_dir + "psf_" + str(step) + ".fits")
    other_estimate = iph.create_image_from_ms(ms_name, npixels, cellsize)
    other_psf = iph.create_image_from_ms(ms_name, npixels, cellsize)
    local_filter = other_filter = None

    if step == 0:
        local_filter, other_filter = iph.create_filters(npixels, delta, ell, 1, 1)
    else:
        other_filter, local_filter = iph.create_filters(npixels, delta, ell, 1, 1)

    other_weight = numpy.zeros(1)

    print(weight)

    #share psfs, needed for all major-cycles after the first
    psf_send_start = time.time()
    if step == 0:
        comm.Send(psf.pixels.data, dest=2)
        comm.Recv(other_psf.pixels.data, source=2)
        comm.Send(weight, dest=2)
        comm.Recv(other_weight, source=2)
    elif step == 1:
        comm.Recv(other_psf.pixels.data, source=1)
        comm.Send(psf.pixels.data, dest=1)
        comm.Recv(other_weight, source=1)
        comm.Send(weight, dest=1)

    psf_send_end = time.time()

    total_weight = weight + other_weight
    corrected_local_weight = weight / total_weight
    corrected_other_weight = other_weight / total_weight

    joint_psf = corrected_local_weight * iph.convolve2d(psf["pixels"].data[0, 0, :, :], local_filter) + corrected_other_weight * iph.convolve2d(other_psf["pixels"].data[0, 0, :, :], other_filter)
    iph.tofits(joint_psf, output_dir + "joint_psf" + str(step) + ".fits")

    iph.write_to_csv([num_vis], breakdown_file)
    iph.write_to_csv(weight_timings, breakdown_file)
    iph.write_to_csv(psf_timings, breakdown_file)
    iph.write_to_csv([psf_send_end - psf_send_start], breakdown_file)

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

        residual, resid_timings = iph.compute_residual_bychannel(estimate, ms_name, npixels, cellsize, weighting, robustness, weight_grid, channel_start, channel_end, data_descriptors)

        iph.tofits(residual.pixels.data[0,0,:,:], output_dir + "residual_" + str(step) + "_" + str(i) + ".fits")

        deconv_start = time.time()

        curr_residual = curr_psf = None

        if i > 0:
            constraint = other_estimate.pixels.data[0, 0, :, :] - estimate.pixels.data[0, 0, :, :]
            other_residual = iph.convolve2d(constraint, other_psf["pixels"].data[0, 0, :, :])
            #local_resid_var = iph.compute_windowed_var(residual["pixels"].data[0, 0, :, :], var_window)
            #other_resid_var = iph.compute_windowed_var(other_residual, var_window)

            curr_residual = corrected_local_weight * residual["pixels"].data[0, 0, :, :] + corrected_other_weight * iph.convolve2d(other_residual, other_filter)
            curr_psf = joint_psf
        else:
            curr_residual = residual.pixels.data[0,0,:,:]
            curr_psf = psf["pixels"].data[0, 0, :, :]

        iph.tofits(curr_psf, output_dir + "curr_psf_" + str(step) + "_" + str(i) + ".fits")
        iph.tofits(curr_residual, output_dir + "curr_residual_" + str(step) + "_" + str(i) + ".fits")

        deconvolved, _ = msclean(curr_residual, curr_psf, None, sens, gain, thresh, msc_niter if i > 0 or step == 1 else msc_niter // 2, scales if i > 0 or step == 0 else small_scales, fracthresh)

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
elif rank == 1 or rank == 2:
    recon(rank - 1)