from mpi4py import MPI
import json
import ip_helpers as iph
from rascil.processing_components import create_visibility_from_ms
import numpy

wavelet_type_dict = {"daubechies" : 0, "iuwt" : 1}

def combine(reconstructions, npixels):
    final_recon = numpy.zeros((npixels, npixels))

    for recon in reconstructions:
        while len(recon.shape) > 2:
            recon = recon[0]

        final_recon += recon

    final_recon /= len(reconstructions)

    return final_recon

#master node, responsible for reading config, sending this to individual nodes, and gathering and combining reconstructed images
def master():
    with open("setup.config") as f: 
        data = f.read() 

    config = json.loads(data)
    
    config = comm.bcast(config, root=0)

    for i in range(config["nmajcyc"]):
        #for now doing this once every major cycle, but in practice only needed at the end of the reconstruction
        reconstructions = comm.gather(None, root=0)
        curr_recon = combine(reconstructions[1:], config["npixels"])
        iph.tofits(curr_recon, config["output_filename"] + "_" + str(i) + ".fits")
        pass

#reconstruction nodes, responsible for reconstructing the individual full resolution images. For now assumes only 2 partitions
def recon(step):
    config = comm.bcast(None, root=0)

    msname = config["lowres-dataset"] if step == 0 else config["highres-dataset"]
    npixels = config["npixels"]
    cellsize = config["cellsize"]
    weighting = config["weighting"]
    wavelet_idx = wavelet_type_dict[config["wavelet_dict"]]

    [vis] = create_visibility_from_ms(msname)

    vis = iph.compute_weights(vis, npixels, cellsize, weighting)

    estimate = iph.create_empty_image(vis, npixels, cellsize)
    other_estimate = iph.create_empty_image(vis, npixels, cellsize)

    psf = iph.compute_psf(vis, npixels, cellsize)

    init_lambda = config["init_lambda_low"] if step == 0 else config["init_lambda_high"]
    lambda_mul = config["lambda_mul_low"] if step == 0 else config["lambda_mul_high"]

    for i in range(config["nmajcyc"]):
        print(i)
        #We use constraints after the first major cycle, which are injected into our objective function. These constraints are sent to and obtained from the other reconstruction node
        if i > 0:
            if step == 0:
                comm.Send(estimate.pixels.data, dest=2)
                comm.Recv(other_estimate.pixels.data, source=2)
            elif step == 1:
                comm.Recv(other_estimate.pixels.data, source=1)
                comm.Send(estimate.pixels.data, dest=1)

        residual = iph.compute_residual(estimate, vis, npixels, cellsize)

        constraint = [estimate.pixels.data[0,0,:,:], other_estimate.pixels.data[0,0,:,:]] if step == 0 else [other_estimate.pixels.data[0,0,:,:], estimate.pixels.data[0,0,:,:]]
        deconvolved = iph.deconvolve(step, residual.pixels.data[0,0,:,:], psf.pixels.data[0,0,:,:], constraint, config["nfistaiter"], wavelet_idx, i, init_lambda, lambda_mul, \
            config["sep_center"], config["sep_hw"], config["visvar_window"], config["reconvar_factor"])

        estimate = iph.add_to_image(estimate, deconvolved)
        comm.gather(estimate.pixels.data, root=0)


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    master()
elif rank == 1 or rank == 2:
    recon(rank - 1)