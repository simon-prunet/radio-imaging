#helper file for interleaved prototype. This may contain some duplicated code from the various python scripts in ../notebooks_and_pyscripts. This is done primarily to have a localized standalone prototype

import astropy
from astropy.io import fits

from ska_sdp_func_python.visibility import subtract_visibility
from ska_sdp_func_python.imaging import invert_ng, predict_ng, create_image_from_visibility, advise_wide_field
from ska_sdp_datamodels.gridded_visibility import create_griddata_from_image
from ska_sdp_func_python.grid_data import grid_visibility_weight_to_griddata,griddata_visibility_reweight, griddata_merge_weights
from rascil.processing_components import create_visibility_from_ms, generate_baselines
from ska_sdp_func_python.visibility import convert_visibility_to_stokesI
from ska_sdp_datamodels.configuration.config_model import Configuration
from ska_sdp_datamodels.science_data_model.polarisation_model import (
    ReceptorFrame,
    PolarisationFrame,
)
from ska_sdp_datamodels.visibility.vis_model import Visibility
import pandas
import numpy
import os
import time
import gc
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.units import Quantity
from astropy import units as u
import csv

def write_to_csv(data, filename):
    with open(filename, 'a+', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(data)

def tofits(data, filename):
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(filename, overwrite=True)
    hdulist.close()

def fromfits(filename):
    dat = fits.open(filename)[0].data
    while len(dat.shape) > 2:
        dat = dat[0]

    return dat

def compute_snr(gt, recon):
    difnorm = numpy.linalg.norm(gt-recon)
    if difnorm == 0:
        return 0
        
    return 20 * numpy.log10(numpy.linalg.norm(gt) / difnorm)

def compute_rmse(gt, recon):
    return numpy.sqrt(numpy.mean((gt-recon) ** 2))

def compute_windowed_var(image, window):
    estimated_variance = numpy.zeros(image.shape)

    #convolving initial signal with a variance estimation kernel with a size 2xwindow_hsize+1
    for y in range(image.shape[1]):
        for x in range(image.shape[0]):
            start_x = max(x - window, 0)
            end_x = min(x + window + 1, image.shape[0])
            start_y = max(y - window, 0)
            end_y = min(y + window + 1, image.shape[1])

            estimated_variance[x, y] = numpy.var(image[start_x:end_x, start_y:end_y])
            
    return estimated_variance

#copy of rascil's create_visibility_from_ms function, but with a fix to the channel allocation so that it works with large datasets
#also handy as the latest version of rascil has removed this function
def create_visibility_from_ms2(
    msname,
    channum=None,
    start_chan=None,
    end_chan=None,
    ack=False,
    datacolumn="DATA",
    selected_sources=None,
    selected_dds=None,
    average_channels=False,
):
    """Minimal MS to Visibility converter

    The MS format is much more general than the RASCIL Visibility so we cut many corners.
    This requires casacore to be installed. If not an exception ModuleNotFoundError is raised.

    Creates a list of Visibility's, split by field and spectral window

    Reading of a subset of channels is possible using either start_chan and end_chan or channnum. Using start_chan
    and end_chan is preferred since it only reads the channels required. Channum is more flexible and can be used to
    read a random list of channels.

    :param msname: File name of MS
    :param channum: range of channels e.g. range(17,32), default is None meaning all
    :param start_chan: Starting channel to read
    :param end_chan: End channel to read
    :param ack: Ask casacore to acknowledge each table operation
    :param datacolumn: MS data column to read DATA, CORRECTED_DATA, or MODEL_DATA
    :param selected_sources: Sources to select
    :param selected_dds: Data descriptors to select
    :param average_channels: Average all channels read
    :return: List of Visibility

    For example::

        selected_sources = ['1302+5748', '1252+5634']
        bvis_list = create_visibility_from_ms('../../data/3C277.1_avg.ms', datacolumn='CORRECTED_DATA',
                                           selected_sources=selected_sources)
        sources = numpy.unique([bv.source for bv in bvis_list])
        print(sources)
        ['1252+5634' '1302+5748']

    """
    try:
        from casacore.tables import table  # pylint: disable=import-error
    except ModuleNotFoundError:
        raise ModuleNotFoundError("casacore is not installed")
    try:
        from rascil.processing_components.visibility import msv2
    except ModuleNotFoundError:
        raise ModuleNotFoundError("cannot import msv2")

    tab = table(msname, ack=ack)

    if selected_sources is None:
        fields = numpy.unique(tab.getcol("FIELD_ID"))
    else:
        fieldtab = table("%s/FIELD" % msname, ack=False)
        sources = fieldtab.getcol("NAME")
        fields = list()
        for field, source in enumerate(sources):
            if source in selected_sources:
                fields.append(field)
        assert len(fields) > 0, "No sources selected"

    if selected_dds is None:
        dds = numpy.unique(tab.getcol("DATA_DESC_ID"))
    else:
        dds = selected_dds

    total_vis = 0

    vis_list = list()
    for field in fields:
        ftab = table(msname, ack=ack).query("FIELD_ID==%d" % field, style="")
        assert ftab.nrows() > 0, "Empty selection for FIELD_ID=%d" % (field)
        for dd in dds:
            # Now get info from the subtables
            ddtab = table("%s/DATA_DESCRIPTION" % msname, ack=False)
            spwid = ddtab.getcol("SPECTRAL_WINDOW_ID")[dd]
            polid = ddtab.getcol("POLARIZATION_ID")[dd]
            ddtab.close()

            meta = {"MSV2": {"FIELD_ID": field, "DATA_DESC_ID": dd}}
            ms = ftab.query("DATA_DESC_ID==%d" % dd, style="")
            assert (
                ms.nrows() > 0
            ), "Empty selection for FIELD_ID=%d and DATA_DESC_ID=%d" % (field, dd)
            # The TIME column has descriptor:
            # {'valueType': 'double', 'dataManagerType': 'IncrementalStMan', 'dataManagerGroup': 'TIME',
            # 'option': 0, 'maxlen': 0, 'comment': 'Modified Julian Day',
            # 'keywords': {'QuantumUnits': ['s'], 'MEASINFO': {'type': 'epoch', 'Ref': 'UTC'}}}
            otime = ms.getcol("TIME")
            datacol = ms.getcol(datacolumn, nrow=1)
            datacol_shape = list(datacol.shape)
            channels = datacol.shape[-2]
            if channum is None:
                if start_chan is not None and end_chan is not None:
                    try:
                        blc = [start_chan, 0]
                        trc = [end_chan, datacol_shape[-1] - 1]
                        channum = range(start_chan, end_chan + 1)
                        ms_vis = ms.getcolslice(datacolumn, blc=blc, trc=trc)
                        ms_flags = ms.getcolslice("FLAG", blc=blc, trc=trc)
                        ms_weight = ms.getcol("WEIGHT")

                    except IndexError:
                        raise IndexError("channel number exceeds max. within ms")

                else:
                    try:
                        channum = range(channels)
                        ms_vis = ms.getcol(datacolumn)[:, channum, :]
                        ms_weight = ms.getcol("WEIGHT")
                        ms_flags = ms.getcol("FLAG")[:, channum, :]
                        channum = range(channels)
                    except IndexError:
                        raise IndexError("channel number exceeds max. within ms")
            else:
                channum = range(channels)
                try:
                    ms_vis = ms.getcol(datacolumn)[:, channum, :]
                    ms_flags = ms.getcol("FLAG")[:, channum, :]
                    ms_weight = ms.getcol("WEIGHT")[:, :]
                except IndexError:
                    raise IndexError("channel number exceeds max. within ms")

            if average_channels:
                weight = ms_weight[:, numpy.newaxis, :] * (1.0 - ms_flags)
                ms_vis = numpy.sum(weight * ms_vis, axis=-2)[..., numpy.newaxis, :]
                sumwt = numpy.sum(weight, axis=-2)[..., numpy.newaxis, :]
                ms_vis[sumwt > 0.0] = ms_vis[sumwt > 0] / sumwt[sumwt > 0.0]
                ms_vis[sumwt <= 0.0] = 0.0 + 0.0j
                ms_flags = sumwt
                ms_flags[ms_flags <= 0.0] = 1.0
                ms_flags[ms_flags > 0.0] = 0.0

            total_vis += ms_vis.shape[0]

            uvw = -1 * ms.getcol("UVW")
            antenna1 = ms.getcol("ANTENNA1")
            antenna2 = ms.getcol("ANTENNA2")
            integration_time = ms.getcol("INTERVAL")

            time = otime - integration_time / 2.0

            start_time = numpy.min(time) / 86400.0
            end_time = numpy.max(time) / 86400.0

            spwtab = table("%s/SPECTRAL_WINDOW" % msname, ack=False)
            cfrequency = numpy.array(spwtab.getcol("CHAN_FREQ")[spwid][channum])
            cchannel_bandwidth = numpy.array(
                spwtab.getcol("CHAN_WIDTH")[spwid][channum]
            )
            nchan = cfrequency.shape[0]
            if average_channels:
                cfrequency = numpy.array([numpy.average(cfrequency)])
                cchannel_bandwidth = numpy.array([numpy.sum(cchannel_bandwidth)])
                nchan = cfrequency.shape[0]
            else:
                nchan = len(channum)

            # Get polarisation info
            poltab = table("%s/POLARIZATION" % msname, ack=False)
            corr_type = poltab.getcol("CORR_TYPE")[polid]
            corr_type = sorted(corr_type)
            # These correspond to the CASA Stokes enumerations
            if numpy.array_equal(corr_type, [1, 2, 3, 4]):
                polarisation_frame = PolarisationFrame("stokesIQUV")
                npol = 4
            elif numpy.array_equal(corr_type, [1, 2]):
                polarisation_frame = PolarisationFrame("stokesIQ")
                npol = 2
            elif numpy.array_equal(corr_type, [1, 4]):
                polarisation_frame = PolarisationFrame("stokesIV")
                npol = 2
            elif numpy.array_equal(corr_type, [5, 6, 7, 8]):
                polarisation_frame = PolarisationFrame("circular")
                npol = 4
            elif numpy.array_equal(corr_type, [5, 8]):
                polarisation_frame = PolarisationFrame("circularnp")
                npol = 2
            elif numpy.array_equal(corr_type, [9, 10, 11, 12]):
                polarisation_frame = PolarisationFrame("linear")
                npol = 4
            elif numpy.array_equal(corr_type, [9, 12]):
                polarisation_frame = PolarisationFrame("linearnp")
                npol = 2
            elif numpy.array_equal(corr_type, [9]) or numpy.array_equal(corr_type, [1]):
                npol = 1
                polarisation_frame = PolarisationFrame("stokesI")
            else:
                raise KeyError("Polarisation not understood: %s" % str(corr_type))

            # Get configuration
            anttab = table("%s/ANTENNA" % msname, ack=False)
            names = numpy.array(anttab.getcol("NAME"))

            ant_map = list()
            actual = 0
            # This assumes that the names are actually filled in!
            for i, name in enumerate(names):
                if name != "":
                    ant_map.append(actual)
                    actual += 1
                else:
                    ant_map.append(-1)
            # assert actual > 0, "Dish/station names are all blank - cannot load"
            if actual == 0:
                ant_map = list(range(len(names)))
                names = numpy.repeat("No name", len(names))

            mount = numpy.array(anttab.getcol("MOUNT"))[names != ""]
            diameter = numpy.array(anttab.getcol("DISH_DIAMETER"))[names != ""]
            xyz = numpy.array(anttab.getcol("POSITION"))[names != ""]
            offset = numpy.array(anttab.getcol("OFFSET"))[names != ""]
            stations = numpy.array(anttab.getcol("STATION"))[names != ""]
            names = numpy.array(anttab.getcol("NAME"))[names != ""]
            nants = len(names)

            antenna1 = list(map(lambda i: ant_map[i], antenna1))
            antenna2 = list(map(lambda i: ant_map[i], antenna2))

            baselines = pandas.MultiIndex.from_tuples(
                generate_baselines(nants), names=("antenna1", "antenna2")
            )
            nbaselines = len(baselines)

            location = EarthLocation(
                x=Quantity(xyz[0][0], "m"),
                y=Quantity(xyz[0][1], "m"),
                z=Quantity(xyz[0][2], "m"),
            )

            configuration = Configuration.constructor(
                name="",
                location=location,
                names=names,
                xyz=xyz,
                mount=mount,
                frame="ITRF",
                receptor_frame=ReceptorFrame("linear"),
                diameter=diameter,
                offset=offset,
                stations=stations,
            )
            # Get phasecentres
            fieldtab = table("%s/FIELD" % msname, ack=False)
            pc = fieldtab.getcol("PHASE_DIR")[field, 0, :]
            source = fieldtab.getcol("NAME")[field]
            phasecentre = SkyCoord(
                ra=pc[0] * u.rad, dec=pc[1] * u.rad, frame="icrs", equinox="J2000"
            )

            time_index_row = numpy.zeros_like(time, dtype="int")
            time_last = time[0]
            time_index = 0
            for row, _ in enumerate(time):
                if time[row] > time_last + 0.5 * integration_time[row]:
                    assert (
                        time[row] > time_last
                    ), "MS is not time-sorted - cannot convert"
                    time_index += 1
                    time_last = time[row]
                time_index_row[row] = time_index

            ntimes = time_index + 1

            assert ntimes == len(
                numpy.unique(time_index_row)
            ), "Error in finding data times"

            bv_times = numpy.zeros([ntimes])
            bv_vis = numpy.zeros([ntimes, nbaselines, nchan, npol]).astype("complex")
            bv_flags = numpy.zeros([ntimes, nbaselines, nchan, npol]).astype("int")
            bv_weight = numpy.zeros([ntimes, nbaselines, nchan, npol])
            bv_uvw = numpy.zeros([ntimes, nbaselines, 3])
            bv_integration_time = numpy.zeros([ntimes])

            for row, _ in enumerate(time):
                ibaseline = baselines.get_loc((antenna1[row], antenna2[row]))
                time_index = time_index_row[row]
                bv_times[time_index] = time[row]
                bv_vis[time_index, ibaseline, ...] = ms_vis[row, ...]
                bv_flags[time_index, ibaseline, ...][
                    ms_flags[row, ...].astype("bool")
                ] = 1
                bv_weight[time_index, ibaseline, :, ...] = ms_weight[
                    row, numpy.newaxis, ...
                ]
                bv_uvw[time_index, ibaseline, :] = uvw[row, :]
                bv_integration_time[time_index] = integration_time[row]

            vis_list.append(
                Visibility.constructor(
                    uvw=bv_uvw,
                    baselines=baselines,
                    time=bv_times,
                    frequency=cfrequency,
                    channel_bandwidth=cchannel_bandwidth,
                    vis=bv_vis,
                    flags=bv_flags,
                    weight=bv_weight,
                    integration_time=bv_integration_time,
                    configuration=configuration,
                    phasecentre=phasecentre,
                    polarisation_frame=polarisation_frame,
                    source=source,
                    meta=meta,
                )
            )
        tab.close()

    return vis_list, total_vis

def compute_residual(sky_estimate, vis, npixel, cellsize):
    #degridding, get visibilities of sky estimate
    vest = vis.copy(deep=True)
    vest = predict_ng(vest, sky_estimate, context='ng')

    #subtraction to obtain residual vis
    vres = subtract_visibility(vis, vest)

    #obtain dirty image and psf
    model = create_image_from_visibility(vres,cellsize=cellsize,npixel=npixel, polarisation_frame=vres.visibility_acc.polarisation_frame)

    dirty, sumwt = invert_ng(vres, model, context='ng')

    return dirty

def create_image_from_ms(ms_name, npixels, cellsize):
    [vis], _ = create_visibility_from_ms2(ms_name, start_chan=0, end_chan=0, selected_dds=[0])
    vis = convert_visibility_to_stokesI(vis)

    return create_image_from_visibility(vis, cellsize=cellsize, npixel=npixels, polarisation_frame=vis.visibility_acc.polarisation_frame)


#computes residual piecemeal channel by channel, this is so that we can handle very large datasets and not be bound by memory
#assumes that all channels are treated together, and that we only deal with Stokes I polarization
def compute_residual_bychannel(sky_estimate, ms_name, npixel, cellsize, weighting, robustness, weight_grid, channel_start, channel_end, data_descriptors, algorithm='ng'):
    final_residual = None

    read_from_disk_timings = []
    convert_polarization_timings = []
    weight_timings = []
    malloc_timings = []
    predict_timings = []
    subtract_timings = []
    invert_timings = []
    add_timings = []

    channels = range(channel_start, channel_end + 1)

    for dd in data_descriptors:
        for curr_channel in channels:
            read_start = time.time()
            [measured_vis], _ = create_visibility_from_ms2(ms_name, start_chan=curr_channel, end_chan=curr_channel, selected_dds=[dd])
            polarization_start = time.time()
            measured_vis = convert_visibility_to_stokesI(measured_vis)
            weight_start = time.time()
            measured_vis = griddata_visibility_reweight(measured_vis, weight_grid[0], weighting=weighting, robustness=robustness, sumwt=weight_grid[1])

            allocate_start = time.time()
            estimated_vis = measured_vis.copy(deep=True)
            predict_start = time.time()
            estimated_vis = predict_ng(estimated_vis, sky_estimate, context=algorithm)
            subtract_start = time.time()
            residual_vis = subtract_visibility(measured_vis, estimated_vis)

            if final_residual is None:
                final_residual = create_empty_image(measured_vis, npixel, cellsize)

            invert_start = time.time()
            channel_residual, sumwt = invert_ng(residual_vis, final_residual, context=algorithm)

            add_start = time.time()
            final_residual = add_to_image(final_residual, channel_residual.pixels.data[0,0,:,:])
            channel_end = time.time()

            gc.collect()

            read_from_disk_timings.append(polarization_start - read_start)
            convert_polarization_timings.append(weight_start - polarization_start)
            weight_timings.append(allocate_start - weight_start)
            malloc_timings.append(predict_start - allocate_start)
            predict_timings.append(subtract_start - predict_start)
            subtract_timings.append(invert_start - subtract_start)
            invert_timings.append(add_start - invert_start)
            add_timings.append(channel_end - add_start)

    

    read_from_disk_total = sum(read_from_disk_timings)
    convert_polarization_total = sum(convert_polarization_timings)
    weight_total = sum(weight_timings)
    malloc_total = sum(malloc_timings)
    predict_total = sum(predict_timings)
    subtract_total = sum(subtract_timings)
    invert_total = sum(invert_timings)
    add_total = sum(add_timings)

    gc.collect()

    return final_residual, [read_from_disk_total, convert_polarization_total, weight_total, malloc_total, predict_total, subtract_total, invert_total, add_total]


def compute_psf(vis, npixel, cellsize):
    model = create_image_from_visibility(vis,cellsize=cellsize,npixel=npixel, polarisation_frame=vis.visibility_acc.polarisation_frame)
    psf, sumwt = invert_ng(vis, model, context='ng', dopsf=True)

    return psf

#computes psf piecemeal channel by channel, this is so that we can handle very large datasets and not be bound by memory
#assumes that all channels are treated together, and that we only deal with Stokes I polarization
def compute_psf_by_channel(ms_name, npixel, cellsize, weight_grid, weighting, robustness, channel_start, channel_end, data_descriptors, algorithm='ng'):
    model = None
    psf = None

    read_timings = []
    polarization_timings = []
    weight_timings = []
    invert_timings = []
    add_timings = []

    channels = range(channel_start, channel_end + 1)

    for dd in data_descriptors:
        for curr_channel in channels:
            read_start = time.time()
            [vis], _ = create_visibility_from_ms2(ms_name, start_chan=curr_channel, end_chan=curr_channel, selected_dds=[dd])
            polar_start = time.time()
            vis = convert_visibility_to_stokesI(vis)
            weight_start = time.time()
            vis = griddata_visibility_reweight(vis, weight_grid[0], weighting=weighting, robustness=robustness, sumwt=weight_grid[1])

            if model is None:
                model = create_image_from_visibility(vis, cellsize=cellsize, npixel=npixel, polarisation_frame=vis.visibility_acc.polarisation_frame)

            invert_start = time.time()
            curr_psf, sumwt = invert_ng(vis, model, context=algorithm, dopsf=True)

            add_start = time.time()
            if psf is None:
                psf = curr_psf
            else:
                psf = add_to_image(psf, curr_psf.pixels.data[0,0,:,:])
            channel_end = time.time()

            gc.collect()

            read_timings.append(polar_start - read_start)
            polarization_timings.append(weight_start - polar_start)
            weight_timings.append(invert_start - weight_start)
            invert_timings.append(add_start - invert_start)
            add_timings.append(channel_end - add_start)


    read_total = sum(read_timings)
    polarization_total = sum(polarization_timings)
    weight_total = sum(weight_timings)
    invert_total = sum(invert_timings)
    add_total = sum(add_timings)

    gc.collect()

    return psf, model, [read_total, polarization_total, weight_total, invert_total, add_total]

def compute_weights(vis, npixel, cellsize, weighting, robustness=0.0):
    if (weighting != "natural"):
        model = create_image_from_visibility(vis, cellsize=cellsize, npixel=npixel, polarisation_frame=vis.visibility_acc.polarisation_frame)
        grid_weights = create_griddata_from_image(model, polarisation_frame=model.image_acc.polarisation_frame)
        grid_weights = grid_visibility_weight_to_griddata(vis, grid_weights)
        vis = griddata_visibility_reweight(vis, grid_weights[0], weighting=weighting, robustness=robustness, sumwt=grid_weights[1])
    else:
        vis = griddata_visibility_reweight(vis, None, weighting=weighting)

    return vis

#create griddata piecemeal by channel. This is essentially a tally per grid cell for all the visibilities that fall into it, and is needed for computing visibility weights
#piecemeal, as it cannot be done in one go for large datasets
def compute_weights_griddata_by_channel(ms_name, npixel, cellsize, channel_start, channel_end, data_descriptors):
    total_grid = None
    model = None

    read_timings = []
    polarization_timings = []
    grid_timings = []

    channels = range(channel_start, channel_end + 1)

    total_vis = 0

    for dd in data_descriptors:
        for curr_channel in channels:
            read_start = time.time()
            [vis], num_vis = create_visibility_from_ms2(ms_name, start_chan=curr_channel, end_chan=curr_channel, selected_dds=[dd])
            total_vis += num_vis

            pol_start = time.time()
            vis = convert_visibility_to_stokesI(vis)

            if model is None:
                model = create_image_from_visibility(vis, cellsize=cellsize, npixel=npixel, polarisation_frame=vis.visibility_acc.polarisation_frame)

            grid_start = time.time()
            curr_grid_weights = create_griddata_from_image(model, polarisation_frame=model.image_acc.polarisation_frame)
            curr_grid_weights = grid_visibility_weight_to_griddata(vis, curr_grid_weights)
            
            if total_grid is None:
                total_grid = curr_grid_weights
            else:
                total_grid = griddata_merge_weights([total_grid, curr_grid_weights])

            channel_end = time.time()

            gc.collect()

            read_timings.append(pol_start - read_start)
            polarization_timings.append(grid_start - pol_start)
            grid_timings.append(channel_end - grid_start)
            print(str(dd) + ": " + str(curr_channel))

    read_total = sum(read_timings)
    polarization_total = sum(polarization_timings)
    weight_total = sum(polarization_timings)

    gc.collect()

    return total_grid, [read_total, polarization_total, weight_total], total_vis


def create_empty_image(vis, npixel, cellsize):
    return create_image_from_visibility(vis, npixel=npixel, cellsize=cellsize, polarisation_frame=vis.visibility_acc.polarisation_frame)

def deconvolve_single(dirty, psf, niter, wavelet_type_idx, curr_maj_iter, initial_lambda, lambda_mul):
    res = numpy.array(dirty)
    np_psf = numpy.array(psf)

    curr_lambda = initial_lambda * (lambda_mul ** curr_maj_iter)
    curr_lambda *= numpy.linalg.norm(dirty)

    tmp_psf_name = "tmp_psf.fits"
    tmp_res_name = "tmp_residual.fits"
    tmp_output_name = "tmp_output.fits"

    tofits(psf, tmp_psf_name)
    tofits(dirty, tmp_res_name)

    os.system("julia julia/make_fullres.jl " + str(curr_lambda) + " " + tmp_psf_name + " " + tmp_res_name + " " + str(wavelet_type_idx) + " " + str(niter) + " " + tmp_output_name)

    deconvolved = fromfits(tmp_output_name)

    return deconvolved

#interleaved deconvolution, currently assumes 2 partitions
def deconvolve(step, dirty, psf, prev_estimates, niter, wavelet_type_idx, curr_maj_iter, initial_lambda, lambda_mul, cut_center, cut_halfwidth, variance_window, recon_variance_factor):
    res = numpy.array(dirty)
    np_psf = numpy.array(psf)

    curr_lambda = 0

    if curr_maj_iter > 1:
        curr_lambda = initial_lambda * (lambda_mul ** (curr_maj_iter - 1))

    tmp_psf_name = "tmp_psf_" + str(step) + ".fits"
    tmp_res_name = "tmp_residual_" + str(step) + ".fits"
    tmp_constraint_name = "tmp_constraint_" + str(step) + ".fits"
    tmp_output_name = "tmp_output_" + str(step) + ".fits"

    vis_variance = numpy.mean(compute_windowed_var(dirty, variance_window))

    low_variance = vis_variance if step == 0 else vis_variance / recon_variance_factor
    high_variance = vis_variance / recon_variance_factor if step == 0 else vis_variance

    tofits(psf, tmp_psf_name)
    tofits(dirty, tmp_res_name)

    constraint = prev_estimates[1] - prev_estimates[0] if step == 0 else prev_estimates[0] - prev_estimates[1]

    curr_lambda *= (numpy.linalg.norm(dirty) + numpy.linalg.norm(constraint))

    tofits(constraint, tmp_constraint_name)

    os.system("julia julia/make_multistep_interleaved.jl " + str(curr_lambda) + " " + tmp_psf_name + " " + tmp_res_name + " " + tmp_constraint_name + " " + str(wavelet_type_idx) + " " + str(niter) + " " + \
            str(low_variance) + " " + str(high_variance) + " " + str(cut_center) + " " + str(cut_halfwidth)  + " " + str(step) + " " + str(curr_maj_iter) + " " + tmp_output_name)

    deconvolved = fromfits(tmp_output_name)

    return deconvolved


#stub for now, implemented later
def add_to_image(image, nparr):
    image.pixels.data[0,0,:,:] += nparr

    return image
