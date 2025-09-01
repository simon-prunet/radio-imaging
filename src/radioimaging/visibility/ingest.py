#!/usr/bin/env python

"""\
This file contains code that reads from a measurement set.
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import pandas
import numpy

from astropy.coordinates import SkyCoord, EarthLocation
from astropy.units import Quantity
from astropy import units as u

from ska_sdp_func_python.imaging import create_image_from_visibility
from ska_sdp_datamodels.configuration.config_model import Configuration
from ska_sdp_datamodels.science_data_model.polarisation_model import (
    ReceptorFrame,
    PolarisationFrame,
)
from ska_sdp_datamodels.visibility.vis_model import Visibility
from ska_sdp_datamodels.image.image_model import Image
from ska_sdp_func_python.visibility import convert_visibility_to_stokesI

def generate_fake_baselines(nbaselines):
    """
    generate_fake_baselines is a helper function for create_visibility_from_ms to help deal with the flattened visibility structure

    :param nbaselines: number of fake baselines to generated
    :return: List of fake baselines
    """
    for i in range(0, nbaselines):
        yield 0, 0


def create_visibility_from_ms(
    msname,
    channum=None,
    start_chan=None,
    end_chan=None,
    ack=False,
    datacolumn="DATA",
    selected_sources=None,
    selected_dds=None,
    average_channels=False,
    use_weight_spec=True,
):
    """
    create_visibility_from_ms is an updated version of RASCIL's create_visibility_from_ms function, made to more efficient, both memory and processing-wise
    Flattens the visibility structure to [1, time x baselines, freq, pol]. For this reason, this set of visibilities shouldn't be used for baseline or time
    specific operations such as BDA, and rather only for computing residuals

    :msname: File name of MS
    :channum: range of channels e.g. range(17,32), default is None meaning all
    :start_chan: Starting channel to read
    :end_chan: End channel to read
    :ack: Ask casacore to acknowledge each table operation
    :datacolumn: MS data column to read DATA, CORRECTED_DATA, or MODEL_DATA
    :selected_sources: Sources to select
    :selected_dds: Data descriptors to select
    :average_channels: Average all channels read
    :use_weight_spec: Use WEIGHT_SPECTRUM column instead of WEIGHT column for visibility weights
    :return: RASCIL visibility object
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
                        ms_weight = ms.getcol("WEIGHT_SPECTRUM") if use_weight_spec else ms.getcol("WEIGHT")

                    except IndexError:
                        raise IndexError("channel number exceeds max. within ms")

                else:
                    try:
                        channum = range(channels)
                        ms_vis = ms.getcol(datacolumn)[:, channum, :]
                        ms_weight = ms.getcol("WEIGHT_SPECTRUM")[:, channum, :] if use_weight_spec else ms.getcol("WEIGHT")
                        ms_flags = ms.getcol("FLAG")[:, channum, :]
                        channum = range(channels)
                    except IndexError:
                        raise IndexError("channel number exceeds max. within ms")
            else:
                channum = range(channels)
                try:
                    ms_vis = ms.getcol(datacolumn)[:, channum, :]
                    ms_flags = ms.getcol("FLAG")[:, channum, :]
                    ms_weight = ms.getcol("WEIGHT_SPECTRUM")[:, channum, :] if use_weight_spec else ms.getcol("WEIGHT")[:, :]
                except IndexError:
                    raise IndexError("channel number exceeds max. within ms")

            if average_channels:
                weight = ms_weight * (1.0 - ms_flags) if use_weight_spec else ms_weight[:, numpy.newaxis, :] * (1.0 - ms_flags)
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
                generate_fake_baselines(ms_vis.shape[0]), names=("antenna1", "antenna2")
                #generate_baselines(nants), names=("antenna1", "antenna2")
                #generate_baselines(1), names=("antenna1", "antenna2")
            )

            #nbaselines = len(baselines)
            nbaselines = ms_vis.shape[0]

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

            #ntimes = ms_vis.shape[0]
            ntimes = 1

            bv_times = numpy.zeros([ntimes])
            bv_vis = numpy.zeros([ntimes, nbaselines, nchan, npol]).astype("complex")
            bv_flags = numpy.zeros([ntimes, nbaselines, nchan, npol]).astype("int")
            bv_weight = numpy.zeros([ntimes, nbaselines, nchan, npol])
            bv_uvw = numpy.zeros([ntimes, nbaselines, 3])
            bv_integration_time = numpy.zeros([ntimes])

            for row, _ in enumerate(time):
                ibaseline = row

                time_index = 0

                bv_times[time_index] = time[row]
                bv_vis[time_index, ibaseline, ...] = ms_vis[row, ...]
                bv_flags[time_index, ibaseline, ...][
                    ms_flags[row, ...].astype("bool")
                ] = 1
                bv_weight[time_index, ibaseline, :, ...] = ms_weight[row, :, ...] if use_weight_spec else ms_weight[row, numpy.newaxis, ...]
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


def create_image_from_ms(ms_name, npixels, cellsize):
    """
    create_image_from_ms creates a dirty image from a given measurement set

    :msname: File name of MS
    :npixels: image dimensions, the same is used for both x and y axis
    :cellsize: angular resolution of each pixel in radians
    :return: dirty image as a RASCIL Image
    """
    [vis], _ = create_visibility_from_ms(ms_name, start_chan=0, end_chan=0, selected_dds=[0])
    vis = convert_visibility_to_stokesI(vis)

    return create_image_from_visibility(vis, cellsize=cellsize, npixel=npixels, polarisation_frame=vis.visibility_acc.polarisation_frame)