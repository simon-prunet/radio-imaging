#!/usr/bin/env python

"""\
Code for performing baseline dependent averaging, as described in Wijnholds et al 2018 but also taking into account frequency
decorrelation.
all ha/ra/dec values are assumed to be in radians (except for the initial call to bda, which is in degrees), frequencies in hertz (1/s)
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import astropy
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, Angle
from astropy import constants as const
import astropy.units as u
from astropy.time import Time

import matplotlib.pyplot as plt
import matplotlib
import numpy
import math
import pandas
import time

import csv

from scipy import signal
from scipy import stats

import os

from ska_sdp_datamodels.gridded_visibility import create_griddata_from_image
from ska_sdp_datamodels.science_data_model.polarisation_model import (
    ReceptorFrame,
    PolarisationFrame,
)
from ska_sdp_datamodels.visibility.vis_model import Visibility
from ska_sdp_datamodels.configuration.config_model import Configuration

import numpy

import casacore.tables.tableutil as pt
from casacore.tables import (
    makescacoldesc,
    makearrcoldesc,
    table,
    maketabdesc,
    tableexists,
    tableiswritable,
    tableinfo,
    tablefromascii,
    tabledelete,
    makecoldesc,
    msconcat,
    removeDerivedMSCal,
    taql,
    tablerename,
    tablecopy,
    tablecolumn,
    addDerivedMSCal,
    removeImagingColumns,
    addImagingColumns,
    required_ms_desc,
    tabledefinehypercolumn,
    default_ms,
    makedminfo,
    default_ms_subtable,
)
from rascil.processing_components.visibility.msv2fund import Antenna, Stand
from rascil.processing_components.visibility import msv2

from radioimaging.visibility import ingest

#earth rot speed, in rads per second
omega = 7.2925e-5

#PS. 

def create_config(vis):
    """
    create_config creates a pseudo config for the visibility constructor, 
    doesn't really matter what is put here since we never plan to save the visibilities

    :vis: RASCIL visibility structure
    :return: pseudo configuration
    """

    return Configuration.constructor(
        name=vis.configuration.name,
        location=vis.configuration.location,
        names=vis.configuration.names,
        xyz=vis.configuration.xyz,
        mount=vis.configuration.mount,
        frame=vis.configuration.frame,
        receptor_frame=vis.configuration.receptor_frame,
        diameter=vis.configuration.diameter,
        offset=vis.configuration.offset,
        stations=vis.configuration.stations,
    )

def get_baseline_coords(antenna_coords):
    """
    get_baseline_coords obtains the baseline coords for some given antenna coordinates. Assumes ECEF but euclidian
    frames also work since the coordinates are relative. Assumes all pairs of antennas are correlated.

    :antenna_coords: list of coordinates for antennas
    :return: list of baseline coordinates
    """

    nants = antenna_coords.shape[0]
    nbaselines = int(nants * (nants + 1) / 2)
    baseline_coords = numpy.zeros((nbaselines, 3))

    a1idx_accum = 0
    for a1 in range(nants):
        for a2 in range(a1, nants):
            a1coord = antenna_coords[a1,:]
            a2coord = antenna_coords[a2,:]

            baseline_coords[a1idx_accum + a2 - a1,:] = a1coord - a2coord

        a1idx_accum += (nants - a1)

    return baseline_coords

def deg2rad2d(coords):
    """
    deg2rad2d converts degrees to radians. Assumes 2d coordinates

    :coords: coordinates in degrees
    :return: coordinates in radians
    """
    return ((coords[0] * u.deg).to(u.rad).value, (coords[1] * u.deg).to(u.rad).value)

def dudt(lx, ly, ha, lambd):
    """
    dudt calculates dudt of visibility coordinates

    :lx: x baseline coordinate
    :ly: y baseline coordinate
    :ha: hour angle of observed object
    :lambd: freq/c
    :return: dudt
    """
    return (lx * numpy.cos(ha) - ly * numpy.sin(ha)) * omega / lambd

def dvdt(lx, ly, ha, dec, lambd):
    """
    dvdt calculates dvdt of visibility coordinates

    :lx: x baseline coordinate
    :ly: y baseline coordinate
    :ha: hour angle of observed object
    :dec: declination of observed object
    :lambd: freq/c
    :return: dvdt
    """
    return (lx * numpy.sin(dec) * numpy.sin(ha) + ly * numpy.sin(dec) * numpy.cos(ha)) * omega / lambd


def radec_to_ecefunit(ra, dec, t):
    """
    radec_to_ecefunit converts ra-dec to a unit vector in ecef

    :ra: right ascension
    :dec: declination
    :t: time of observation as an astropy time object
    :return: ecef unit vector pointing to ra-dec coordinates
    """
    gst = t.sidereal_time('apparent', 'greenwich')
    gst_rad = (gst.hour * 15 * u.deg).to(u.rad).value

    return (numpy.cos(dec) * numpy.cos(ra + gst_rad), numpy.cos(dec) * numpy.sin(ra + gst_rad), numpy.sin(dec))

def get_maxdecorr_obstime_source(phasecenter, half_fov, vis, telescope_coords):
    """
    get_maxdecorr_obstime_source gets the time of the observation that has the maximum decorrelation. This would be when the source is
    the lowest on the horizon, ie. largest absolute hour angle

    :phasecenter: pointing phase center
    :half_fov: half of the fov of the observed field
    :vis: RASCIL visibility structure containing observation
    :telescope_coords: telescope coordinates in ecef (typically center of array)
    :return: time and position corresponding to the maximum decorrelation timestamp
    """

    sources = [(phasecenter[0] + half_fov, phasecenter[1] + half_fov),
                (phasecenter[0] + half_fov, phasecenter[1] - half_fov),
                (phasecenter[0] - half_fov, phasecenter[1] + half_fov),
                (phasecenter[0] - half_fov, phasecenter[1] - half_fov)]

    maxha = 0
    maxtime = None
    source_idx = 0

    observer_location = EarthLocation(lon=vis.configuration.location.lon.value, lat=vis.configuration.location.lat.value)

    for datetime in vis.datetime:
        t = Time(datetime, scale='utc', location=observer_location)
        lst = t.sidereal_time('mean').to(u.hourangle)

        for i, src in enumerate(sources):
            src_ra = (src[0] * u.deg).to(u.hourangle)
            ha = numpy.abs((src_ra - lst).value)

            if ha > maxha:
                maxha = ha
                source_idx = i
                maxtime = t

    return t, sources[source_idx]

#freq decorrelation
def freq_decorr(channel_width, baseline_coords_ecef, source_radec, obs_datetime):
    """
    freq_decorr calculates frequency decorrelation

    :channel_width: bandwidth of channel
    :baseline_coords_ecef: baseline coordinates in ecef
    :source_radec: coordinates of source
    :obs_datetime: observation time
    :return: a value between 0 and 1 denoting frequency decorrelation, with 1 being no decorrelation
    """
    src_coords = radec_to_ecefunit(source_radec[0], source_radec[1], obs_datetime)
    tau_g = (src_coords[0] * baseline_coords_ecef[0] + src_coords[1] * baseline_coords_ecef[1] + src_coords[2] * baseline_coords_ecef[2]) / const.c

    #we don't multiply by pi since numpy's sinc does this already, 1 is no decorrelation
    return numpy.abs(numpy.sinc((channel_width * tau_g).value / 2))

def time_decorr(baseline_coord, total_int_time, source_radec, freq, phase_center, obs_datetime):
    """
    time_decorr estimates the maximum residual time decorrelation assuming perfect tracking. Uses the taylor approximation from 
    wijnholds 2018 paper eq 41

    :baseline_coord: baseline coordinates in ecef
    :total_int_time: total integration time after averaging
    :source_radec: coordinates of source
    :freq: em frequency
    :phase_center: phase center of pointing
    :obs_datetime: observation time
    :return: a value between 0 and 1 denoting residual time decorrelation, with 1 being no decorrelation, and a change rate
    variable which can be used to calculate maximum integration time given some maximum decorrelation
    """
    lambd = const.c / freq

    lst = obs_datetime.sidereal_time("mean")
    lst_rad = (lst.hour * 15 * u.deg).to(u.rad).value
    ha = lst_rad - source_radec[0]
    dec = source_radec[1]

    du = dudt(baseline_coord[0], baseline_coord[1], ha, lambd)
    dv = dvdt(baseline_coord[0], baseline_coord[1], ha, dec, lambd)

    #small angles approximation, this may need to change if we are simultaneously observing large fields, needed because ra converges towards celestial poles
    l = (source_radec[0] - phase_center[0]) * numpy.cos(phase_center[1])
    m = source_radec[1] - phase_center[1]

    change_rate = ((l * du + m * dv) ** 2).value

    # eq 41 of wijnholds paper
    decorr_frac = 1 - (numpy.pi**2 * total_int_time ** 2) * change_rate / 6.0 
    
    #1 is no decorrelation
    return decorr_frac, change_rate

#computes the number of subdivisions we can perform while maintaining a max decorrelation below what is specified
#assumes that the channel bandwidths and visibility integration times are all the same
def averaging_levels(max_decorr, baseline_coord, visint_time, freq, channel_width, source_radec, phase_center, obs_datetime, only_residual_decorr):
    """
    averaging_levels calculates the amount of averaging (by power of two) that one can perform given some maximum decorrelation

    :max_decorr: max decorrelation
    :baseline_coord: baseline coordinate in ecef
    :visint_time: pre-averaging (but post fs) visibility integration time
    :freq: em frequency
    :channel_width: channel bandwidth
    :source_radec: source coordinates
    :phase_center: pointing phase center
    :obs_datetime: time of observation
    :only_residual_decorr: if set to true, frequency decorrelation is ignored
    :return: an integer x, where the maximum amount of averaging that can be done is 2^x
    """

    curr_freq_decorr = 1 if only_residual_decorr else freq_decorr(channel_width, baseline_coord, source_radec, obs_datetime)
    curr_time_decorr, change_rate = time_decorr(baseline_coord, visint_time, source_radec, freq, phase_center, obs_datetime)

    if change_rate > 0:
        max_time_decorr = numpy.minimum(1.0, max_decorr / curr_freq_decorr)
        max_visint_time = numpy.maximum(visint_time, numpy.sqrt((6 - 6 * max_time_decorr) / (numpy.pi**2 * change_rate)))

        return numpy.floor(numpy.log2(max_visint_time / visint_time))
    else:
        return 9999

def compute_averaged_vis_count(baseline_coords, max_decorr, source_radec, phase_center, visint_time, freqs, channel_width, obs_datetime, num_time_samples, only_residual_decorr):
    """
    compute_averaged_vis_count calculates the number of visibilities post averaging without
    performing averaging. This is needed to know how much memory to allocate

    :baseline_coords: baseline coordinates in ecef
    :max_decorr: max decorrelation
    :source_radec: source coordinates
    :phase_center: pointing phase center
    :visint_time: pre-averaging (but post fs) visibility integration time
    :freqs: em frequencies
    :channel_width: channel bandwidth, assumed to be the same for all channels
    :obs_datetime: time of observation
    :num_time_samples: total number of samples on the time axis
    :only_residual_decorr: if set to true, frequency decorrelation is ignored
    :return: total number of post averaging visibilities, the number of time samples per baseline, the baselines corresponding
    to each number of time samples
    """
    total_vis = 0
    samples_per_baseline = numpy.zeros((baseline_coords.shape[0]))
    baselines_per_sample = {}

    for i, baseline_coord in enumerate(baseline_coords):
        min_levels = numpy.ceil(numpy.log2(num_time_samples))

        #for ease of data organization, we average by the smallest amount across all frequencies of the baseline
        for freq in freqs:
            levels = averaging_levels(max_decorr, baseline_coord, visint_time, freq, channel_width, source_radec, phase_center, obs_datetime, only_residual_decorr)
            min_levels = numpy.minimum(levels, min_levels)

        num_samples = int(numpy.ceil(num_time_samples / (2 ** min_levels)))
        
        samples_per_baseline[i] = num_samples
        if num_samples in baselines_per_sample:
            baselines_per_sample[num_samples].append(i)
        else:
            baselines_per_sample[num_samples] = [i]

        total_vis += num_samples

    return total_vis, samples_per_baseline, baselines_per_sample

#performs bda for all baselines, assumes visibilities are taken in succession with no gap (ie. from a continuous observation)
def bda(vis, max_decorr, half_fov, only_residual_decorr=False):
    """
    bda performs baseline dependent averaging along the time axis

    :vis: RASCIL visibilities to average
    :max_decorr: max decorrelation
    :half_fov: half of the field of view in radians of observed field
    :only_residual_decorr: if set to true, frequency decorrelation is ignored
    :return: RASCIL visibility structure containing the averaged visibilities. The structure is flattened to a 3D structure
    so that we still have a nice block of contiguous memory, with fake baselines being used, thus should not be used for
    any additional baseline-dependent analysis
    """
    baseline_coords = get_baseline_coords(vis.configuration.xyz)
    phase_center = (vis.attrs["phasecentre"].ra.value, vis.attrs["phasecentre"].dec.value)
    freqs = vis.frequency.data
    telescope_coords = (vis.configuration.location.lon.value, vis.configuration.location.lat.value)
    channel_width = vis.channel_bandwidth.data[0]
    visint_time = vis["integration_time"].data[0]

    obs_datetime, source_radec = get_maxdecorr_obstime_source(phase_center, half_fov, vis, telescope_coords)

    source_radec_rad = deg2rad2d(source_radec)
    phase_center_rad = deg2rad2d(phase_center)

    num_output_vis, _, _ = compute_averaged_vis_count(baseline_coords, max_decorr, source_radec_rad, phase_center_rad, visint_time, freqs, channel_width, obs_datetime, vis.vis.shape[0], only_residual_decorr)

    ntimes, nbaselines, nfreqs, npols = vis.vis.shape

    #initialise visibility structure here
    bv_times = numpy.zeros([1])
    bv_vis = numpy.zeros([1, num_output_vis, nfreqs, npols]).astype("complex")
    bv_flags = numpy.zeros([1, num_output_vis, nfreqs, npols]).astype("int")
    bv_weight = numpy.zeros([1, num_output_vis, nfreqs, npols])
    bv_uvw = numpy.zeros([1, num_output_vis, 3])
    bv_integration_time = numpy.array([visint_time])
    bv_baselines = pandas.MultiIndex.from_tuples(
        ingest.generate_fake_baselines(num_output_vis), names=("antenna1", "antenna2")
    )

    config = create_config(vis)

    #tranpose time and baseline because we iterate first by baseline, this makes things much faster
    #copies are created as I want to avoid modifying the input
    t_flags = numpy.copy(vis.flags)
    t_vis = numpy.copy(vis.vis)
    t_weight = numpy.copy(vis.weight)
    t_uvw = numpy.copy(vis.uvw)

    t_flags = numpy.swapaxes(t_flags, 0, 1)
    t_vis = numpy.swapaxes(t_vis, 0, 1)
    t_weight = numpy.swapaxes(t_weight, 0, 1)
    t_uvw = numpy.swapaxes(t_uvw, 0, 1)    


    curr_vis_idx = 0
    for baseline_idx in range(nbaselines):

        min_levels = numpy.ceil(numpy.log2(vis.vis.shape[0]))

        for freq in freqs:
            levels = averaging_levels(max_decorr, baseline_coords[baseline_idx], visint_time, freq, channel_width, source_radec_rad, phase_center_rad, obs_datetime, only_residual_decorr)
            min_levels = numpy.minimum(levels, min_levels)

        vis_to_average = 2 ** int(numpy.round(min_levels))

        for time_idx in range(0, ntimes, vis_to_average):
            vis_to_average_clamped = numpy.minimum(vis_to_average, ntimes - time_idx)

            for freq_idx in range(nfreqs):
                for pol_idx in range(npols):
                    non_flagged = 0

                    for i in range(vis_to_average_clamped):
                        curr_flag = t_flags[baseline_idx, time_idx+i, freq_idx, pol_idx]

                        #we only average non-flagged visibilities
                        if curr_flag == 0:
                            non_flagged += 1

                            bv_vis[0, curr_vis_idx, freq_idx, pol_idx] += t_vis[baseline_idx, time_idx+i, freq_idx, pol_idx]
                            bv_weight[0, curr_vis_idx, freq_idx, pol_idx] += t_weight[baseline_idx, time_idx+i, freq_idx, pol_idx]

                    if non_flagged > 0:
                        bv_vis[0, curr_vis_idx, freq_idx, pol_idx] /= non_flagged
                        bv_flags[0, curr_vis_idx, freq_idx, pol_idx] = 0
                    else:
                        bv_flags[0, curr_vis_idx, freq_idx, pol_idx] = 1

            for i in range(vis_to_average_clamped):
                bv_uvw[0, curr_vis_idx, :] += t_uvw[baseline_idx, time_idx+i, :] / vis_to_average_clamped

            curr_vis_idx += 1
    
    return Visibility.constructor(
                uvw=bv_uvw,
                baselines=bv_baselines,
                time=bv_times,
                frequency=freqs,
                channel_bandwidth=[channel_width]*len(freqs),
                vis=bv_vis,
                flags=bv_flags,
                weight=bv_weight,
                integration_time=bv_integration_time,
                configuration=config,
                phasecentre=vis.phasecentre,
                polarisation_frame=vis.visibility_acc.polarisation_frame,
                source=vis.source,
                meta=vis.meta,
            )


def generate_telescope_histogram(vis, show_partitions = False, num_partitions = 1, show_pixel_hist=False, n_pixels=512, pixel_size=0.001, output_file=None):
    """
    generate_telescope_histogram plots a visibility histogram by uv distsance to the center for some given set of visibilities

    :vis: RASCIL visibilities to plot
    :show_partitions: show approximate ideal partitioning
    :num_partitions: number of partitions
    :show_pixel_hist: use pixel units, otherwise lambda units will be used
    :n_pixels: number of pixels of image along one dimension, assumed to be square
    :pixel_size: pixel angular resolution in radians
    :output_file: output file for plot to be saved
    :return: approximate partition line locations
    """

    uvlambdas = numpy.array(vis.visibility_acc.uvw_lambda[..., 0, 0:2])
    uvlambdas = uvlambdas.reshape((int(uvlambdas.shape[0] * uvlambdas.shape[1]), 2))
    distances = numpy.zeros(uvlambdas.shape[0])
    for i, uvlambda in enumerate(uvlambdas):
        distances[i] = numpy.sqrt(uvlambda[0]**2 + uvlambda[1]**2)

    distances = numpy.sort(distances)

    if show_pixel_hist:
        distances = distances * pixel_size * n_pixels

    hist_stats = plt.hist(distances, bins=n_pixels)

    counts = hist_stats[0]
    counts_sum = numpy.sum(counts)
    buckets = hist_stats[1]

    ymin = 0
    ymax = numpy.max(counts)

    divider_step = counts_sum / num_partitions - 1
    acc = 0
    vlines = [0]

    for i, count in enumerate(counts):
        acc += count
        if acc > divider_step:
            acc -= divider_step
            vlines.append(buckets[i])

    if show_partitions:
        plt.vlines(vlines, ymin, ymax, colors='r')

    plt.xlabel('Dist to center (pix)')
    plt.ylabel('Num vis')

    if output_file is not None:
        plt.savefig(output_file, pad_inches=0.0, bbox_inches='tight')

    return vlines


def get_polarization(vis):
    """
    get_polarization returns a list of polarizations for some RASCIL visibility structure

    :vis: RASCIL visibilities to plot
    :return: list of polarization channels
    """

    if vis.visibility_acc.polarisation_frame.type == "linear":
        polarization = ["XX", "XY", "YX", "YY"]
    elif vis.visibility_acc.polarisation_frame.type == "linearnp":
        polarization = ["XX", "YY"]
    elif vis.visibility_acc.polarisation_frame.type == "stokesI":
        polarization = ["I"]
    elif vis.visibility_acc.polarisation_frame.type == "circular":
        polarization = ["RR", "RL", "LR", "LL"]
    elif vis.visibility_acc.polarisation_frame.type == "circularnp":
        polarization = ["RR", "LL"]
    elif vis.visibility_acc.polarisation_frame.type == "stokesIQUV":
        polarization = ["I", "Q", "U", "V"]
    elif vis.visibility_acc.polarisation_frame.type == "stokesIQ":
        polarization = ["I", "Q"]
    elif vis.visibility_acc.polarisation_frame.type == "stokesIV":
        polarization = ["I", "V"]
    else:
        raise ValueError(
            "Unknown visibility polarisation %s"
            % (vis.visibility_acc.polarisation_frame.type)
        )

    return polarization


def get_antennas(vis):
    """
    get_antennas gets the antennas from a RASCIL visibility structure

    :vis: RASCIL visibilities to plot
    :return: list of antennas, being their index and coordinates
    """
    antennas = []
    names = vis.configuration.names.data
    xyz = vis.configuration.xyz.data
    for i in range(len(names)):
        antennas.append(
            Antenna(i, Stand(names[i], xyz[i, 0], xyz[i, 1], xyz[i, 2]))
        )

    return antennas


def bda_and_export_to_ms2(msname, vis, max_decorr, half_fov, only_residual_decorr=False):
    """
    bda_and_export_to_ms2 performs baseline dependent averaging along the time axis and saves the output to a measurement set.
    This function uses a modified version of RASCIL's writing to ms function as the original does not support the non 4d block
    structure of the averaged visibilities

    :msname: name of output measurement set
    :vis: RASCIL visibilities to average
    :max_decorr: max decorrelation
    :half_fov: half of the field of view in radians of observed field
    :only_residual_decorr: if set to true, frequency decorrelation is ignored
    """

    baseline_coords = get_baseline_coords(vis.configuration.xyz)
    phase_center = (vis.attrs["phasecentre"].ra.value, vis.attrs["phasecentre"].dec.value)
    freqs = vis.frequency.data
    telescope_coords = (vis.configuration.location.lon.value, vis.configuration.location.lat.value)
    channel_width = vis.channel_bandwidth.data[0]
    visint_time = vis["integration_time"].data[0]

    obs_datetime, source_radec = get_maxdecorr_obstime_source(phase_center, half_fov, vis, telescope_coords)

    source_radec_rad = deg2rad2d(source_radec)
    phase_center_rad = deg2rad2d(phase_center)

    _, samples_per_baseline, baselines_per_samples = compute_averaged_vis_count(baseline_coords, max_decorr, source_radec_rad, phase_center_rad, visint_time, freqs, channel_width, obs_datetime, vis.vis.shape[0], only_residual_decorr)

    #tranpose time and baseline because we iterate first by baseline, this makes things much faster
    #copies are created as I want to avoid modifying the input
    t_flags = numpy.copy(vis.flags)
    t_vis = numpy.copy(vis.vis)
    t_weight = numpy.copy(vis.weight)
    t_uvw = numpy.copy(vis.uvw)
    t_time = numpy.copy(vis.time)

    t_flags = numpy.swapaxes(t_flags, 0, 1)
    t_vis = numpy.swapaxes(t_vis, 0, 1)
    t_weight = numpy.swapaxes(t_weight, 0, 1)
    t_uvw = numpy.swapaxes(t_uvw, 0, 1)    

    ntimes, nbaselines, nfreqs, npols = vis.vis.shape
    polarization = get_polarization(vis)
    source_name = vis.source
    obs_start = Time(vis.datetime[0].data).datetime
    antennas = get_antennas(vis)


    tbl = msv2.Ms(
        msname,
        ref_time=obs_start,
        source_name=source_name,
        frame=vis.configuration.attrs["frame"],
        if_delete=True,
    )

    tbl.set_stokes(polarization)
    tbl.set_frequency(vis["frequency"].data, vis["channel_bandwidth"].data)
    tbl.set_geometry(vis.configuration, antennas)

    k = 0

    nants = len(vis.configuration.names.data)

    for samples in baselines_per_samples:
        baseline_indices = baselines_per_samples[samples]
        baselines_per_level = len(baseline_indices)

        vis_to_average = int(numpy.ceil(ntimes / samples))
        curr_level_visint_time = visint_time * vis_to_average

        #maybe start with transposed instead and then flip the arrays, test later
        bv_times = []
        bv_vis = numpy.zeros([samples, baselines_per_level, nfreqs, npols]).astype("complex")

        bv_flags = numpy.zeros([samples, baselines_per_level, nfreqs, npols]).astype("int")
        bv_weight = numpy.zeros([samples, baselines_per_level, nfreqs, npols])
        bv_uvw = numpy.zeros([samples, baselines_per_level, 3])

        baseline_pairs = []
        for baseline_idx in baseline_indices:
            pair = vis.baselines[baseline_idx].data.tolist()
            baseline_pairs.append((antennas[pair[0]], antennas[pair[1]]))

        for output_baseline_idx, baseline_idx in enumerate(baseline_indices):
            output_time_idx = 0
            for time_idx in range(0, ntimes, vis_to_average):
                vis_to_average_clamped = numpy.minimum(vis_to_average, ntimes - time_idx)
                #output_time_idx = time_idx // vis_to_average

                for freq_idx in range(nfreqs):
                    for pol_idx in range(npols):
                        non_flagged = 0

                        for i in range(vis_to_average_clamped):
                            curr_flag = t_flags[output_baseline_idx, time_idx+i, freq_idx, pol_idx]

                            #we only average non-flagged visibilities
                            if curr_flag == 0:
                                non_flagged += 1

                                bv_vis[output_time_idx, output_baseline_idx, freq_idx, pol_idx] += t_vis[baseline_idx, time_idx+i, freq_idx, pol_idx]
                                bv_weight[output_time_idx, output_baseline_idx, freq_idx, pol_idx] += t_weight[baseline_idx, time_idx+i, freq_idx, pol_idx]

                        if non_flagged > 0:
                            bv_vis[output_time_idx, output_baseline_idx, freq_idx, pol_idx] /= non_flagged
                            bv_flags[output_time_idx, output_baseline_idx, freq_idx, pol_idx] = 0
                        else:
                            bv_flags[output_time_idx, output_baseline_idx, freq_idx, pol_idx] = 1

                for i in range(vis_to_average_clamped):
                    bv_uvw[output_time_idx, output_baseline_idx, :] += t_uvw[baseline_idx, time_idx+i, :] / vis_to_average_clamped

                bv_times.append(t_time[time_idx])

                output_time_idx += 1

        for itime in range(samples):
            for ipol, pol in enumerate(polarization):
                tbl.add_data_set(
                    bv_times[itime],
                    curr_level_visint_time,
                    baseline_pairs,
                    bv_vis[itime, ..., ipol],
                    weights=bv_weight[itime, ..., ipol],
                    flags=bv_flags[itime, ..., ipol],
                    pol=pol,
                    source=source_name,
                    phasecentre=vis.phasecentre,
                    uvw=bv_uvw[itime, :, :],
                )
        


    tbl.write()
    tbl.close()