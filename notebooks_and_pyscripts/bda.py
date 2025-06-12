#The BDA implementation here was largely taken from https://github.com/HERA-Team/baseline_dependent_averaging and adapted to be used with RASCIL. 
#This itself is based on the Wijnholds 2018 paper "Baseline-dependent averaging in radio interferometry"
#It is probably a good idea to in the future adopt this to be performed while doing the ms ingestion

import astropy
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, Angle
from astropy import constants as const
import astropy.units as u
from astropy.time import Time

import matplotlib.pyplot as plt
import matplotlib
import numpy

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

import numpy

#earth rot speed, in rads per second
omega = 7.2925e-5

#PS. all ha/ra/dec values are assumed to be in radians (except for the initial call to bda, which is in degrees), frequencies in hertz (1/s)

#pseudo config for the visibility constructor, doesn't really matter what is put here since we never plan to save the visibilities
def create_config(vis):
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
    nants = antenna_coords.shape[0]
    nbaselines = nants * (nants + 1) / 2
    baseline_coords = numpy.zeros((nbaselines, 3))

    for a1 in range(nants):
        for a2 in range(nants):
            a1coord = antenna_coords[a1,:]
            a2coord = antenna_coords[a2,:]

            baseline_coords[a1 * nants + a2,:] = a1coord - a2coord

    return baseline_coords

#assumes 2d coords
def deg2rad2d(coords):
    return ((coords[0] * u.deg).to(u.rad).value, (coords[1] * u.deg).to(u.rad).value)

def generate_fake_baselines(nbaselines):
    for i in range(0, nbaselines):
        yield 0, 0

#dudt taken from wijnholds 2018 paper equation 42
def dudt(lx, ly, ha, omega, lambd):
    return (lx * numpy.cos(ha) - ly * numpy.sin(ha)) * omega / lambd

def dvdt(lx, ly, ha, dec, lambd):
    return (lx * numpy.sin(dec) * numpy.sin(ha) + ly * numpy.sin(dec) * numpy.cos(ha)) * omega / lambd

def radec_to_ecefunit(ra, dec, t):
    gst = t.sidereal_time('apparent', 'greenwich')
    gst_rad = (gst.hour * 15 * u.deg).to(u.rad).value

    return (cos(dec) * cos(ra + gst_rad), cos(dec) * sin(ra + gst_rad), sin(dec))

#max decorrelation is when ha=0, so we need to compute some observation time where this is the case
def get_maxdecorr_obstime(source_radec, telescope_coords):
    target_ra = (source_radec[0] * u.deg).to(u.hourangle).value
    observer_location = EarthLocation(lon=telescope_coords[0], lat=observer_lat[1])

    #arbitrary time can be used, this doesn't matter so much
    initial_time = Time('2020-06-12 12:00:00', scale='utc', location=observer_location)

    for i in range(5):
        lst = initial_time.sidereal_time('mean')
        ha = target_ra - lst
        time_correction = (ha.to(u.hour).value) / 1.00273791
        initial_time += time_correction * u.hour

    return initial_time

#freq decorrelation
def freq_decorr(channel_width, baseline_coords_ecef, source_radec, obs_datetime):
    src_coords = radec_to_ecefunit(source_radec[0], source_radec[1], obs_datetime)
    tau_g = (src_coords[0] * baseline_coords_ecef[0] + src_coords[1] * baseline_coords_ecef[1] + src_coords * baseline_coords_ecef[2]) / const.c

    #we don't multiply by pi since numpy's sinc does this already, 1 is no decorrelation
    return numpy.sinc(channel_width * tau_g)

#time decorrelation from wijnholds 2018 paper eq 41, with the approximation being used instead of the full equation which includes the w-term
#phase center is a 2d tuple given as ra, dec
#computes the maximum decorrelation, which is in a corner away from the center of the image
def time_decorr(baseline_coord, total_int_time, source_radec, freq, phase_center, obs_datetime):
    lambd = const.c / freq

    lst = obs_datetime.sidereal_time("mean")
    lst_rad = (lst.hour * 15 * u.deg).to(u.rad).value
    ha = lst_rad - source_radec[0]
    dec = source_radec[1]

    du = dudt(baseline_coord[0], baseline_coord[1], ha, omega, lambd)
    dv = dvdt(baseline_coord[0], baseline_coord[1], ha, dec, omega, lambd)

    #small angles approximation, this may need to change if we are simultaneously observing large fields, needed because ra converges towards celestial poles
    l = (source_radec[0] - phase_center[0]) * cos(phase_center[1])
    m = source_radec[1] - phase_center[1]

    rfac = (l * du + m * dv) ** 2

    # eq 41 of wijnholds paper
    decorr_frac = 1 - (numpy.pi**2 * total_int_time ** 2) / 6.0 * rfac
    
    #1 is no decorrelation
    return decorr_frac, rfac

#computes the number of subdivisions we can perform while maintaining a max decorrelation below what is specified
#assumes that the channel bandwidths and visibility integration times are all the same
def averaging_levels(max_decorr, baseline_coord, visint_time, freq, channel_width, source_radec, phase_center, obs_datetime):
    max_freq_decorr = freq_decorr(channel_width, baseline_coord, source_radec, obs_datetime)
    max_time_decorr, max_rfac = time_decorr(baseline_coord, visint_time, source_radec, freq, phase_center, obs_datetime)
    max_total_decorr = 1 - pre_fs_decorr * curr_post_fs_decorr

    max_visint_time = numpy.sqrt(6 * post_fs_decorr / (numpy.pi**2 * max_rfac))

    return numpy.floor(numpy.log2(max_visint_time / visint_time))

#computes the number of output visibilities. This needs to be done first so that we know how much memory to allocate for the output structure
def compute_averaged_vis_count(baseline_coords, max_decorr, source_radec, phase_center, visint_time, freqs, channel_width, obs_datetime, num_time_samples):
    total_vis = 0

    for baseline_coord in baseline_coords:
        min_levels = 9999

        #for ease of data organization, we average by the smallest amount across all frequencies of the baseline
        for freq in freqs:
            levels = averaging_levels(max_decorr, baseline_coord, visint_time, freq, channel_width, source_radec, phase_center, obs_datetime)
            min_levels = numpy.min(levels, min_levels)

        num_samples = int(numpy.ceil(num_time_samples / (2 ** min_levels)))
        total_vis += num_samples

    return total_vis

#performs bda for all baselines, assumes visibilities are taken in succession with no gap (ie. from a continuous observation)
#-only averaging over time for now as this is what was provided in wijnholds 2018
#-this function also assumes that visibilities are taken at even integration times one after another with no gap in between, if this is not the case, averaging should be performed seperately for each
#observation
#-the input vis is the expected rascil format, with a 4d array ordering of [time, baseline, frequency, pol]
#-half_fov is given in degrees, and should be the half fov of the imaged field. It is used to calculate the max decorrelation in 4 cardinal directions
#-fs_time is fringe stopping time and ignored for now, visint_time is the integration time per visibility, assumed to be the same for each visibility
#-channel width is in 
#-baseline_lx is the difference between antenna positions in ecef coordinates, [baseline, ]
#-freqs is an array of frequencies with the same number of elements as the frequency axis
#-the output vis is a flattened structure of [1, 1, all_vis_flattened, pol], this is because time, baseline, frequency no longer form a 
#block structure (to be fair on real data this is typically not the case anyways). One could also simplify this to a [1, baseline_and_time, frequency, pol] if the lowest frequency is used for each
#baseline when computing the decorrelation
def bda(vis, max_decorr, phase_center, half_fov, visint_time, freqs, channel_width, telescope_coords, pol):
    baseline_coords = get_baseline_coords(vis.configuration.xyz)

    source_radec = (phase_center[0] + half_fov, phase_center[1] + half_fov)
    obs_datetime = get_maxdecorr_obstime(source_radec, telescope_coords)

    source_radec_rad = deg2rad2d(source_radec)
    phase_center_rad = deg2rad2d(phase_center)

    num_output_vis = compute_averaged_vis_count(baseline_coords, max_decorr, half_fov, source_radec_rad, visint_time, freqs, channel_width, obs_datetime)

    ntimes, nbaselines, nfreq, npol = vis.vis.shape

    #initialise visibility structure here
    bv_times = numpy.zeros([1])
    bv_vis = numpy.zeros([1, num_output_vis, nfreq, npol]).astype("complex")
    bv_flags = numpy.zeros([1, num_output_vis, nfreq, npol]).astype("int")
    bv_weight = numpy.zeros([1, num_output_vis, nfreq, npol])
    bv_uvw = numpy.zeros([1, num_output_vis, 3])
    bv_integration_time = numpy.array([visint_time])
    bv_baselines = generate_fake_baselines(num_output_vis)
    
    for baseline in range(nbaselines):
        min_levels = 9999

        for freq in range(nfreqs):
            levels = averaging_levels(max_decorr, baseline_coord, visint_time, freq, channel_width, source_radec_rad, source_radec_rad, obs_datetime)
            min_levels = numpy.min(levels, min_levels)

        for time in range(ntimes):

    curr_vis_idx = 0
    for baseline in range(nbaselines):
        min_levels = 9999

        for freq in freqs:
            levels = averaging_levels(max_decorr, baseline_coords[baseline], visint_time, freqs[freq], channel_width, source_radec_rad, source_radec_rad, obs_datetime)
            min_levels = numpy.min(levels, min_levels)

        vis_to_average = 2 ** min_levels

        for time in range(0, ntimes, vis_to_average):
            vis_to_average_clamped = numpy.min(vis_to_average, ntimes - time)

            for freq in range(nfreqs):
                for pol in range(npol):
                    non_flagged = 0

                    for i in range(vis_to_average_clamped):
                        curr_flag = vis.flags[time+i, baseline, freq, pol]

                        #we only average non-flagged visibilities
                        if curr_flag > 0:
                            non_flagged += 1

                            bv_vis[1, curr_vis_idx, freq, pol] += vis.vis[time+i, baseline, freq, pol]
                            bv_weight[1, curr_vis_idx, freq, pol] += vis.weight[time+i, baseline, freq, pol]

                    if non_flagged > 0:
                        bv_vis[1, curr_vis_idx, freq, pol] /= non_flagged
                        bv_flags[1, curr_vis_idx, freq, pol] = 1

            for i in range(vis_to_average_clamped):
                bv_uvw[0, curr_vis_idx, :] += vis.vis[time+i, baseline, :] / vis_to_average_clamped


    config = create_config(vis)

    return Visibility.constructor(
                uvw=bv_uvw,
                baselines=bv_baselines,
                time=bv_times,
                frequency=freqs,
                channel_bandwidth=channel_width,
                vis=bv_vis,
                flags=bv_flags,
                weight=bv_weight,
                integration_time=bv_integration_time,
                configuration=config,
                phasecentre=phase_center,
                polarisation_frame=pol,
                source=vis.vis.source,
                meta=vis.vis.meta,
            )