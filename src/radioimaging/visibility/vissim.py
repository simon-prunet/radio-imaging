#!/usr/bin/env python

"""\
Code for simulating visibilities from various telescope configurations. This should be updated in the future
to account for various other features, such as subarrays, proper noise calculations, incorporating the primary
beams, and allowing for simulating very large amounts of visibilities and writing directly to an .ms, as it
is currently limited to the amount of memory available
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import numpy
from ska_sdp_datamodels.configuration import create_named_configuration
from ska_sdp_datamodels.visibility import create_visibility
from astropy.coordinates import SkyCoord
from ska_sdp_datamodels.science_data_model.polarisation_model import PolarisationFrame



def generate_visibilities(phasecentre, ha_interval, integration_time=120., tel='MEERKAT',
                          rmax=None,frequencies=None,channel_bandwidths=None, elevation_limit=None):

    """
    generate_visibilities creates a RASCIL visibility structure according to the provided observation parameters. This
    is all done in memory currently, so the number of visibilities is somewhat limited

    :phasecentre: phase center of pointing
    :ha_interval: hour angle interval of phase center
    :integration_time: post fringe stopping integration time per visibility sample
    :tel: telescope configuration
    :rmax: maximum radius of antennas to consider. This somewhat approximates subarrays but should be changed to actual
    subarrays
    :frequencies: list of EM frequencies
    :channel_bandwidths: list of bandwidths for each frequency
    :elevation_limit: visibilities observed when phase center is below this elevation are discarded
    :return: RASCIL visibility structure containing visibilities of observation. The actual values of the visibilities
    are empty and an image needs to be degridded to them to populate this
    """
    
    if not isinstance(phasecentre, SkyCoord):
        print("phasecentre should be a SkyCoord instance")
        return

    config = None

    if tel=="MID_AASTAR" or tel=="LOW_AASTAR":
        from ska_ost_array_config.array_config import LowSubArray, MidSubArray
        if tel == "MID_AASTAR":
            config = MidSubArray(subarray_type="AA*").array_config
        elif tel == "LOW_AASTAR":
            config = LowSubArray(subarray_type="AA*").array_config
    else:
        config = create_named_configuration(tel,rmax=rmax)

    if frequencies is None:
        frequencies = numpy.array([1.e9])
    if channel_bandwidths is None:
        channel_bandwidths = numpy.array([1.e6])

    # Now compute number of integration times and corresponding HAs
    dtime_hr = integration_time / 3600.
    ntimes = int((ha_interval[1]-ha_interval[0])/dtime_hr)
    # Centered w.r.t. transit, in radian
    times = numpy.linspace(ha_interval[0]+dtime_hr/2., ha_interval[1]-dtime_hr/2.,ntimes) *numpy.pi / 12.0 
    vt = create_visibility(config, times, frequencies, channel_bandwidth=channel_bandwidths, 
                           weight=1.0, phasecentre=phasecentre, polarisation_frame=PolarisationFrame('stokesI'),
                           elevation_limit=elevation_limit)

    return(vt)


def add_noise_to_vis(vis, perc, real_deviation=-1, imag_deviation=-1):
    """
    add_noise_to_vis adds noise to visibilities based on a percentage of standard deviation of the visibilities, or as a percentage
    of some hard set deviations

    :vis: noiseless visibilities
    :perc: percentage of noise
    :real_deviation: hardset real deviation
    :imag_deviation: hardset imaginary deviation
    :return: visibilities with noise
    """
    real_deviation = numpy.std(vis.vis.data.real.flatten()) if real_deviation < 0 else real_deviation
    imag_deviation = numpy.std(vis.vis.data.imag.flatten()) if imag_deviation < 0 else real_deviation

    real_deviation *= (perc / 100)
    imag_deviation *= (perc / 100)

    noise_real = numpy.random.normal(loc=0, scale=real_deviation, size=vis.vis.shape)
    noise_imag = numpy.random.normal(loc=0, scale=imag_deviation, size=vis.vis.shape)

    noise = numpy.vectorize(complex)(noise_real, noise_imag)
    vis_with_noise = vis.vis + noise

    nvis = vis.copy(deep=True)
    nvis["vis"].data = vis_with_noise

    return nvis
