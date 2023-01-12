import os
import sys
sys.path.append(os.path.join('..', '..'))

import numpy as np
from ska_sdp_datamodels.configuration import create_named_configuration
from ska_sdp_datamodels.visibility import create_visibility

# results_dir = '/tmp/'

# from matplotlib import pylab

# pylab.rcParams['figure.figsize'] = (8.0, 8.0)
# pylab.rcParams['image.cmap'] = 'rainbow'

from astropy.coordinates import SkyCoord
from astropy import units as u

# from matplotlib import pyplot as plt
from ska_sdp_datamodels.science_data_model.polarisation_model import PolarisationFrame
from ska_sdp_func_python.image import deconvolve_cube, restore_cube 
from ska_sdp_func_python.imaging import invert_ng, predict_ng, \
     create_image_from_visibility, advise_wide_field 
from rascil.processing_components import show_image, create_test_image, \
    plot_uvcoverage, plot_visibility

import logging

log = logging.getLogger()
log.setLevel(logging.DEBUG)
log.addHandler(logging.StreamHandler(sys.stdout))

def generate_visibilities(phase_center, ha_interval, integration_time=120., tel='MEERKAT',
                          rmax=None,frequencies=None,channel_bandwidths=None):

	'''
	This routine will use rascil tools to create a visibility structure, for observations taken from a given telescope, in the
	direction of a (SkyCoord) phase_center.
	By default times are expressed as HA in radians. 
	frequencies: array of values in Hz. If frequencies is not specified, a single value of 1e9 Hz is assumed.
	channel_bandwidths: array of values in Hz. If not specified, a single value of 1e6 Hz is assumed.
	'''
	
	if not isinstance(phase_center, SkyCoord):
		print("phase_center should be a SkyCoord instance")
		return
	config = create_named_configuration(tel,rmax=rmax)
	if frequencies is None:
		frequencies = np.array([1.e9])
	if channel_bandwidths is None:
		channel_bandwidths = np.array([1.e6])

	# Now compute number of integration times and corresponding HAs
	ntimes = (ha_interval[1]-ha_interval[0])/integration_time
	times = np.linspace(ha_interval[0]+integration_time/2., ha_interval[1]-integration_time/2.,ntimes) # Centered w.r.t. transit
	vt = create_visibility(config, times, frequencies, channel_bandwidth=channel_bandwidths, 
	                       weight=1.0, phasecentre=phasecentre, polarisation_frame=PolarisationFrame('stokesI'))

	return(vt)




