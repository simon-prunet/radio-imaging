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

from images import create_image_from_fits

import logging

log = logging.getLogger()
log.setLevel(logging.DEBUG)
log.addHandler(logging.StreamHandler(sys.stdout))

def generate_visibilities(phasecentre, ha_interval, integration_time=120., tel='MEERKAT',
                          rmax=None,frequencies=None,channel_bandwidths=None, elevation_limit=None):

	'''
	This routine will use rascil tools to create a visibility structure, for observations taken from a given telescope, in the
	direction of a (SkyCoord) phase_center.
	By default times are expressed as HA in radians. 
	frequencies: array of values in Hz. If frequencies is not specified, a single value of 1e9 Hz is assumed.
	channel_bandwidths: array of values in Hz. If not specified, a single value of 1e6 Hz is assumed.
	'''
	
	if not isinstance(phasecentre, SkyCoord):
		print("phasecentre should be a SkyCoord instance")
		return
	config = create_named_configuration(tel,rmax=rmax)
	if frequencies is None:
		frequencies = np.array([1.e9])
	if channel_bandwidths is None:
		channel_bandwidths = np.array([1.e6])

	# Now compute number of integration times and corresponding HAs
	dtime_hr = integration_time / 3600.
	ntimes = int((ha_interval[1]-ha_interval[0])/dtime_hr)
	# Centered w.r.t. transit, in radian
	times = np.linspace(ha_interval[0]+dtime_hr/2., ha_interval[1]-dtime_hr/2.,ntimes) *np.pi / 12.0 
	vt = create_visibility(config, times, frequencies, channel_bandwidth=channel_bandwidths, 
	                       weight=1.0, phasecentre=phasecentre, polarisation_frame=PolarisationFrame('stokesI'),
	                       elevation_limit=elevation_limit)

	return(vt)

def select_visibilities(vt, uvmin=0, uvmax=np.inf):

	'''
	Make a new copy of the visibility structure,
	flag baselines according to uv radius
	'''
	nvt = vt.copy(deep=True)
	nvt.visibility_acc.select_uv_range(uvmin=uvmin,uvmax=uvmax)
	return(nvt)


def visibilities_from_image(vt,fitsfile,scale_factor=1.0,return_cellsize=True, return_image=False):
	'''
	Load an image from fits file, and create an Image structure, with 
	corresponding wcs info
	'''
	advice = advise_wide_field(vt, guard_band_image=3.0, delA=0.1, facets=1, 
		oversampling_synthesised_beam=4.0)
	cellsize = scale_factor*advice['cellsize']

	im = create_image_from_fits(fitsfile,frequency=vt.frequency.data,cellsize=cellsize,phasecentre=vt.phasecentre)
	ivt = predict_ng(vt,im,context='2d')
	if return_cellsize and return_image:
		return(ivt,cellsize,im)
	elif (return_cellsize and not return_image):
		return(ivt,cellsize)
	elif (return_image and not return_cellsize):
		return(ivt,im)
	else:
		return(ivt)

def dirty_psf_from_visibilities(vt,cellsize,npix=512):

	'''
	Now that visibility data corresponds to Nifty-Gridder sampling of Fourier plane of data
	(using visibilities_from_images), visibilities are used to create dirty image and corresponding psf
	:param vt: visibility data (not just baseline info)
	:param cellsize: pixel size in radian
	:param npix: number of pixels along each axis of the output images (dirty, psf)
	'''

	# First create empty rascil Image instance from visibilities
	model = create_image_from_visibility(vt,cellsize=cellsize,npixel=npix)
	print ("Model image plate scale (arcsec) is %e"%np.abs((model.image_acc.wcs.wcs.cdelt[0]*3600)))
	dirty, sumwt = invert_ng(vt, model, context='2d')
	psf, sumwt   = invert_ng(vt, model, context='2d', dopsf=True)

	return (dirty,psf)










