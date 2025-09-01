#!/usr/bin/env python

"""\
notebook specific code used for split analysis
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import numpy
import matplotlib

from scipy import signal

from radioimaging.util import util
from radioimaging.evaluation import evaluation

def bandpass(image, low, high):
    if len(image.shape) < 2:
        print("Error: Image must have 2 dimensions or more")
        return None

    #ensure image is 2d
    while len(image.shape) > 2:
        image = image[0]

    center_x = image.shape[0] / 2# + 1 if image.shape[0] % 2 == 0 else image.shape[0] / 2 
    center_y = image.shape[1] / 2# + 1 if image.shape[1] % 2 == 0 else image.shape[1] / 2

    bandpass_filter_f = numpy.zeros(image.shape)

    low2 = low * low
    high2 = high * high

    for y in range(bandpass_filter_f.shape[1]):
        for x in range(bandpass_filter_f.shape[0]):
            centered_x = x - center_x
            centered_y = y - center_y

            dist2 = centered_x * centered_x + centered_y * centered_y

            if dist2 >= low2 and dist2 <= high2:
                bandpass_filter_f[x, y] = 1


    return numpy.real(numpy.fft.fftshift(numpy.fft.ifft2(numpy.fft.ifftshift(bandpass_filter_f))))

def banded_psnrs(path, gt, bands, nmaj_iter):
    data = []
    deconv_images = []

    prev_deconv = None
    for maj_iter in range(0, nmaj_iter):
        filename = path + "/deconv_iteration_" + str(maj_iter) + "_channel_0.fits"
        deconv = util.fromfits(filename)
        curr_deconv = None
        
        if prev_deconv is not None:
            curr_deconv = deconv + prev_deconv
        else:
            curr_deconv = deconv
        
        prev_deconv = curr_deconv
        deconv_images.append(curr_deconv)


    for band in bands:
        curr_data = ["low=" + str(band[0]) + ", high=" + str(band[1])]
        curr_band_filter = bandpass(gt, band[0], band[1])
        curr_gt = signal.fftconvolve(gt, curr_band_filter, mode='same')
        
        for maj_iter in range(0, nmaj_iter):
            deconv = deconv_images[maj_iter]
            filtered_deconv = signal.fftconvolve(deconv, curr_band_filter, mode='same')
            psnr = evaluation.compute_snr(curr_gt, filtered_deconv)
            curr_data.append(psnr)

        data.append(curr_data)
    
    return data

def banded_psnrs_lowresdirty(path, bands, cut, hw, nmaj_iter):
    deconv_images = []

    psf = util.fromfits(path + "/psf_channel0_iteration0.fits")
    lowres_constraint = util.fromfits(path + "/lowres_iteration_0_channel_0.fits")
    dirty = util.fromfits(path + "/dirty_iteration_0_channel_0.fits")

    prev_deconv = None
    for maj_iter in range(0, nmaj_iter):
        deconv_filename = path + "/deconv_iteration_" + str(maj_iter) + "_channel_0.fits"
        deconv = util.fromfits(deconv_filename)
        curr_deconv = None
        
        if prev_deconv is not None:
            curr_deconv = deconv + prev_deconv
        else:
            curr_deconv = deconv
        
        prev_deconv = curr_deconv
        deconv_images.append(curr_deconv)


    data_low = []
    data_high = []
        
    for band in bands:
        curr_data_low = ["low=" + str(band[0]) + ", high=" + str(band[1])]
        curr_data_high = ["low=" + str(band[0]) + ", high=" + str(band[1])]
        
        curr_band_filter = bandpass(psf, band[0], band[1])
        
        low_gt = signal.fftconvolve(lowres_constraint, curr_band_filter, mode='same')
        high_gt = signal.fftconvolve(dirty, curr_band_filter, mode='same')
        idfilter = bandpass(high_gt, 0, 500)

        #needed as fftconvolve with mode same creates a shift, which is a problem when convolving the deconvolved image with the psf
        high_gt = signal.fftconvolve(high_gt, idfilter, mode='same')
        
        for maj_iter in range(0, 5):
            deconv = deconv_images[maj_iter]
            
            filtered_deconv = signal.fftconvolve(deconv, curr_band_filter, mode='same')
            
            if band[1] <= (cut + hw):
                psnr = evaluation.compute_snr(low_gt, filtered_deconv)
                curr_data_low.append(psnr)
                
            if band[0] >= (cut - hw):
                convolved_deconv = signal.fftconvolve(filtered_deconv, psf, mode='same')
                psnr = evaluation.compute_snr(high_gt, convolved_deconv)
                curr_data_high.append(psnr)
        
        if band[1] <= (cut + hw):
            data_low.append(curr_data_low)
        
        if band[0] >= (cut - hw):
            data_high.append(curr_data_high)

    return data_low, data_high