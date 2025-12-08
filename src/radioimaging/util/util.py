#!/usr/bin/env python

"""\
Utility code used by various different modules and notebooks
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import matplotlib.pyplot as plt
import astropy
from astropy.io import fits
import numpy
import csv

from mpl_toolkits.axes_grid1 import make_axes_locatable

def write_to_csv(data, filename):
    """
    write_to_csv writes a list to a csv line. This function needs to be called multiple times should one wish to write many lines

    :param data: list of data to write
    :param filename: filename to append data to
    """
    with open(filename, 'a+', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(data)

def read_csv(filename, separate_rows=False):
    """
    read_csv reads a csv into a list of lists

    :param filename: csv file to read
    :param separate_rows: if set to false, returns a flattened list
    """
    data = []
    with open(filename, newline='') as file:
        reader = csv.reader(file, delimiter=',')
        for row in reader:
            if separate_rows:
                data.append(row)
            else:
                data += row

    return data

def tofits(data, filename):
    """
    tofits writes a 2d numpy array to a .fits file

    :param data: 2d numpy array containing data to write
    :param filename: filename of output fits file
    """
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(filename, overwrite=True)
    hdulist.close()

def fromfits(filename):
    """
    fromfits reads a .fits file to a 2d numpy array. It throws away any additional data if the stored data has more than 2 dimensions

    :param filename: filename of .fits file
    return: 2d numpy array
    """
    dat = fits.open(filename)[0].data
    while len(dat.shape) > 2:
        dat = dat[0]

    return dat

def convolve2d(signal, kernel, linear=False):
	"""
    convolve2d performs a 2d convolution between two given images. It assumes that the images have the same dimensionality for now

    :param img1: first image
    :param img2: second image
    :linear: perform linear convolution, if set to false, circular convolution is performed
    return: convolved image
    """
	if not linear:
		return numpy.fft.ifft2(numpy.fft.fft2(signal) * numpy.fft.fft2(numpy.fft.ifftshift(kernel))).real
	else:
		pad_length = (signal.shape[0] // 2, signal.shape[1] // 2)
		sig_padded = numpy.pad(signal, pad_length)
		kernel_padded = numpy.pad(kernel, pad_length)

		return (numpy.fft.ifft2(numpy.fft.fft2(sig_padded) * numpy.fft.fft2(numpy.fft.ifftshift(kernel_padded))))[pad_length[0]:signal.shape[0]+pad_length[0], pad_length[1]:signal.shape[1]+pad_length[1]].real

def exp_growth(x, low, high, steepness = 2):
    """
    exp_growth gives an exponential scaling of some linear input

    :x: linear input
    :low: lower bound of exponential
    :high: upper bound of exponential
    :steepness: steepness of exponential
    return: array of new exponential values
    """
    return low * (high / low) ** (x ** steepness)


def plot1D(x, y, xlabel, ylabel):
    """
    plot1D plots a 1d set of data

    :param x: x axis
    :param y: y axis
    :param xlabel: x axis label
    :param ylabel: y axis label
    """
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

def plotNImages(images, names, cmap, show_cbar=True, same_scale=False, scale_mul=1.0, output_file=None, additional_scale_imgs=None, hide_ticks=False, colorbar_location="bottom", cbar_labelsize=None, logNorm=False, vpadding=0):
    """
    plotNImages plots N images

    :param images: images to plot
    :param names: titles to each plot
    :param cmap: colour map
    :param same_scale: images share the same scale
    :param scale_mul: multiplier for image scale
    :param output_file: output filename to save plot to
    :param additional_scale_imgs: additional images that are not plotted but you want to take into account their scale for same_scale
    :param hide_ticks: hide scale ticks
    :param colorbar_location: colourbar location
    :param cbar_labelsize: colourbar label size
    :param logNorm: log scaling
    :param vpadding: vertical padding between images
    """

    num_images = len(images)

    fig, axes = plt.subplots(1, num_images)

    im = None

    vmin = 999999999999
    vmax = -999999999999

    if same_scale:
        for img in images:
            vmin = min(vmin, numpy.min(img))
            vmax = max(vmax, numpy.max(img))
        if additional_scale_imgs is not None:
            for img in additional_scale_imgs:
                vmin = min(vmin, numpy.min(img))
                vmax = max(vmax, numpy.max(img))

    vmin *= scale_mul
    vmax *= scale_mul
    vmin -= vpadding
    vmax += vpadding



    for i, img in enumerate(images):
        while(len(img.shape) > 2):
            img = img[0]

        if num_images > 1:
            axes[i].set_title(names[i])
            if same_scale:
                im = axes[i].imshow(img, norm="log" if logNorm else "linear",cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)
            else:
                im = axes[i].imshow(img, norm="log" if logNorm else "linear", cmap=cmap, origin='lower')

            if show_cbar:
                divider = make_axes_locatable(axes[i])
                cax = divider.append_axes(colorbar_location, size="5%", pad=0.25)

                cb = fig.colorbar(im, orientation='horizontal', cax=cax)
                #cb.formatter.set_powerlimits((-10, 10))
                cb.ax.locator_params(nbins=5)
                if cbar_labelsize is not None:
                    cb.ax.tick_params(labelsize=cbar_labelsize)
        else:
            axes.set_title(names[i])
            if same_scale:
                im = axes.imshow(img, norm="log" if logNorm else "linear", cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)
            else:
                im = axes.imshow(img, norm="log" if logNorm else "linear", cmap=cmap, origin='lower')

            if show_cbar:
                divider = make_axes_locatable(axes)
                cax = divider.append_axes(colorbar_location, size="5%", pad=0.25)

                cb = fig.colorbar(im, orientation='horizontal', cax=cax)
                cb.formatter.set_powerlimits((-10, 10))
                cb.ax.locator_params(nbins=5)
                if cbar_labelsize is not None:
                    cb.ax.tick_params(labelsize=cbar_labelsize)

    if hide_ticks:
        axes.set_xticks([])
        axes.set_yticks([])


    if output_file is not None:
        plt.savefig(output_file, pad_inches=0.0, bbox_inches='tight')
    else:
        plt.show()

def plotGDP(gt, dirty, psf, cmap):
    """
    plotGDP plots groundtruth, dirty, and psf

    :param gt: ground truth
    :param dirty: dirty
    :param psf: psf
    :param cmap: colour map
    """
    plotNImages([gt, dirty, psf], ["True Sky", "Dirty Image", "PSF"], cmap)

def plot1Dscatter(x, label):
    """
    plot1Dscatter plots a 1D scatter plot

    :param x: x values of data
    :param label: plot label
    """
    f = plt.figure()
    f.set_figheight(1)

    plt.xlabel(label)
    plt.tick_params(left = False, labelleft = False) 
    plt.scatter(x, [0] * len(x))

    plt.show()


def shift_img(img, x_offset, y_offset, zero_centered = True):
    """
    shift_img shifts an image through multiplication by a phasor

    :param img: image to shift
    :param x_offset: x offset
    :param y_offset: y offset
    :param zero_centered: if image has been fft shifted already
    :return: shifted image with the same size as the input
    """
    x_dim = img.shape[0]
    y_dim = img.shape[1]

    fimg = numpy.fft.fft2(numpy.fft.ifftshift(img)) if zero_centered else numpy.fft.fft2(img)

    u = numpy.fft.fftfreq(x_dim)
    v = numpy.fft.fftfreq(y_dim)
    U, V = numpy.meshgrid(u, v)
    fsimg = fimg * numpy.exp(-2j * numpy.pi * (x_offset * U + y_offset * V))

    return numpy.real(numpy.fft.fftshift(numpy.fft.ifft2(fsimg)) if zero_centered else numpy.fft.ifft2(fsimg))

def image_histogram_equalization(image, number_bins=256):
    """
    image_histogram_equalization performs histogram equalization on some input image

    :param image: input image
    :param number_bins: number of histogram bins
    :return: histogram equalized image
    """
    image_histogram, bins = numpy.histogram(image.flatten(), number_bins, density=True)
    cdf = image_histogram.cumsum() # cumulative distribution function
    cdf = (number_bins-1) * cdf / cdf[-1] # normalize
    image_equalized = numpy.interp(image.flatten(), bins[:-1], cdf)

    return image_equalized.reshape(image.shape)