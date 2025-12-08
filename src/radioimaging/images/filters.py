#!/usr/bin/env python

"""\
Code dealing with filters
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from astropy.coordinates import SkyCoord
from astropy import units as u

def filter_mstep(x, deltas, ells, sigma2s):
    """
    filter_mstep creates a set of 1d filters used in original multi-step paper, using sins for the transition area.
    Unlike the other filter functions, this one already takes into account the normalization by sigma.

    :x: range of filter in a list, e.g. [-2, -1, 0, 1, 2]
    :deltas: halfwidths of transitions areas
    :ells: middle of transition areas
    :sigma2s: square of sigmas
    :return: list of 1d mstep filters
    """

    x = np.abs(x)

    filter_idx = 0
    outer_lower = inner_lower = outer_upper = inner_upper = 0

    filter_vals = [0] * len(sigma2s)

    for i, sigma2 in enumerate(sigma2s):
        outer_lower = 0 if i == 0 else ells[i - 1] - deltas[i - 1]
        outer_upper =  10 * ells[i-1] if i == len(sigma2s) - 1 else ells[i] + deltas[i]
        inner_lower = 0 if i == 0 else ells[i - 1] + deltas[i - 1]
        inner_upper =  10 * ells[i-1] if i == len(sigma2s) - 1 else ells[i] - deltas[i]

        #in overlap region with lower frequency
        if x >= outer_lower and x < inner_lower:
            curr_ell = ells[i-1]
            curr_delta = deltas[i-1]
            lower_sigma2 = sigma2s[i-1]
            upper_sigma2 = sigma2s[i]

            low = 0.5*(1 - np.sin(2*np.pi*(x - curr_ell)/(4*curr_delta)))
            high = 0.5*(1 + np.sin(2*np.pi*(x - curr_ell)/(4*curr_delta)))
            tmp = np.sqrt(upper_sigma2*high**2 + lower_sigma2*low**2)
            filter_vals[i] = high / tmp

        #in non-overlapping region
        elif x >= inner_lower and x < inner_upper:
            filter_vals[i] = 1.0/np.sqrt(sigma2)
        #in overlap with higher frequency
        elif x >= inner_upper and x < outer_upper:
            curr_ell = ells[i]
            curr_delta = deltas[i]
            lower_sigma2 = sigma2s[i]
            upper_sigma2 = sigma2s[i+1]

            low = 0.5*(1 - np.sin(2*np.pi*(x - curr_ell)/(4*curr_delta)))
            high = 0.5*(1 + np.sin(2*np.pi*(x - curr_ell)/(4*curr_delta)))
            tmp = np.sqrt(upper_sigma2*high**2 + lower_sigma2*low**2)
            filter_vals[i] = low / tmp
        else:
            filter_vals[i] = 0

    return filter_vals


def create_filters_mstep(x_max, deltas, ells, sigma2s):
    """
    create_filters_mstep wrapper function for filters_mstep to prepare the various parameters

    :x_max: integer specifiying maximum range for filter
    :deltas: halfwidths of transitions areas
    :ells: middle of transition areas
    :sigma2s: square of sigmas
    :return: list of 1d filters, the list specifying the filter range, and the sum of the filters squared to check if they normalize
    correctly
    """
    xvals = list(np.arange(0, x_max))
    filter_vals = [[] for x in sigma2s]
    filter_sum = []

    s2snparr = np.array(sigma2s)

    for x in xvals:
        curr_vals = filter_mstep(x, deltas, ells, sigma2s)
        for i, f in enumerate(curr_vals):
            filter_vals[i].append(f)

        filter_sum.append(np.sum(s2snparr * curr_vals * curr_vals))

    return filter_vals, xvals, filter_sum


def apply_window(signal, window_size, window_type, kaiser_beta=1):
    """
    apply_window applies a windowing function to an input signal

    :signal: input signal
    :window_size: size of window
    :window_type: type of window
    :kaiser_beta: beta parameter of kaiser window
    :return: the windowed signal and the window
    """
    signal_length = signal.shape[0]

    if window_type == "blackman":
        window = np.blackman(window_size)
    elif window_type == "hanning":
        window = np.hanning(window_size)
    elif window_type == "hamming":
        window = np.hamming(window_size)
    elif window_type == "rectangular":
        window = np.ones(window_size)
    elif window_type == "kaiser":
        window = np.kaiser(window_size, kaiser_beta)
    elif window_type == "prolate_spheroidal":
        window = sig.windows.dpss(window_size, )
    else:
        print("unknown window, defaulting to rectangular")
        window = np.ones(window_size)

    if window_size < signal_length:
        truncated_window = np.pad(window, (signal_length - window_size)//2)
    elif window_size > signal_length:
        window_middle = window_size // 2 + 1 if window_size % 2 == 0 else window_size // 2
        upper = window_middle + signal_length // 2
        lower = window_middle - signal_length // 2 if signal_length % 2 == 0 else window_middle - signal_length // 2 - 1

        truncated_window = window[lower:upper]
    else:
        truncated_window = window

    spatial_w = signal * truncated_window

    return spatial_w, truncated_window


def windowed_sinc_lp(cutoff, signal_length, window_size, window_type, kaiser_beta=1):
    """
    windowed_sinc_lp obtains a lowpass prototype of the windowed sinc filter

    :cutoff: ideal cutoff for filter
    :signal_length: length of signal to apply filter
    :window_size: size of window
    :window_type: type of window to apply
    :kaiser_beta: beta parameter of kaiser window
    :return: non-windowed spatial and fourier transfer functions, windowed fourier and spatial transfer functions, windowing function
    """
    half_signal_length = signal_length // 2 if signal_length % 2 == 0 else signal_length // 2 + 1
    tf_nw = np.array([1 if x < cutoff else 0 for x in range(0, half_signal_length)])
    tf_nw = np.concatenate((tf_nw[::-1], tf_nw if signal_length % 2 == 0 else tf_nw[1:]))
    spatial_nw = np.fft.ifftshift(np.fft.ifft(np.fft.fftshift(tf_nw))).real

    spatial_w, window = apply_window(spatial_nw, window_size, window_type, kaiser_beta)
    tf_w = np.fft.ifftshift(np.fft.fft(np.fft.fftshift(spatial_w)))

    return tf_nw[signal_length // 2:], spatial_nw, tf_w[signal_length // 2:], spatial_w, window[signal_length//2-signal_length//2:signal_length//2+signal_length//2 + signal_length % 2]


def parks_mcclellan_lp(passband_range, transition_hw, signal_length, sigma2):
    """
    parks_mcclellan_lp obtains a lowpass prototype using the parks-mcclellan algorithm. Takes sigma into account

    :passband_range: list containing ranges for pass and stop bands
    :transition_hw: ideal transition halfwidth
    :signal_length: length of signal to apply filter
    :sigma2: sigma squared
    :return: filter in fourier and spatial domain
    """
    bands = []
    desired = []

    if passband_range[1] == signal_length//2:
        bands = [0, passband_range[0] - transition_hw, passband_range[0] + transition_hw, 0.5*signal_length]
        desired = [0, 1/sigma2]
    elif passband_range[0] == 0:
        bands = [0, passband_range[1] - transition_hw, passband_range[1] + transition_hw, 0.5*signal_length]
        desired = [1/sigma2, 0]
    else:
        bands = [0, passband_range[0] - transition_hw, passband_range[0] + transition_hw, passband_range[1] - transition_hw, passband_range[1] + transition_hw, 0.5*signal_length]
        desired = [0, 1/sigma2, 0]

    spatial = sig.remez(signal_length-1, bands, desired, fs=signal_length)
    spatial_padded = np.pad(spatial, (0, 1), mode='constant')
    freqs = np.fft.ifftshift(np.fft.fft(np.fft.fftshift(spatial_padded), signal_length))

    return freqs[signal_length//2:], spatial_padded


def lp_to_bp(lp_spatial, bp_center):
    """
    lp_to_bp obtains a bandpass filter from a lowpass prototype using a phase shift

    :lp_spatial: lowpass prototype
    :bp_center: bandpass center
    :return: bandpass filter in fourier and spatial
    """
    signal_length = lp_spatial.shape[0]
    n = np.arange(signal_length) / signal_length
    shift = np.exp(1j * 2 * np.pi * bp_center * n)
    shifted_filter_spatial = shift * np.fft.fftshift(lp_spatial)
    shifted_filter_fourier = np.fft.ifftshift(np.fft.fft(shifted_filter_spatial))

    return shifted_filter_fourier[signal_length//2:signal_length], np.fft.ifftshift(shifted_filter_spatial)[signal_length-signal_length//2:signal_length+signal_length//2]

def lp_to_hp(lp_spatial):
    """
    lp_to_hp obtains a highpass filter from a lowpass prototype using spectral inversion

    :lp_spatial: lowpass prototype
    :return: highpass filter in fourier and spatial
    """
    signal_length = lp_spatial.shape[0]
    hp_spatial = -lp_spatial
    hp_spatial[signal_length//2] += 1
    hp_fourier = np.fft.ifftshift(np.fft.fft(np.fft.fftshift(hp_spatial)))[signal_length//2:]

    return hp_fourier, hp_spatial


def normalize_filters(filters_freq, sigma2s, cutoffs=None):
    """
    normalize_filters normalizes a set of filters taking into account some sigma values so that their squares total to one.
    Assumes the filters are the same length

    :filters_freq: list of filters in frequency domain
    :sigma2s: list of sigma squares corresponding to each filter
    :cutoffs: list of cutoffs for filters if one wishes for an exact support
    :return: list of normalized filters and the sum of the filters squared to verify that the normalization is correct
    """
    signal_length = filters_freq[0].shape[0]
    normalized_filters = [x.copy() for x in filters_freq]
    filter_sum = np.zeros(signal_length, dtype=np.complex128)

    if cutoffs is not None:
        for i, filt in enumerate(normalized_filters):
            low = cutoffs[i][0]
            high = cutoffs[i][1]

            if low > 0:
                normalized_filters[i][0:low] = 0
            if high < signal_length - 1:
                normalized_filters[i][high:signal_length] = 0

    for i, coeff in enumerate(filters_freq[0]):
        total = 0
        for j, sigma2 in enumerate(sigma2s):
            total += sigma2 * np.abs(normalized_filters[j][i]) ** 2

        total = np.sqrt(total)
        curr_coeff_sum = 0

        for j, sigma2 in enumerate(sigma2s):
            normalized_filters[j][i] /= total
            curr_coeff_sum += sigma2 * normalized_filters[j][i] ** 2
        
        filter_sum[i] = curr_coeff_sum

    return normalized_filters, filter_sum


def compute_spatial_filters(filters_freq):
    """
    compute_spatial_filters translates a set of frequency filters to spatial.
    assumes the frequencies go from zero, thus the actual signal will be 2x the length of the frequency response of the filter

    :filters_freq: list of filters in frequency domain
    :return: list of filters in the spatial domain
    """
    filters_spatial = []
    filters_freq_cutoff = []
    signal_length = filters_freq[0].shape[0]
    for i, filt in enumerate(filters_freq):
        expanded_filt = np.concatenate((filt[::-1], filt))
        filt_spatial = np.fft.ifftshift(np.fft.ifft(np.fft.fftshift(expanded_filt)))
        filters_spatial.append(filt_spatial)

    return filters_spatial


def compute_cutoffs(ells, hw, signal_length):
    """
    compute_cutoffs obtains a list of ideal supports for filters defined by a set of transition centers and halfwidth

    :ells: transition area centers
    :hw: halfwidth, assumed to be constant
    :signal_length: length of signal to apply filters
    :return: list of support ranges
    """
    cutoffs = []
    cutoffs.append((0, ells[0] + hw))

    for i, ell in enumerate(ells):
        if i == len(ells) - 1:
            cutoffs.append([ells[i] - hw, signal_length])
        else:
            cutoffs.append([ells[i] - hw, ells[i+1] + hw])

    return cutoffs

def plot_filters(filt_names, filter_vals, filter_sum, colours, plot_mask, xmin=None, xmax=None, dbscale=False, spatial=False, output_filename=None):
    """
    plot_filters plots a set of given filters

    :filt_names: list of filter names
    :filter_vals: list of filters with their values
    :filter_sum: sum of filters, can be None
    :colours: list of filter colours
    :plot_mask: list of booleans if one doesn't want to plot all filters
    :xmin: minimum x value to plot
    :xmax: maximum x value to plot
    :dbscale: flag specifying if the filters should be plotted in db, if not, their actual response is used
    :spatial: flag specifying if the filters are spatial. If not, they are assumed to be frequency responses, and the
    absolute values are plotted instead
    """

    for i, filt in enumerate(filt_names):
        if plot_mask[i]:
            if dbscale:
                curr_vals = 20 * np.log10(np.abs(filter_vals[i]))
            else:
                curr_vals = filter_vals[i] if spatial else np.abs(filter_vals[i])

            plt.plot(np.array(curr_vals), label=filt, c=colours[i], alpha=0.6)

    if filter_sum is not None and not spatial:
        if dbscale:
            rescaled_sum = 20 * np.log10(np.abs(filter_sum))
        else:
            rescaled_sum = filter_sum

        plt.plot(np.abs(rescaled_sum), label="$\sum_n \sigma^2 f_n(x)^2 $", c="black")

    plt.legend(loc="lower right", fontsize=10)

    if xmax is None:
        xmax = filter_vals[0].shape[0]
    if xmin is None:
        xmin = 0

    plt.xlim(xmin, xmax)

    if output_filename is not None:
        plt.savefig(output_filename, pad_inches=0.0, bbox_inches='tight')

    plt.show()

def find_support(func, threshold, offset=0):
    """
    find_support finds the set of supports of a function above a given threshold

    :func: discretized function values
    :threshold: threshold of support. Set to zero to find the true supports
    :offset: value to offset support array by, used for visualization when plotting so that the supports don't overlap
    :return: binary array illustrating support, and list of ranges where the support is above threshold
    """
    support_arr = np.zeros(func.shape[0])
    support_idx = []

    prev_zero = True
    curr_support_start = -1
    for i, val in enumerate(func):
        curr_nonzero = np.abs(val) > threshold
        support_arr[i] = 1 + offset if curr_nonzero else 0
        if curr_nonzero and prev_zero:
            curr_support_start = i 
            prev_zero = False
        elif not curr_nonzero and not prev_zero:
            support_idx.append((curr_support_start, i + 1))
            prev_zero = True

    if not prev_zero:
        support_idx.append((curr_support_start, func.shape[0]))

    return support_arr, support_idx

#assumes square filters, and that the 1d frequency response starts from the center and stops at signal length / 2 (i.e. not mirrored along the center)
def freq1d_to_radial2d(freq1d, signal_length):
    """
    freq1d_to_radial2d converts a 1d filter in frequency to a 2d square radial filter

    :freq1d: frequency response of filter
    :signal_length: length of signal to apply filter
    :return: filter in frequency and spatial domains
    """
    filtfreq = np.zeros((signal_length, signal_length), dtype=np.complex128)

    center = signal_length//2 + 1
    filt1d_length = freq1d.shape[0]

    for y in range(signal_length):
        for x in range(signal_length):
            xd = np.abs(x - center)
            yd = np.abs(y - center)
            dist = np.sqrt(xd**2 + yd**2)

            #linear interpolation as the sample-rate on the 1d filter is not the same as the radial pixel centers
            t = dist - int(dist)
            low = int(dist)
            high = low + 1

            val = None

            if high >= filt1d_length:
                filtfreq[y, x] = freq1d[-1]
            else:
                filtfreq[y, x] = (1-t)*freq1d[low] + t*freq1d[high]

    filt_spatial = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(filtfreq))).real

    return filtfreq, filt_spatial