import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from visibilities import generate_visibilities
from astropy.coordinates import SkyCoord
from astropy import units as u

def filter_mstep(x, deltas, ells, sigma2s):
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

#assumes signal length and window length are either both odd or even
def windowed_sinc_lp(cutoff, signal_length, window_size, window_type, kaiser_beta=1):
    half_signal_length = signal_length // 2 if signal_length % 2 == 0 else signal_length // 2 + 1
    tf_nw = np.array([1 if x < cutoff else 0 for x in range(0, half_signal_length)])
    tf_nw = np.concatenate((tf_nw[::-1], tf_nw if signal_length % 2 == 0 else tf_nw[1:]))
    spatial_nw = np.fft.ifftshift(np.fft.ifft(np.fft.fftshift(tf_nw))).real

    spatial_w, window = apply_window(spatial_nw, window_size, window_type, kaiser_beta)
    tf_w = np.fft.ifftshift(np.fft.fft(np.fft.fftshift(spatial_w)))

    return tf_nw[signal_length // 2:], spatial_nw, tf_w[signal_length // 2:], spatial_w, window[signal_length//2-signal_length//2:signal_length//2+signal_length//2 + signal_length % 2]

#assumes signal_length is even
def parks_mcclellan_lp(passband_range, transition_hw, signal_length, sigma2):
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
    signal_length = lp_spatial.shape[0]
    n = np.arange(signal_length) / signal_length
    shift = np.exp(1j * 2 * np.pi * bp_center * n)
    shifted_filter_spatial = shift * np.fft.fftshift(lp_spatial)
    shifted_filter_fourier = np.fft.ifftshift(np.fft.fft(shifted_filter_spatial))

    return shifted_filter_fourier[signal_length//2:signal_length], np.fft.ifftshift(shifted_filter_spatial)[signal_length-signal_length//2:signal_length+signal_length//2]

def lp_to_hp(lp_spatial):
    signal_length = lp_spatial.shape[0]
    hp_spatial = -lp_spatial
    hp_spatial[signal_length//2] += 1
    hp_fourier = np.fft.ifftshift(np.fft.fft(np.fft.fftshift(hp_spatial)))[signal_length//2:]

    return hp_fourier, hp_spatial

#assumes same length filters, and that there is a sigma2 for each filter
#can optionally provide a cutoff to the filters so that they completely adhere to some given supports
def normalize_filters(filters_freq, sigma2s, cutoffs=None):
    signal_length = filters_freq[0].shape[0]
    normalized_filters = [x.copy() for x in filters_freq]
    filter_sum = np.zeros(signal_length, dtype=np.complex)

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

#assumes the frequencies go from zero, thus the actual signal will be 2x the length of the frequency response of the filter
def compute_spatial_filters(filters_freq):
    filters_spatial = []
    filters_freq_cutoff = []
    signal_length = filters_freq[0].shape[0]
    for i, filt in enumerate(filters_freq):
        expanded_filt = np.concatenate((filt[::-1], filt))
        filt_spatial = np.fft.ifftshift(np.fft.ifft(np.fft.fftshift(expanded_filt)))
        filters_spatial.append(filt_spatial)

    return filters_spatial

def compute_cutoffs(ells, hw, signal_length):
    cutoffs = []
    cutoffs.append((0, ells[0] + hw))

    for i, ell in enumerate(ells):
        if i == len(ells) - 1:
            cutoffs.append([ells[i] - hw, signal_length])
        else:
            cutoffs.append([ells[i] - hw, ells[i+1] + hw])

    return cutoffs

def plot_filters(filt_names, filter_vals, filter_sum, colours, plot_mask, xmin=None, xmax=None, dbscale=False, spatial=False):
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
    plt.show()

def find_support(func, threshold, offset=0):
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


def generate_telescope_histogram(dec, telescope, freq, signal_length, time_res, ha_interval = [-0.5, 0.5], \
    show_partitions = False, num_partitions = 1, show_pixel_hist=False, n_pixels=512, pixel_size=0.001, ignore_autocorrelation=False):

    phasecentre = SkyCoord(ra=+56.0 * u.deg, dec=dec * u.deg, frame='icrs', equinox='J2000')

    vis = generate_visibilities(phasecentre, ha_interval, integration_time=3600/time_res, tel=telescope, frequencies=[freq])

    uvlambdas = vis.visibility_acc.uvw_lambda[..., 0, 0:2]
    uvlambdas = uvlambdas.reshape((int(uvlambdas.shape[0] * uvlambdas.shape[1]), 2))
    distances = np.zeros(uvlambdas.shape[0])
    for i, uvlambda in enumerate(uvlambdas):
        distances[i] = np.sqrt(uvlambda[0]**2 + uvlambda[1]**2)

    distances = np.sort(distances)

    if ignore_autocorrelation:
        idx = 0

        while distances[idx] == 0:
            idx += 1

        distances = distances[idx:]

    if show_pixel_hist:
        distances = distances * pixel_size * n_pixels

    hist_stats = plt.hist(distances, bins=signal_length)

    counts = hist_stats[0]
    counts_sum = np.sum(counts)
    buckets = hist_stats[1]

    if show_partitions:
        ymin = 0
        ymax = np.max(counts)

        divider_step = counts_sum / num_partitions - 1
        acc = 0
        vlines = []

        for i, count in enumerate(counts):
            acc += count
            if acc > divider_step:
                acc -= divider_step
                vlines.append(buckets[i])

        plt.vlines(vlines, ymin, ymax, colors='r')

        print(vlines)

    plt.xlabel('Dist to center (pix)')
    plt.ylabel('Num vis')

#assumes square filters, and that the 1d frequency response starts from the center and stops at signal length / 2 (i.e. not mirrored along the center)
def freq1d_to_radial2d(freq1d, signal_length):
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