#!/usr/bin/env python

"""\
notebook specific code used for analysis of parallel method
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import numpy
import matplotlib

from scipy import signal
import matplotlib.pyplot as plt

from radioimaging.util import util
from radioimaging.evaluation import evaluation
import split

def plot_snr_across_bands(title, path, cuts, bands, datasets, colours, linestyles, labels):
    for i, dataset in enumerate(datasets):
        gt = util.fromfits("../data/" + dataset + "_full_gt.fits")
        for j, cut in enumerate(cuts):
            reconstructed = util.fromfits(path + dataset + "_" + str(cut) + ".fits")
            snrs = []
            for band in bands:
                curr_band_filter = split.bandpass(reconstructed, band[0], band[1])
                
                curr_recon = signal.fftconvolve(reconstructed, curr_band_filter, mode='same')
                curr_gt = signal.fftconvolve(gt, curr_band_filter, mode='same')
                
                curr_snr = evaluation.compute_snr(curr_gt, curr_recon)
                snrs.append(curr_snr)
                
            plt.plot([x[0] for x in bands], snrs, c=colours[i], linestyle=linestyles[j], \
                     label=labels[i] + " " + str(cut), lw=3)

    plt.legend(loc="upper right", fontsize=10)

    plt.title(title)

    plt.xlabel("Frequency band", fontsize=12)
    plt.ylabel("SNR", fontsize=12)
    plt.show()

def plot_snr_across_bands_allvis(title, path, bands, datasets, colours, linestyles, labels):
    for i, dataset in enumerate(datasets):
        gt = util.fromfits("../data/" + dataset + "_full_gt.fits")

        reconstructed = util.fromfits(path + dataset + "_55.fits")
        snrs = []
        for band in bands:
            curr_band_filter = split.bandpass(reconstructed, band[0], band[1])
            
            curr_recon = signal.fftconvolve(reconstructed, curr_band_filter, mode='same')
            curr_gt = signal.fftconvolve(gt, curr_band_filter, mode='same')
            
            curr_snr = evaluation.compute_snr(curr_gt, curr_recon)
            snrs.append(curr_snr)
            
        plt.plot([x[0] for x in bands], snrs, c=colours[i], linestyle=linestyles[0], \
                 label=labels[i], lw=3)

    plt.legend(loc="upper right", fontsize=10)

    plt.title(title)

    plt.xlabel("Frequency band", fontsize=12)
    plt.ylabel("SNR", fontsize=12)
    plt.show()

def plot_snr_across_bands_and_cases(dataset, cut, dataset_title, bands, cases, case_paths, colours):
    gt = util.fromfits("../data/" + dataset + "_full_gt.fits")

    for i, case in enumerate(cases):
        filename = case_paths[i] + dataset + "_" + str(cut) + ".fits" if case.lower() != "single step" else case_paths[i] + dataset + "_55.fits"
        reconstructed = util.fromfits(filename)
        snrs = []
        for band in bands:
            curr_band_filter = split.bandpass(reconstructed, band[0], band[1])
            
            curr_recon = signal.fftconvolve(reconstructed, curr_band_filter, mode='same')
            curr_gt = signal.fftconvolve(gt, curr_band_filter, mode='same')
            
            curr_snr = evaluation.compute_snr(curr_gt, curr_recon)
            snrs.append(curr_snr)
            
        plt.plot([x[0] for x in bands], snrs, c=colours[i], label=cases[i], lw=3)

    plt.legend(loc="upper right", fontsize=10)

    plt.title("SNRs by strategies for " + dataset_title)

    plt.xlabel("Frequency band", fontsize=12)
    plt.ylabel("SNR", fontsize=12)
    plt.show()