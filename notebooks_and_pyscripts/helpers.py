import astropy
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.units import Quantity
from astropy import units as u

import matplotlib.pyplot as plt
import matplotlib
import numpy

import csv
from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy import signal
from scipy import stats

from ska_sdp_func_python.visibility import subtract_visibility
import os

from ska_sdp_func_python.visibility import subtract_visibility
from ska_sdp_func_python.imaging import invert_ng, predict_ng, create_image_from_visibility, advise_wide_field
from ska_sdp_datamodels.gridded_visibility import create_griddata_from_image
from ska_sdp_func_python.grid_data import grid_visibility_weight_to_griddata,griddata_visibility_reweight, griddata_merge_weights
from rascil.processing_components import create_visibility_from_ms, generate_baselines
from ska_sdp_func_python.visibility import convert_visibility_to_stokesI
from ska_sdp_datamodels.configuration.config_model import Configuration
from ska_sdp_datamodels.science_data_model.polarisation_model import (
    ReceptorFrame,
    PolarisationFrame,
)
from ska_sdp_datamodels.visibility.vis_model import Visibility

import pandas
import numpy
import gc

def write_to_csv(data, filename):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(data)

def write_nparr_to_fits(data, filename):
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(filename, overwrite=True)
    hdulist.close()

def t_test_image(recon_resid, target, window_size):
    assert(recon_resid.shape == target.shape)
    assert(len(recon_resid.shape) == 2)

    output_pvals = numpy.zeros(recon_resid.shape)

    for y in range(output_pvals.shape[1]):
        for x in range(output_pvals.shape[0]):
            start_x = max(x - window_size, 0)
            end_x = min(x + window_size + 1, output_pvals.shape[0])
            start_y = max(y - window_size, 0)
            end_y = min(y + window_size + 1, output_pvals.shape[1])

            recon_window = recon_resid[start_x:end_x, start_y:end_y]
            target_window = target[start_x:end_x, start_y:end_y]

            #ttest_res = stats.ttest_ind(recon_window.flatten(), target_window.flatten(), equal_var=False)
            #ttest_res = stats.anderson_ksamp([recon_window.flatten(), target_window.flatten()])
            ttest_res = stats.wasserstein_distance(recon_window.flatten(), target_window.flatten())
            #ttest_res = stats.mannwhitneyu(recon_window.flatten(), target_window.flatten())

            #output_pvals[x, y] = ttest_res.pvalue
            output_pvals[x, y] = ttest_res 

    return output_pvals, numpy.linalg.norm(output_pvals)


def compute_snr(gt, recon):
    difnorm = numpy.linalg.norm(gt-recon)
    if difnorm == 0:
        return 0
        
    return 20 * numpy.log10(numpy.linalg.norm(gt) / difnorm)

def compute_rmse(gt, recon):
    return numpy.sqrt(numpy.mean((gt-recon) ** 2))

def compute_maxabserr(gt, recon):
    return numpy.max(numpy.abs(gt - recon))

def compute_errstd(gt, recon):
    return numpy.std(gt - recon)

def compute_ssim(gt, recon):
    gt_mean = numpy.mean(gt)
    recon_mean = numpy.mean(recon)
    #gt_variance = numpy.variance(gt)
    #recon_variance = numpy.variance(recon)
    covariance_mat = numpy.cov([gt.flatten(), recon.flatten()])

    gt_variance = covariance_mat[0, 0]
    recon_variance = covariance_mat[1, 1]
    covariance = covariance_mat[0, 1]

    dynamic_range = numpy.max(gt) - numpy.min(gt) * 2

    c1 = (0.01 * dynamic_range) ** 2
    c2 = (0.03 * dynamic_range) ** 2
    c3 = c2 / 2

    l = (2 * gt_mean * recon_mean + c1) / (gt_mean ** 2 + recon_mean ** 2 + c1)
    c = (2 * numpy.sqrt(gt_variance) * numpy.sqrt(recon_variance) + c2) / (gt_variance + recon_variance + c2)
    s = (covariance + c3) / (numpy.sqrt(gt_variance) * numpy.sqrt(recon_variance) + c3)

    return (l ** 1) * (c ** 1) * (s ** 3)

    #return ((2 * gt_mean * recon_mean + c1) * (2 * covariance + c2)) / ((gt_mean ** 2 + recon_mean ** 2 + c1) * (gt_variance + recon_variance + c2))


def exp_growth(x, low, high, steepness = 2):
    return low * (high / low) ** (x ** steepness)

def readFits(filename):
    dat = fits.open(filename)[0].data
    while len(dat.shape) > 2:
        dat = dat[0]

    return dat

def plot1D(x, y, xlabel, ylabel):
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

def plotNImages(images, names, cmap, same_scale=False, scale_mul=1.0, output_file=None, additional_scale_imgs=None, hide_ticks=False, colorbar_location="bottom", cbar_labelsize=None, logNorm=False, vpadding=0):
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

            #divider = make_axes_locatable(axes)
            #cax = divider.append_axes(colorbar_location, size="5%", pad=0.25)

            #cb = fig.colorbar(im, orientation='horizontal', cax=cax)
            #cb.formatter.set_powerlimits((-10, 10))
            #cb.ax.locator_params(nbins=5)
            #if cbar_labelsize is not None:
            #    cb.ax.tick_params(labelsize=cbar_labelsize)

    if hide_ticks:
        axes.set_xticks([])
        axes.set_yticks([])


    if output_file is not None:
        plt.savefig(output_file, pad_inches=0.0, bbox_inches='tight')
    else:
        plt.show()

def plotGDP(gt, dirty, psf, cmap):
    plotNImages([gt, dirty, psf], ["True Sky", "Dirty Image", "PSF"], cmap)

def plotSNRvsSSIM(lambdas, path, snr_idx, gt, cmap, same_scale = False):
    snr_fn = path + "lambda_" + str(lambdas[snr_idx]) + ".fits"
    snrfile = readFits(snr_fn)
    err_snr = numpy.abs(gt - snrfile)
    snr_title = "Best SNR Image $\lambda = " + "{:.2f}".format(lambdas[snr_idx]) + "$"

    plotNImages([snrfile, err_snr], [snr_title, "SNR absolute error"], cmap, same_scale = same_scale)


def read_csv(filename, separate_rows=False):
    data = []
    with open(filename, newline='') as file:
        reader = csv.reader(file, delimiter=',')
        for row in reader:
            if separate_rows:
                data.append(row)
            else:
                data += row

    return data

def computeErrorForSavedResults(gt, lambdas, path, ofilename, errorMetric):
    errors = [0] * len(lambdas)
    for i, val in enumerate(lambdas):
        curr_filename = path + "lambda_" + str(val) + ".fits"
        curr_file = readFits(curr_filename)
        errors[i] = errorMetric(gt, curr_file)

    write_to_csv(errors, path + ofilename)

def circularConv(gt, psf):
    psf_fft = numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.ifftshift(psf)))
    gt_fft = numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.ifftshift(gt)))

    return numpy.real(numpy.fft.ifftshift(numpy.fft.ifft2(numpy.fft.fftshift(numpy.conj(psf_fft) * gt_fft))))

def linearConv(gt, psf):
    psf_padded = numpy.zeros((int(psf.shape[0] * 2), int(psf.shape[1] * 2)))
    psf_padded[int(psf.shape[0] / 2) : int(psf.shape[0] / 2 + psf.shape[0]), int(psf.shape[1] / 2) : int(psf.shape[1] / 2 + psf.shape[1])] = psf
    gt_padded = numpy.zeros((int(gt.shape[0] * 2), int(gt.shape[1] * 2)))
    gt_padded[int(gt.shape[0] / 2) : int(gt.shape[0] / 2 + gt.shape[0]), int(gt.shape[1] / 2) : int(gt.shape[1] / 2 + gt.shape[1])] = gt

    psf_fft = numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.ifftshift(psf)))
    gt_fft = numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.ifftshift(gt)))
    psf_fft_padded = numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.ifftshift(psf_padded)))
    gt_fft_padded = numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.ifftshift(gt_padded)))

    return numpy.real(numpy.fft.ifftshift(numpy.fft.ifft2(numpy.fft.fftshift(numpy.conj(psf_fft_padded) * gt_fft_padded)))[int(gt.shape[0] / 2) : int(gt.shape[0] / 2 + gt.shape[0]), int(gt.shape[1] / 2) : int(gt.shape[0] / 2 + gt.shape[0])])

def linearConvScipy(gt, psf):
    return signal.fftconvolve(gt, psf, mode='same')

def addNoiseToVis(vis, perc, real_deviation=-1, imag_deviation=-1):
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

def plot1Dscatter(x, label):
    f = plt.figure()
    f.set_figheight(1)

    plt.xlabel(label)
    plt.tick_params(left = False, labelleft = False) 
    plt.scatter(x, [0] * len(x))

    plt.show()

def subtractVis(recon, model):
    return subtract_visibility(recon, model)


def compute_windowed_var(image, window):
    estimated_variance = numpy.zeros(image.shape)

    #convolving initial signal with a variance estimation kernel with a size 2xwindow_hsize+1
    for y in range(image.shape[1]):
        for x in range(image.shape[0]):
            start_x = max(x - window, 0)
            end_x = min(x + window + 1, image.shape[0])
            start_y = max(y - window, 0)
            end_y = min(y + window + 1, image.shape[1])

            estimated_variance[x, y] = numpy.var(image[start_x:end_x, start_y:end_y])
            
    return estimated_variance


#bandpass filter in pixels, filters in the range of [low, high]
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

def lambdatests_allvis(lambdas, path, dirty, psf, gt, wavelet_type_idx, niter, runtests, lp_filter=None, stats_type = "psnr"):
    if lp_filter is not None:
        gt = signal.fftconvolve(gt, lp_filter, mode='same')

    if runtests:
        snrs = [0] * len(lambdas)

        tmp_dirty_fn = "tmp_dirty.fits"
        tmp_psf_fn = "tmp_psf.fits"

        write_nparr_to_fits(dirty, tmp_dirty_fn)
        write_nparr_to_fits(psf, tmp_psf_fn)

        for i, curr_lambda in enumerate(lambdas):
            curr_output_name = path + "lambda_" + str(curr_lambda) + ".fits"
            os.system("julia ../julia_rascil_scripts/make_fullres.jl " + str(curr_lambda) + " " + tmp_psf_fn + " " \
                      + tmp_dirty_fn + " " + str(wavelet_type_idx) + " " + str(niter) + " " + curr_output_name)
            recon = fits.open(curr_output_name)[0].data

            if lp_filter is not None:
                recon = signal.fftconvolve(recon, lp_filter, mode='same')

            if stats_type == "psnr":
                snrs[i] = compute_snr(gt, recon)
            elif stats_type == "ispace_resid_norm":
                recon = signal.fftconvolve(recon, psf, mode='same')
                identity_filt = bandpass(dirty, 0, dirty.shape[0] * 2)
                dirty = signal.fftconvolve(dirty, identity_filt, mode='same')

                snrs[i] = numpy.linalg.norm(recon - dirty)
            else:
                snrs[i] = compute_snr(gt, recon)
            
        write_to_csv(snrs, path + "snr.dat")

    snrs = read_csv(path + "snr.dat")
    snrs = [float(x) for x in snrs]

    return snrs


def lambdatests_highvis(lambdas, path, dirty, psf, gt, lowres, wavelet_type_idx, niter, cut_center, cut_hw, runtests):
    if runtests:
        snrs = [0] * len(lambdas)

        tmp_dirty_fn = "tmp_dirty.fits"
        tmp_psf_fn = "tmp_psf.fits"
        tmp_low_fn = "tmp_low.fits"

        write_nparr_to_fits(dirty, tmp_dirty_fn)
        write_nparr_to_fits(psf, tmp_psf_fn)
        write_nparr_to_fits(lowres, tmp_low_fn)

        vis_noise = numpy.mean(compute_windowed_var(dirty, 5))
        recon_noise = vis_noise / 1000

        for i, curr_lambda in enumerate(lambdas):
            curr_output_name = path + "lambda_" + str(curr_lambda) + ".fits"
            os.system("julia ../julia_rascil_scripts/make_multistep.jl " + str(curr_lambda) + " " + tmp_psf_fn + \
                      " " + tmp_dirty_fn + " " + tmp_low_fn + " " + str(wavelet_type_idx) + " " + str(niter) + " " + \
                      str(recon_noise) + " " + str(vis_noise) + " " + str(cut_center) + " " + str(cut_hw) + " " + curr_output_name)
            recon = fits.open(curr_output_name)[0].data
            snrs[i] = compute_snr(gt, recon)
        
        write_to_csv(snrs, path + "snr.dat")

    snrs = read_csv(path + "snr.dat")
    snrs = [float(x) for x in snrs]

    return snrs

def banded_psnrs(path, gt, bands, nmaj_iter):
    data = []
    deconv_images = []

    prev_deconv = None
    for maj_iter in range(0, nmaj_iter):
        filename = path + "/deconv_iteration_" + str(maj_iter) + "_channel_0.fits"
        deconv = readFits(filename)
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
            psnr = compute_snr(curr_gt, filtered_deconv)
            curr_data.append(psnr)

        data.append(curr_data)
    
    return data

def banded_psnrs_lowresdirty(path, bands, cut, hw, nmaj_iter):
    deconv_images = []

    psf = readFits(path + "/psf_channel0_iteration0.fits")
    lowres_constraint = readFits(path + "/lowres_iteration_0_channel_0.fits")
    dirty = readFits(path + "/dirty_iteration_0_channel_0.fits")

    prev_deconv = None
    for maj_iter in range(0, nmaj_iter):
        deconv_filename = path + "/deconv_iteration_" + str(maj_iter) + "_channel_0.fits"
        deconv = readFits(deconv_filename)
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
                psnr = compute_snr(low_gt, filtered_deconv)
                curr_data_low.append(psnr)
                
            if band[0] >= (cut - hw):
                convolved_deconv = signal.fftconvolve(filtered_deconv, psf, mode='same')
                psnr = compute_snr(high_gt, convolved_deconv)
                curr_data_high.append(psnr)
        
        if band[1] <= (cut + hw):
            data_low.append(curr_data_low)
        
        if band[0] >= (cut - hw):
            data_high.append(curr_data_high)

    return data_low, data_high


def deconvolve(dirty, psf, lowres, niter, wavelet_type_idx, curr_maj_iter, initial_lambda, lambda_mul, cut_center=20, cut_halfwidth=5, variance_window=5, recon_variance_factor=1, lowin=None, vis_variance=None, recon_variance=None):
    curr_lambda = initial_lambda * (lambda_mul ** curr_maj_iter)
    curr_lambda *= numpy.linalg.norm(dirty)

    tmp_psf_name = "tmp_psf.fits"
    tmp_res_name = "tmp_residual.fits"
    tmp_lowin_name = "tmp_lowin.fits"
    tmp_output_name = "tmp_output.fits"

    helpers.write_nparr_to_fits(psf, tmp_psf_name)
    helpers.write_nparr_to_fits(dirty, tmp_res_name)

    #low resolution step if no constraint
    if lowin is None:
        os.system("julia ../julia_rascil_scripts/make_lowres.jl " + str(curr_lambda) + " " + tmp_psf_name + " " + tmp_res_name + " " + str(wavelet_type_idx) + " " + str(niter) + " " + tmp_output_name)
    else:
        #auto calculate sigma and eta if none provided
        if vis_variance is None:
            vis_variance = numpy.mean(helpers.compute_windowed_var(dirty, variance_window))
            recon_variance = vis_variance / recon_variance_factor

        helpers.write_nparr_to_fits(lowres, tmp_lowin_name)

        os.system("julia ../julia_rascil_scripts/make_multistep.jl " + str(curr_lambda) + " " + tmp_psf_name + " " + tmp_res_name + " " + tmp_lowin_name + " " + str(wavelet_type_idx) + " " + str(niter) + " " + \
                str(recon_variance) + " " + str(vis_variance) + " " + str(cut_center) + " " + str(cut_halfwidth) + " " + tmp_output_name)

    deconvolved = helpers.readFits(tmp_output_name)

    return deconvolved

def create_empty_image(vis, npixel, cellsize):
    return create_image_from_visibility(vis, npixel=npixel, cellsize=cellsize, polarisation_frame=vis.visibility_acc.polarisation_frame)

def generate_fake_baselines(nbaselines):
    for i in range(0, nbaselines):
        yield 0, 0

#copy of rascil's create_visibility_from_ms function, but with a fix to the channel allocation so that it works with large datasets
#also handy as the latest version of rascil has removed this function
def create_visibility_from_ms2(
    msname,
    channum=None,
    start_chan=None,
    end_chan=None,
    ack=False,
    datacolumn="DATA",
    selected_sources=None,
    selected_dds=None,
    average_channels=False,
):
    """Minimal MS to Visibility converter

    The MS format is much more general than the RASCIL Visibility so we cut many corners.
    This requires casacore to be installed. If not an exception ModuleNotFoundError is raised.

    Creates a list of Visibility's, split by field and spectral window

    Reading of a subset of channels is possible using either start_chan and end_chan or channnum. Using start_chan
    and end_chan is preferred since it only reads the channels required. Channum is more flexible and can be used to
    read a random list of channels.

    :param msname: File name of MS
    :param channum: range of channels e.g. range(17,32), default is None meaning all
    :param start_chan: Starting channel to read
    :param end_chan: End channel to read
    :param ack: Ask casacore to acknowledge each table operation
    :param datacolumn: MS data column to read DATA, CORRECTED_DATA, or MODEL_DATA
    :param selected_sources: Sources to select
    :param selected_dds: Data descriptors to select
    :param average_channels: Average all channels read
    :return: List of Visibility

    For example::

        selected_sources = ['1302+5748', '1252+5634']
        bvis_list = create_visibility_from_ms('../../data/3C277.1_avg.ms', datacolumn='CORRECTED_DATA',
                                           selected_sources=selected_sources)
        sources = numpy.unique([bv.source for bv in bvis_list])
        print(sources)
        ['1252+5634' '1302+5748']

    """
    try:
        from casacore.tables import table  # pylint: disable=import-error
    except ModuleNotFoundError:
        raise ModuleNotFoundError("casacore is not installed")
    try:
        from rascil.processing_components.visibility import msv2
    except ModuleNotFoundError:
        raise ModuleNotFoundError("cannot import msv2")

    tab = table(msname, ack=ack)

    if selected_sources is None:
        fields = numpy.unique(tab.getcol("FIELD_ID"))
    else:
        fieldtab = table("%s/FIELD" % msname, ack=False)
        sources = fieldtab.getcol("NAME")
        fields = list()
        for field, source in enumerate(sources):
            if source in selected_sources:
                fields.append(field)
        assert len(fields) > 0, "No sources selected"

    if selected_dds is None:
        dds = numpy.unique(tab.getcol("DATA_DESC_ID"))
    else:
        dds = selected_dds

    total_vis = 0

    vis_list = list()
    for field in fields:
        ftab = table(msname, ack=ack).query("FIELD_ID==%d" % field, style="")
        assert ftab.nrows() > 0, "Empty selection for FIELD_ID=%d" % (field)
        for dd in dds:
            # Now get info from the subtables
            ddtab = table("%s/DATA_DESCRIPTION" % msname, ack=False)
            spwid = ddtab.getcol("SPECTRAL_WINDOW_ID")[dd]
            polid = ddtab.getcol("POLARIZATION_ID")[dd]
            ddtab.close()

            meta = {"MSV2": {"FIELD_ID": field, "DATA_DESC_ID": dd}}
            ms = ftab.query("DATA_DESC_ID==%d" % dd, style="")
            assert (
                ms.nrows() > 0
            ), "Empty selection for FIELD_ID=%d and DATA_DESC_ID=%d" % (field, dd)
            # The TIME column has descriptor:
            # {'valueType': 'double', 'dataManagerType': 'IncrementalStMan', 'dataManagerGroup': 'TIME',
            # 'option': 0, 'maxlen': 0, 'comment': 'Modified Julian Day',
            # 'keywords': {'QuantumUnits': ['s'], 'MEASINFO': {'type': 'epoch', 'Ref': 'UTC'}}}
            otime = ms.getcol("TIME")
            datacol = ms.getcol(datacolumn, nrow=1)
            datacol_shape = list(datacol.shape)
            channels = datacol.shape[-2]
            if channum is None:
                if start_chan is not None and end_chan is not None:
                    try:
                        blc = [start_chan, 0]
                        trc = [end_chan, datacol_shape[-1] - 1]
                        channum = range(start_chan, end_chan + 1)
                        ms_vis = ms.getcolslice(datacolumn, blc=blc, trc=trc)
                        ms_flags = ms.getcolslice("FLAG", blc=blc, trc=trc)
                        ms_weight = ms.getcol("WEIGHT")

                    except IndexError:
                        raise IndexError("channel number exceeds max. within ms")

                else:
                    try:
                        channum = range(channels)
                        ms_vis = ms.getcol(datacolumn)[:, channum, :]
                        ms_weight = ms.getcol("WEIGHT")
                        ms_flags = ms.getcol("FLAG")[:, channum, :]
                        channum = range(channels)
                    except IndexError:
                        raise IndexError("channel number exceeds max. within ms")
            else:
                channum = range(channels)
                try:
                    ms_vis = ms.getcol(datacolumn)[:, channum, :]
                    ms_flags = ms.getcol("FLAG")[:, channum, :]
                    ms_weight = ms.getcol("WEIGHT")[:, :]
                except IndexError:
                    raise IndexError("channel number exceeds max. within ms")

            if average_channels:
                weight = ms_weight[:, numpy.newaxis, :] * (1.0 - ms_flags)
                ms_vis = numpy.sum(weight * ms_vis, axis=-2)[..., numpy.newaxis, :]
                sumwt = numpy.sum(weight, axis=-2)[..., numpy.newaxis, :]
                ms_vis[sumwt > 0.0] = ms_vis[sumwt > 0] / sumwt[sumwt > 0.0]
                ms_vis[sumwt <= 0.0] = 0.0 + 0.0j
                ms_flags = sumwt
                ms_flags[ms_flags <= 0.0] = 1.0
                ms_flags[ms_flags > 0.0] = 0.0

            total_vis += ms_vis.shape[0]

            uvw = -1 * ms.getcol("UVW")
            antenna1 = ms.getcol("ANTENNA1")
            antenna2 = ms.getcol("ANTENNA2")
            integration_time = ms.getcol("INTERVAL")

            time = otime - integration_time / 2.0

            start_time = numpy.min(time) / 86400.0
            end_time = numpy.max(time) / 86400.0

            spwtab = table("%s/SPECTRAL_WINDOW" % msname, ack=False)
            cfrequency = numpy.array(spwtab.getcol("CHAN_FREQ")[spwid][channum])
            cchannel_bandwidth = numpy.array(
                spwtab.getcol("CHAN_WIDTH")[spwid][channum]
            )
            nchan = cfrequency.shape[0]
            if average_channels:
                cfrequency = numpy.array([numpy.average(cfrequency)])
                cchannel_bandwidth = numpy.array([numpy.sum(cchannel_bandwidth)])
                nchan = cfrequency.shape[0]
            else:
                nchan = len(channum)

            # Get polarisation info
            poltab = table("%s/POLARIZATION" % msname, ack=False)
            corr_type = poltab.getcol("CORR_TYPE")[polid]
            corr_type = sorted(corr_type)
            # These correspond to the CASA Stokes enumerations
            if numpy.array_equal(corr_type, [1, 2, 3, 4]):
                polarisation_frame = PolarisationFrame("stokesIQUV")
                npol = 4
            elif numpy.array_equal(corr_type, [1, 2]):
                polarisation_frame = PolarisationFrame("stokesIQ")
                npol = 2
            elif numpy.array_equal(corr_type, [1, 4]):
                polarisation_frame = PolarisationFrame("stokesIV")
                npol = 2
            elif numpy.array_equal(corr_type, [5, 6, 7, 8]):
                polarisation_frame = PolarisationFrame("circular")
                npol = 4
            elif numpy.array_equal(corr_type, [5, 8]):
                polarisation_frame = PolarisationFrame("circularnp")
                npol = 2
            elif numpy.array_equal(corr_type, [9, 10, 11, 12]):
                polarisation_frame = PolarisationFrame("linear")
                npol = 4
            elif numpy.array_equal(corr_type, [9, 12]):
                polarisation_frame = PolarisationFrame("linearnp")
                npol = 2
            elif numpy.array_equal(corr_type, [9]) or numpy.array_equal(corr_type, [1]):
                npol = 1
                polarisation_frame = PolarisationFrame("stokesI")
            else:
                raise KeyError("Polarisation not understood: %s" % str(corr_type))

            # Get configuration
            anttab = table("%s/ANTENNA" % msname, ack=False)
            names = numpy.array(anttab.getcol("NAME"))

            ant_map = list()
            actual = 0
            # This assumes that the names are actually filled in!
            for i, name in enumerate(names):
                if name != "":
                    ant_map.append(actual)
                    actual += 1
                else:
                    ant_map.append(-1)
            # assert actual > 0, "Dish/station names are all blank - cannot load"
            if actual == 0:
                ant_map = list(range(len(names)))
                names = numpy.repeat("No name", len(names))

            mount = numpy.array(anttab.getcol("MOUNT"))[names != ""]
            diameter = numpy.array(anttab.getcol("DISH_DIAMETER"))[names != ""]
            xyz = numpy.array(anttab.getcol("POSITION"))[names != ""]
            offset = numpy.array(anttab.getcol("OFFSET"))[names != ""]
            stations = numpy.array(anttab.getcol("STATION"))[names != ""]
            names = numpy.array(anttab.getcol("NAME"))[names != ""]
            nants = len(names)

            antenna1 = list(map(lambda i: ant_map[i], antenna1))
            antenna2 = list(map(lambda i: ant_map[i], antenna2))

            baselines = pandas.MultiIndex.from_tuples(
                generate_fake_baselines(ms_vis.shape[0]), names=("antenna1", "antenna2")
                #generate_baselines(nants), names=("antenna1", "antenna2")
                #generate_baselines(1), names=("antenna1", "antenna2")
            )

            #nbaselines = len(baselines)
            nbaselines = ms_vis.shape[0]

            location = EarthLocation(
                x=Quantity(xyz[0][0], "m"),
                y=Quantity(xyz[0][1], "m"),
                z=Quantity(xyz[0][2], "m"),
            )

            configuration = Configuration.constructor(
                name="",
                location=location,
                names=names,
                xyz=xyz,
                mount=mount,
                frame="ITRF",
                receptor_frame=ReceptorFrame("linear"),
                diameter=diameter,
                offset=offset,
                stations=stations,
            )
            # Get phasecentres
            fieldtab = table("%s/FIELD" % msname, ack=False)
            pc = fieldtab.getcol("PHASE_DIR")[field, 0, :]
            source = fieldtab.getcol("NAME")[field]
            phasecentre = SkyCoord(
                ra=pc[0] * u.rad, dec=pc[1] * u.rad, frame="icrs", equinox="J2000"
            )

            time_index_row = numpy.zeros_like(time, dtype="int")
            time_last = time[0]
            time_index = 0
            for row, _ in enumerate(time):
                if time[row] > time_last + 0.5 * integration_time[row]:
                    assert (
                        time[row] > time_last
                    ), "MS is not time-sorted - cannot convert"
                    time_index += 1
                    time_last = time[row]
                time_index_row[row] = time_index

            ntimes = time_index + 1

            assert ntimes == len(
                numpy.unique(time_index_row)
            ), "Error in finding data times"

            #ntimes = ms_vis.shape[0]
            ntimes = 1

            bv_times = numpy.zeros([ntimes])
            bv_vis = numpy.zeros([ntimes, nbaselines, nchan, npol]).astype("complex")
            bv_flags = numpy.zeros([ntimes, nbaselines, nchan, npol]).astype("int")
            bv_weight = numpy.zeros([ntimes, nbaselines, nchan, npol])
            bv_uvw = numpy.zeros([ntimes, nbaselines, 3])
            bv_integration_time = numpy.zeros([ntimes])

            for row, _ in enumerate(time):
                #ibaseline = baselines.get_loc((antenna1[row], antenna2[row]))
                #ibaseline = 0
                ibaseline = row

                #time_index = time_index_row[row]
                #time_index = row
                time_index = 0

                bv_times[time_index] = time[row]
                bv_vis[time_index, ibaseline, ...] = ms_vis[row, ...]
                bv_flags[time_index, ibaseline, ...][
                    ms_flags[row, ...].astype("bool")
                ] = 1
                bv_weight[time_index, ibaseline, :, ...] = ms_weight[
                    row, numpy.newaxis, ...
                ]
                bv_uvw[time_index, ibaseline, :] = uvw[row, :]
                bv_integration_time[time_index] = integration_time[row]

            vis_list.append(
                Visibility.constructor(
                    uvw=bv_uvw,
                    baselines=baselines,
                    time=bv_times,
                    frequency=cfrequency,
                    channel_bandwidth=cchannel_bandwidth,
                    vis=bv_vis,
                    flags=bv_flags,
                    weight=bv_weight,
                    integration_time=bv_integration_time,
                    configuration=configuration,
                    phasecentre=phasecentre,
                    polarisation_frame=polarisation_frame,
                    source=source,
                    meta=meta,
                )
            )
        tab.close()

    return vis_list, total_vis


def compute_residual_bychannel(sky_estimate, ms_name, npixel, cellsize, weighting, robustness, weight_grid, channel_start, channel_end, data_descriptors, algorithm='ng'):
    final_residual = None

    channels = range(channel_start, channel_end + 1)

    for dd in data_descriptors:
        for curr_channel in channels:
            [measured_vis], _ = create_visibility_from_ms2(ms_name, start_chan=curr_channel, end_chan=curr_channel, selected_dds=[dd])
            measured_vis = convert_visibility_to_stokesI(measured_vis)
            measured_vis = griddata_visibility_reweight(measured_vis, weight_grid[0], weighting=weighting, robustness=robustness, sumwt=weight_grid[1])

            estimated_vis = measured_vis.copy(deep=True)
            estimated_vis = predict_ng(estimated_vis, sky_estimate, context=algorithm)
            residual_vis = subtract_visibility(measured_vis, estimated_vis)

            if final_residual is None:
                final_residual = create_empty_image(measured_vis, npixel, cellsize)

            channel_residual, sumwt = invert_ng(residual_vis, final_residual, context=algorithm)

            final_residual = add_to_image(final_residual, channel_residual.pixels.data[0,0,:,:])

            gc.collect()

    gc.collect()

    return final_residual

def compute_jackknifed_residual_bychannel(sky_estimate, ms_name, npixel, cellsize, weighting, robustness, weight_grid, channel_start, channel_end, data_descriptors, algorithm='ng'):
    final_residual = None

    channels = range(channel_start, channel_end + 1)
    rng = numpy.random.default_rng(42)

    for dd in data_descriptors:
        for curr_channel in channels:
            [measured_vis], _ = create_visibility_from_ms2(ms_name, start_chan=curr_channel, end_chan=curr_channel, selected_dds=[dd])
            measured_vis = convert_visibility_to_stokesI(measured_vis)

            for (i, j, k, l), vis in numpy.ndenumerate(measured_vis.vis):
                if rng.random() > 0.5:
                    measured_vis.vis.data[i, j, k, l] *= -1

            measured_vis = griddata_visibility_reweight(measured_vis, weight_grid[0], weighting=weighting, robustness=robustness, sumwt=weight_grid[1])

            estimated_vis = measured_vis.copy(deep=True)
            estimated_vis = predict_ng(estimated_vis, sky_estimate, context=algorithm)
            residual_vis = subtract_visibility(measured_vis, estimated_vis)

            if final_residual is None:
                final_residual = create_empty_image(measured_vis, npixel, cellsize)

            channel_residual, sumwt = invert_ng(residual_vis, final_residual, context=algorithm)

            final_residual = add_to_image(final_residual, channel_residual.pixels.data[0,0,:,:])

            gc.collect()

    gc.collect()

    return final_residual

def compute_weights_griddata_by_channel(ms_name, npixel, cellsize, channel_start, channel_end, data_descriptors):
    total_grid = None
    model = None

    channels = range(channel_start, channel_end + 1)

    for dd in data_descriptors:
        for curr_channel in channels:
            [vis], num_vis = create_visibility_from_ms2(ms_name, start_chan=curr_channel, end_chan=curr_channel, selected_dds=[dd])
            vis = convert_visibility_to_stokesI(vis)

            if model is None:
                model = create_image_from_visibility(vis, cellsize=cellsize, npixel=npixel, polarisation_frame=vis.visibility_acc.polarisation_frame)

            curr_grid_weights = create_griddata_from_image(model, polarisation_frame=model.image_acc.polarisation_frame)
            curr_grid_weights = grid_visibility_weight_to_griddata(vis, curr_grid_weights)
            
            if total_grid is None:
                total_grid = curr_grid_weights
            else:
                total_grid = griddata_merge_weights([total_grid, curr_grid_weights])

            gc.collect()

    gc.collect()

    return total_grid, model



def plot_snr_across_bands(title, path, cuts, bands, datasets, colours, linestyles, labels):
    for i, dataset in enumerate(datasets):
        gt = readFits("../data/" + dataset + "_full_gt.fits")
        for j, cut in enumerate(cuts):
            reconstructed = readFits(path + dataset + "_" + str(cut) + ".fits")
            snrs = []
            for band in bands:
                curr_band_filter = bandpass(reconstructed, band[0], band[1])
                
                curr_recon = signal.fftconvolve(reconstructed, curr_band_filter, mode='same')
                curr_gt = signal.fftconvolve(gt, curr_band_filter, mode='same')
                
                curr_snr = compute_snr(curr_gt, curr_recon)
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
        gt = readFits("../data/" + dataset + "_full_gt.fits")

        reconstructed = readFits(path + dataset + "_55.fits")
        snrs = []
        for band in bands:
            curr_band_filter = bandpass(reconstructed, band[0], band[1])
            
            curr_recon = signal.fftconvolve(reconstructed, curr_band_filter, mode='same')
            curr_gt = signal.fftconvolve(gt, curr_band_filter, mode='same')
            
            curr_snr = compute_snr(curr_gt, curr_recon)
            snrs.append(curr_snr)
            
        plt.plot([x[0] for x in bands], snrs, c=colours[i], linestyle=linestyles[0], \
                 label=labels[i], lw=3)

    plt.legend(loc="upper right", fontsize=10)

    plt.title(title)

    plt.xlabel("Frequency band", fontsize=12)
    plt.ylabel("SNR", fontsize=12)
    plt.show()

def plot_snr_across_bands_and_cases(dataset, cut, dataset_title, bands, cases, case_paths, colours):
    gt = readFits("../data/" + dataset + "_full_gt.fits")

    for i, case in enumerate(cases):
        filename = case_paths[i] + dataset + "_" + str(cut) + ".fits" if case.lower() != "single step" else case_paths[i] + dataset + "_55.fits"
        reconstructed = readFits(filename)
        snrs = []
        for band in bands:
            curr_band_filter = bandpass(reconstructed, band[0], band[1])
            
            curr_recon = signal.fftconvolve(reconstructed, curr_band_filter, mode='same')
            curr_gt = signal.fftconvolve(gt, curr_band_filter, mode='same')
            
            curr_snr = compute_snr(curr_gt, curr_recon)
            snrs.append(curr_snr)
            
        plt.plot([x[0] for x in bands], snrs, c=colours[i], label=cases[i], lw=3)

    plt.legend(loc="upper right", fontsize=10)

    plt.title("SNRs by strategies for " + dataset_title)

    plt.xlabel("Frequency band", fontsize=12)
    plt.ylabel("SNR", fontsize=12)
    plt.show()

def image_histogram_equalization(image, number_bins=256):
    image_histogram, bins = numpy.histogram(image.flatten(), number_bins, density=True)
    cdf = image_histogram.cumsum() # cumulative distribution function
    cdf = (number_bins-1) * cdf / cdf[-1] # normalize
    image_equalized = numpy.interp(image.flatten(), bins[:-1], cdf)

    return image_equalized.reshape(image.shape)

def shift_img(img, x_offset, y_offset, zero_centered = True):
    x_dim = img.shape[0]
    y_dim = img.shape[1]

    fimg = numpy.fft.fft2(numpy.fft.ifftshift(img)) if zero_centered else numpy.fft.fft2(img)

    u = numpy.fft.fftfreq(x_dim)
    v = numpy.fft.fftfreq(y_dim)
    U, V = numpy.meshgrid(u, v)
    fsimg = fimg * numpy.exp(-2j * numpy.pi * (x_offset * U + y_offset * V))

    return numpy.real(numpy.fft.fftshift(numpy.fft.ifft2(fsimg)) if zero_centered else numpy.fft.ifft2(fsimg))

def add_to_image(image, nparr):
    image.pixels.data[0,0,:,:] += nparr

    return image