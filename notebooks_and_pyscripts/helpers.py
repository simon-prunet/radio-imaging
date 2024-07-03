import astropy
from astropy.io import fits

import matplotlib.pyplot as plt
import matplotlib
import numpy

import csv
from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy import signal

from ska_sdp_func_python.visibility import subtract_visibility
import os

def write_to_csv(data, filename):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(data)

def write_nparr_to_fits(data, filename):
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(filename, overwrite=True)
    hdulist.close()

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
                im = axes[i].imshow(img, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)
            else:
                im = axes[i].imshow(img, cmap=cmap, origin='lower')

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
                im = axes.imshow(img, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)
            else:
                im = axes.imshow(img, cmap=cmap, origin='lower')

            divider = make_axes_locatable(axes)
            cax = divider.append_axes(colorbar_location, size="5%", pad=0.25)

            cb = fig.colorbar(im, orientation='horizontal', cax=cax)
            #cb.formatter.set_powerlimits((-10, 10))
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
    plotNImages([gt, dirty, psf], ["True Sky", "Dirty Image", "PSF"], cmap)

def plotSNRvsSSIM(lambdas, path, snr_idx, gt, cmap, same_scale = False):
    snr_fn = path + "lambda_" + str(lambdas[snr_idx]) + ".fits"
    snrfile = readFits(snr_fn)
    err_snr = numpy.abs(gt - snrfile)
    snr_title = "Best SNR Image $\lambda = " + "{:.2f}".format(lambdas[snr_idx]) + "$"

    plotNImages([snrfile, err_snr], [snr_title, "SNR absolute error"], cmap, same_scale = same_scale)


def read_csv(filename):
    data = []
    with open(filename, newline='') as file:
        reader = csv.reader(file, delimiter=',')
        for row in reader:
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

def lambdatests_allvis(lambdas, path, dirty, psf, gt, wavelet_type_idx, niter, runtests, lp_filter=None):
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