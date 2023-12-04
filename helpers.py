import astropy
from astropy.io import fits

import matplotlib.pyplot as plt
import matplotlib
import numpy

import csv
from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy import signal

from ska_sdp_func_python.visibility import subtract_visibility

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
    return 20 * numpy.log10(numpy.linalg.norm(gt) / numpy.linalg.norm(gt-recon))

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

def plotNImages(images, names, cmap, same_scale=False):
    num_images = len(images)

    fig, axes = plt.subplots(1, num_images)

    im = None

    vmin = 999999999999
    vmax = -999999999999

    if same_scale:
        for img in images:
            vmin = min(vmin, numpy.min(img))
            vmax = max(vmax, numpy.max(img))

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
            cax = divider.append_axes("bottom", size="5%", pad=0.25)

            cb = fig.colorbar(im, orientation='horizontal', cax=cax)
            cb.ax.locator_params(nbins=5)
        else:
            axes.set_title(names[i])
            if same_scale:
                im = axes.imshow(img, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)
            else:
                im = axes.imshow(img, cmap=cmap, origin='lower')

            divider = make_axes_locatable(axes)
            cax = divider.append_axes("bottom", size="5%", pad=0.25)

            cb = fig.colorbar(im, orientation='horizontal', cax=cax)
            cb.ax.locator_params(nbins=5)

    plt.show()

def plotGDP(gt, dirty, psf, cmap):
    plotNImages([gt, dirty, psf], ["True Sky", "Dirty Image", "PSF"], cmap)

def plotSNRvsSSIM(lambdas, path, snr_idx, ssim_idx, gt, cmap):
    snr_fn = path + "lambda_" + str(lambdas[snr_idx]) + ".fits"
    ssim_fn = path + "lambda_" + str(lambdas[ssim_idx]) + ".fits"

    snrfile = readFits(snr_fn)
    ssimfile = readFits(ssim_fn)

    err_snr = numpy.abs(gt - snrfile)
    err_ssim = abs(gt - ssimfile)

    snr_title = "Best SNR Image $\lambda = " + "{:.2f}".format(lambdas[snr_idx]) + "$"
    ssim_title = "Best SSIM Image $\lambda = " + "{:.2f}".format(lambdas[ssim_idx]) + "$"

    plotNImages([snrfile, ssimfile], [snr_title, ssim_title], cmap)
    plotNImages([err_snr, err_ssim], ["SNR Absolute Diff", "SSIM Absolute Diff"], cmap, same_scale=True)


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

# def conv(f, g, p1=-1, p2=-1) :                                                                                                                                                           
#     if p1 > 0 and p2 > 0:                                                                                   
#         p1, p2 = (numpy.r_[f.shape]-g.shape).astype(int)//2                                                                                                                              
#         print(str(p1) + " " + str(p2))
#         gpad = numpy.pad(g,((p1,p1),(p2,p2)),mode='edge')                                                                                                                                                                                                                                          
    
#     print(gpad.shape)
#     print(f.shape)
#     gpad = numpy.fft.ifftshift(gpad)

#     FG = numpy.fft.fft2(f) * numpy.conj(numpy.fft.fft2(gpad))
#     return numpy.real(numpy.fft.ifft2(FG)) 

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
