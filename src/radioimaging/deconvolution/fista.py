#!/usr/bin/env python

"""\
FISTA implementation for deconvolving multi-resolution radio images 
""" 


import numpy as np
import pywt
from joblib import Parallel, delayed
import scipy

from radioimaging.util import util
from radioimaging.util.util import convolve2d

EPSILON = 1e-6

#computes adjoint of a convolution for some given kernel, as with convolve2d, assumes circular convolution
def adj(kernel):
    """
    adj computes adjoint of kernel

    :kernel: input kernel
    :return: kernel adjoint
    """
    shift = (1 if kernel.shape[0] % 2 == 0 else 0, 1 if kernel.shape[1] % 2 == 0 else 0)
    return np.roll(np.flip(kernel), shift=shift, axis=(0, 1)) 


def pseudoinv(kernel, const=1):
    """
    pseudoinv computes the pseudoinverse of a kernel

    :kernel: input kernel
    :const: multiplier to multiply singular values by, used for variance estimations
    :return: kernel pseudoinverse
    """
    fkernel = np.fft.fft2(kernel)
    mask = np.abs(fkernel) > EPSILON

    fkernel[~mask] = 0
    fkernel[mask] = 1 / (fkernel[mask] * const)

    return np.fft.ifft2(fkernel)

class Filter2D:
    def __init__(self, kernel):
        self.kernel = kernel
        self.adj = adj(kernel)
        self.kernel2 = None

    def squared(self, linear_conv=False):
        if self.kernel2 is None:
            self.kernel2 = util.convolve2d(self.kernel, self.adj, linear=linear_conv)

        return self.kernel2

def get_max_wavelet_level(img, wavelets):
    min_size = min(img.shape)

    levels = []

    for wavelet in wavelets:
        levels.append(pywt.dwt_max_level(data_len=min_size, filter_len=wavelet.dec_len))

    return min(levels)

def dwt_single(img, wavelet, level):
    """
    dwt_single computes the wavelet transform of an image for a single wavelet dictionary

    :img: input image
    :wavelet: dictionary
    :return: wavelet coefficients with slices
    """
    return pywt.coeffs_to_array(pywt.wavedec2(img, wavelet=wavelet, mode="periodization", level=level))

def idwt_single(coeffs, slice, wavelet):
    """
    idwt_single computes the image from a given set of wavelet coefficients for some single wavelet dictionary

    :coeffs: wavelet coefficients
    :wavelet: dictionary
    :slice: wavelet slice, provided when computing the wavelet transform
    :return: image
    """
    return pywt.waverec2(pywt.array_to_coeffs(coeffs, slice, output_format="wavedec2"), wavelet=wavelet, mode="periodization")

def dwt(img, wavelets, parallel=False):
    """
    dwt computes the wavelet transform for an overredundant wavelet dictionary

    :img: input image
    :wavelets: list of wavelet dictionaries
    :parallel: multithread, currently is slower than the normal version so leave this false
    :return: wavelet coefficients and slices
    """
    coeffs = np.zeros((len(wavelets),) + img.shape)
    slices = [None for x in range(len(wavelets))]

    level = get_max_wavelet_level(img, wavelets)

    if parallel:
        results = Parallel(n_jobs=len(wavelets))(delayed(dwt_single)(img, w, level) for w in wavelets)

        for i, result in enumerate(results):
            coeffs[i,:,:] = result[0]
            slices[i] = result[1]
    else:
        for i, wavelet in enumerate(wavelets):
            curr_wavelet_coeffs = dwt_single(img, wavelet, level)
            coeffs[i,:,:] = curr_wavelet_coeffs[0]
            slices[i] = curr_wavelet_coeffs[1]

    return coeffs, slices

def idwt(coeffs, slices, wavelets, parallel=False):
    """
    idwt computes an image from an overredundant wavelet dictionary

    :coeffs: wavelet coefficients
    :slices: wavelet slices, obtained from dwt
    :wavelets: list of wavelet dictionaries
    :parallel: multithread, currently is slower than the normal version so leave this false
    :return: image
    """
    img = np.zeros(coeffs[0].shape)

    if parallel:
        results = Parallel(n_jobs=len(wavelets))(delayed(idwt_single)(coeffs[i], slices[i], w) for i, w in enumerate(wavelets))

        for i in results:
            img += i
    else:
        for i, wavelet in enumerate(wavelets):
            img += idwt_single(coeffs[i], slices[i], wavelet)

    return img


def compute_lambda_max(meas_term, wavelets, parallelize_wavelets=False):
    """
    compute_lambda_max computes lambda_max, the smallest lambda for where the solution is zero

    :meas_term: measurement term
    :wavelets: list of wavelet dictionaries
    :parallelize_wavelets: multithread wavelet transform, currently is slower than the normal version so leave this false
    :return: lambda_max
    """
    return np.max(np.abs(2*dwt(meas_term, wavelets, parallel=parallelize_wavelets)[0])) 


def compute_lambda_max_ext(dirty, psf, wavelets, filters=None, constraint_images=None, partition=None, parallelize_wavelets=False, linear_conv=False):
    """
    compute_lambda_max_ext computes lambda_max, the smallest lambda for where the solution is zero. This is meant for external use, thus rather than taking the measurement
    term, it instead takes the individual terms and filters and computes the measurement term

    :dirty: dirty image
    :psf: psf
    :wavelets: list of wavelet dictionaries
    :filters: list of filters for each partition
    :constraint_images: list of additional constraint images, corresponding to previously deconvolved images of each partition
    :partition: current partition
    :parallelize_wavelets: multithread wavelet transform, currently is slower than the normal version so leave this false
    :linear_conv: whether to use linear convolution
    :return: lambda_max
    """
    meas_term = None

    if filters is None:
        meas_term = convolve2d(dirty, psf.adj, linear=linear_conv)
    else:
        curr_filter = filters[partition]

        meas_term = convolve2d(convolve2d(dirty, curr_filter.squared(linear_conv=linear_conv), linear=linear_conv), psf.adj, linear=linear_conv)

        #add additional constraint fidelity terms
        for i, constraint in enumerate(constraint_images):
            if i == partition:
                continue

            curr_filter = filters[i]

            meas_term += convolve2d(constraint_images[i], curr_filter.squared(linear_conv=linear_conv), linear=linear_conv)

    return np.max(np.abs(2*dwt(meas_term, wavelets, parallel=parallelize_wavelets)[0])) 


def compute_step(psf, wavelets, orthowavelets, wavelet_slices, filters=None, partition=None, power_iterations=50, parallelize_wavelets=False, linear_conv=False):
    """
    compute_step computes the gradient step size (1/lipschitz). Problem is quadratic so uses largest eigenvalue, either from
    a closed form calculation when the operator is easily diagonalizable, or power iteration when not

    :psf: point spread function
    :wavelets: list of wavelet dictionaries
    :orthowavelets: whether individual wavelet dictionaries are orthogonal
    :wavelet_slices: wavelet slices obtained from dwt
    :filters: filters for different partitions, or None when deconvolving full resolution image
    :partition: current partition, or None when deconvolving full resolution image
    :power_iterations: number of iterations for power iteration
    :parallelize_wavelets: multithread wavelet transform, currently is slower than the normal version so leave this false
    :linear_conv: perform linear convolution, rather than circular
    :return: gradient step size
    """
    if orthowavelets:
        lipschitz = None

        if filters is None:
            lipschitz = np.max(np.abs(np.fft.fft2(psf.squared(linear_conv=linear_conv))))
        else:
            fpsf = np.fft.fft2(psf.kernel)
            ffilt = np.fft.fft2(filters[partition].kernel)
            diag = fpsf*fpsf*ffilt*ffilt
            
            for i, filt in enumerate(filters):
                if i == partition:
                    continue

                diag += np.fft.fft2(filters[i].squared(linear_conv=linear_conv))

            lipschitz = np.max(np.abs(diag))

        return 1 / (2 * len(wavelets) * lipschitz)
    #power iteration
    else:
        operator_kernel = None
        if filters is None:
            operator_kernel = psf.squared(linear_conv=linear_conv)
        else:
            operator_kernel = util.convolve2d(util.convolve2d(psf, curr_filter.squared(linear_conv=linear_conv), linear=linear_conv), psf.adj, linear=linear_conv)

        for i, filt in enumerate(filters):
            if i == partition:
                continue

            operator_kernel += filt.squared(linear_conv=linear_conv)

        img_dims = psf.kernel.shape
        alpha = np.random.rand((len(wavelets),) + img_dims)

        for i in range(power_iterations):
            alpha = dwt(util.convolve2d(idwt(alpha, wavelet_slices, wavelets, parallel=parallelize_wavelets), operator_kernel, linear=linear_conv), wavelets, parallel=parallelize_wavelets)[0]
            alpha = alpha / np.linalg.norm(alpha)

        alpha_new = dwt(util.convolve2d(idwt(alpha, wavelet_slices, wavelets, parallel=parallelize_wavelets), operator_kernel, linear=linear_conv), wavelets, parallel=parallelize_wavelets)[0]

        return np.dot(alpha.T, alpha_new)

def compute_step_cov(psf, cov_pi, wavelets, orthowavelets, wavelet_slices, power_iterations=50, parallelize_wavelets=False, linear_conv=False):
    """
    compute_step_cov computes the gradient step size (1/lipschitz) when taking into account estimated covariance. 
    Problem is quadratic so uses largest eigenvalue, either from a closed form calculation when the operator is easily 
    diagonalizable, or power iteration when not. This is for the full resolution deconvolution, and is still in an experimental phase

    :psf: point spread function
    :cov_pi: pseudo-inverse of covariance matrix (assumed to be diagonal in the fourier)
    :wavelets: list of wavelet dictionaries
    :orthowavelets: whether individual wavelet dictionaries are orthogonal
    :wavelet_slices: wavelet slices obtained from dwt
    :power_iterations: number of iterations for power iteration
    :parallelize_wavelets: multithread wavelet transform, currently is slower than the normal version so leave this false
    :linear_conv: perform linear convolution, rather than circular
    :return: gradient step size
    """

    #analytic
    if orthowavelets:
        lipschitz = np.max(np.abs(np.fft.fft2(util.convolve2d(util.convolve2d(psf.kernel, cov_pi, linear=linear_conv), psf.adj, linear=linear_conv))))

        return 1 / (2 * len(wavelets) * lipschitz)
    #power iteration
    else:
        operator_kernel = util.convolve2d(util.convolve2d(psf.kernel, cov_pi, linear=linear_conv), psf.adj, linear=linear_conv)

        img_dims = psf.kernel.shape
        alpha = np.random.rand((len(wavelets),) + img_dims)

        for i in range(power_iterations):
            alpha = dwt(util.convolve2d(idwt(alpha, wavelet_slices, wavelets, parallel=parallelize_wavelets), operator_kernel, linear=linear_conv), wavelets, parallel=parallelize_wavelets)[0]
            alpha = alpha / np.linalg.norm(alpha)

        alpha_new = dwt(util.convolve2d(idwt(alpha, wavelet_slices, wavelets, parallel=parallelize_wavelets), operator_kernel, linear=linear_conv), wavelets, parallel=parallelize_wavelets)[0]

        return np.dot(alpha.T, alpha_new)

def soft_thresh(coeffs, step):
    """
    soft_thresh performs the soft-thresholding operator on some coefficients

    :coeffs: coefficients
    :step: soft thresholding step
    :return: soft thresholded coefficients
    """
    return np.sign(coeffs) * np.maximum(np.abs(coeffs) - step, 0)

def cost(model, coeffs, psf, dirty, lambd, filters=None, constraint_images=None, partition=None, linear_conv=False):
    """
    compute cost function given some specific image and its wavelet coefficients, needed for the backtracking step of monotone fista
    """
    constraint_term = 0
    
    vis_term_img = dirty - util.convolve2d(model, psf.kernel, linear=False)

    if filters is not None:
        vis_term_img = util.convolve2d(vis_term_img, filters[partition].kernel, linear=False)

        for i, filt in enumerate(filters):
            if i == partition:
                continue

            curr_constraint_term = np.linalg.norm(util.convolve2d(constraint_images[i] - model, filt.kernel, linear=False))
            constraint_term += curr_constraint_term * curr_constraint_term

    vis_term = np.linalg.norm(vis_term_img)

    l1_term = np.sum(np.abs(coeffs)) * lambd

    return vis_term*vis_term + constraint_term + l1_term

def fista(psf, dirty, reg_param, wavelets, niter, orthowavelets=None, filters=None, constraint_images=None, partition=None, lambda_max=None, parallelize_wavelets=False, linear_conv=False):
    """
    fista uses fista to deconvolve an image using l1 regularization. This is for the multi-partition case, and thus uses filters and
    constraint images as additional constraints to the partial resolution input image

    :psf: point spread function
    :dirty: dirty image
    :reg_param: multiplier for lambda_max
    :wavelets: list of wavelet dictionaries
    :niter: number of fista iterations
    :orthowavelets: whether individual wavelet dictionaries are orthogonal
    :filters: list of filters for each partition, leave as none for full resolution case
    :constraint_images: list of constraint images for each partition, leave as none for full resolution case
    :lambda_max: a hard set lambda_max, leave as None if it is to be calculated
    :parallelize_wavelets: multithread wavelet transform, currently is slower than the normal version so leave this false
    :linear_conv: perform linear convolution, rather than circular
    :return: deconvolved image
    """

    img_dims = psf.kernel.shape
    coeff_dims = (2 ** int(np.ceil(np.log2(img_dims[0]))), 2 ** int(np.ceil(np.log2(img_dims[1]))))
    diff_dims = (coeff_dims[0] - img_dims[0], coeff_dims[1] - img_dims[1])
    #padding = ((diff_dims[0]//2, diff_dims[0]//2), (diff_dims[0]//2, diff_dims[0]//2))
    padding = ((0, diff_dims[0]), (0, diff_dims[1]))
    padding_psf = ((diff_dims[0]//2, diff_dims[0]//2), (diff_dims[1]//2, diff_dims[1]//2))

    #pad psf and dirty to next closest power of 2 if they are originally not, this is because pywavelets does not give the same number of coefficients
    #per dictionary if image is not this size
    if diff_dims[0] > 0 or diff_dims[1] > 0:
        psf = Filter2D(np.pad(psf.kernel, padding_psf, mode="constant"))
        dirty = np.pad(dirty, padding, mode="wrap")

        if filters is not None:
            for i, filt in enumerate(filters):
                filters[i] = Filter2D(np.pad(filt.kernel, padding_psf, mode="constant"))

        if constraint_images is not None:
            for i, constraint in enumerate(constraint_images):
                constraint_images[i] = np.pad(constraint, padding, mode="wrap")

    if orthowavelets is None:
        orthowavelets = True

        for w in wavelets:
            orthowavelets &= w.orthogonal

    beta = np.zeros((len(wavelets),) + coeff_dims)
    alpha = np.zeros((len(wavelets),) + coeff_dims)
    old_beta = np.zeros((len(wavelets),) + coeff_dims)

    #just compute the slice data which will remain the same throughout as the wavelet dictionary and image dimensions don't change
    _, slices = dwt(beta[0], wavelets, parallel=parallelize_wavelets)

    t = 1

    #precomputations for gradient
    meas_term = None
    sol_term = None
    step = None

    if filters is None:
        sol_term = psf.squared(linear_conv=linear_conv)
        meas_term = util.convolve2d(dirty, psf.adj, linear=linear_conv)
        step = compute_step(psf, wavelets, orthowavelets, slices, parallelize_wavelets=parallelize_wavelets, linear_conv=linear_conv)
    else:
        curr_filter = filters[partition]

        #initialize to local fidelity term
        sol_term = util.convolve2d(util.convolve2d(psf.kernel, curr_filter.squared(linear_conv=linear_conv), linear=linear_conv), psf.adj, linear=linear_conv)
        meas_term = util.convolve2d(util.convolve2d(dirty, curr_filter.squared(linear_conv=linear_conv), linear=linear_conv), psf.adj, linear=linear_conv)

        #add additional constraint fidelity terms
        for i, constraint in enumerate(constraint_images):
            if i == partition:
                continue

            curr_filter = filters[i]

            sol_term += curr_filter.squared(linear_conv=linear_conv)
            meas_term += util.convolve2d(constraint_images[i], curr_filter.squared(linear_conv=linear_conv), linear=linear_conv)

        step = compute_step(psf, wavelets, orthowavelets, slices, filters=filters, partition=partition, parallelize_wavelets=parallelize_wavelets, linear_conv=linear_conv)

    if lambda_max is None:
        lambda_max = compute_lambda_max(meas_term, wavelets, parallelize_wavelets=parallelize_wavelets)

    lambd = lambda_max * reg_param

    curr_cost = None

    #iterate
    for i in range(niter):
        curr_img = idwt(alpha, slices, wavelets, parallel=parallelize_wavelets)

        grad_img = util.convolve2d(curr_img, sol_term, linear=linear_conv) - meas_term

        gradient = 2 * dwt(grad_img, wavelets, parallel=parallelize_wavelets)[0]

        old_beta[:,:,:] = beta
        beta[:,:,:] = soft_thresh(alpha - step * gradient, step * lambd)

        #The commented code adds in a backtracking step into fista, enforcing monotonicity, and is a variant called MFISTA, uncomment if there are stability issues causing divergence 

        # new_cost = cost(beta_img, old_beta, psf, dirty, lambd, filters=filters, constraint_images=constraint_images, partition=partition, linear_conv=linear_conv)
        # if curr_cost is None:
        #     curr_cost = new_cost
        # else:
        #     if curr_cost < new_cost:
        #         beta[:,:,:] = old_beta[:,:,:]
        #         step = step * 0.9
        #         print("recomputed step")
        #     else:
        #         curr_cost = new_cost

        new_t = (1 + np.sqrt(1 + 4 * t**2))/2

        alpha[:,:,:] = beta + (t - 1)/new_t * (beta - old_beta)

        t = new_t



    #return idwt(alpha, slices, wavelets, parallel=parallelize_wavelets)[diff_dims[0]//2:coeff_dims[0]-diff_dims[0]//2, diff_dims[1]//2:coeff_dims[1]-diff_dims[1]//2]
    return idwt(alpha, slices, wavelets, parallel=parallelize_wavelets)[0:img_dims[0], 0:img_dims[1]]


def fista_cov(psf, dirty, reg_param, wavelets, niter, vis_var, orthowavelets=None, lambda_max=None, parallelize_wavelets=False):
    """
    fista_cov uses fista to deconvolve an image using l1 regularization while taking into account estimated convariance. 
    only applied to a full-resolution case right now as this is very experimental

    :psf: point spread function
    :dirty: dirty image
    :reg_param: multiplier for lambda_max
    :wavelets: list of wavelet dictionaries
    :niter: number of fista iterations
    :vis_var: visibility variance, assumed to be constant for all visibilities
    :orthowavelets: whether individual wavelet dictionaries are orthogonal
    :lambda_max: a hard set lambda_max, leave as None if it is to be calculated
    :parallelize_wavelets: multithread wavelet transform, currently is slower than the normal version so leave this false
    :return: deconvolved image
    """
    img_dims = psf.kernel.shape

    if orthowavelets is None:
        orthowavelets = True

        for w in wavelets:
            orthowavelets &= w.orthogonal

    beta = np.zeros((len(wavelets),) + img_dims)
    alpha = np.zeros((len(wavelets),) + img_dims)
    old_beta = np.zeros((len(wavelets),) + img_dims)

    #just compute the slice data which will remain the same throughout as the wavelet dictionary and image dimensions don't change
    _, slices = dwt(beta[0], wavelets, parallel=parallelize_wavelets)

    t = 1

    cov_pi = pseudoinv(psf.kernel, 1)
    #cov_pi = np.fft.ifftshift(np.fft.ifft2(np.ones(psf.kernel.shape)))

    #precomputations for gradient
    sol_term = util.convolve2d(util.convolve2d(psf.kernel, cov_pi), psf.adj) 
    meas_term = util.convolve2d(util.convolve2d(dirty, cov_pi), psf.adj)
    step = compute_step_cov(psf, cov_pi, wavelets, orthowavelets, slices, parallelize_wavelets=parallelize_wavelets)

    if lambda_max is None:
        lambda_max = compute_lambda_max(meas_term, wavelets, parallelize_wavelets=parallelize_wavelets)

    lambd = lambda_max * reg_param

    images = []

    images.append(np.sign(np.fft.ifftshift(np.fft.fft2(cov_pi)).real) * np.abs(np.fft.ifftshift(np.fft.fft2(cov_pi)).real) ** 0.1)

    for i in range(niter):
        curr_img = idwt(alpha, slices, wavelets, parallel=parallelize_wavelets)

        images.append(curr_img)
        images.append(np.abs(np.fft.ifftshift(np.fft.fft2(curr_img)).real)**0.1)

        images.append(util.convolve2d(curr_img, sol_term))
        images.append(np.abs(np.fft.ifftshift(np.fft.fft2(util.convolve2d(curr_img, sol_term))).real)**0.1)

        grad_img = util.convolve2d(curr_img, sol_term) - meas_term

        images.append(grad_img)

        gradient = 2 * dwt(grad_img, wavelets, parallel=parallelize_wavelets)[0]

        old_beta[:,:,:] = beta
        beta[:,:,:] = soft_thresh(alpha - step * gradient, step * lambd)

        new_t = (1 + np.sqrt(1 + 4 * t**2))/2

        alpha[:,:,:] = beta + (t - 1)/new_t * (beta - old_beta)

        t = new_t

    return idwt(alpha, slices, wavelets, parallel=parallelize_wavelets), images