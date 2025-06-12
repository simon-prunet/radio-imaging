import numpy as np
import pywt
from joblib import Parallel, delayed
import scipy

EPSILON = 1e-6

#computes adjoint of a convolution for some given kernel, as with convolve2d, assumes circular convolution
def adj(kernel):
	shift = (1 if kernel.shape[0] % 2 == 0 else 0, 1 if kernel.shape[1] % 2 == 0 else 0)
	return np.roll(np.flip(kernel), shift=shift, axis=(0, 1)) 

#2d convolution, circular by default, but can be set to linear. This does change the lipschitz however.
def convolve2d(signal, kernel, linear=False):
	if not linear:
		return np.fft.ifft2(np.fft.fft2(signal) * np.fft.fft2(np.fft.ifftshift(kernel))).real
	else:
		pad_length = (signal.shape[0] // 2, signal.shape[1] // 2)
		sig_padded = np.pad(signal, pad_length)
		kernel_padded = np.pad(kernel, pad_length)

		return (np.fft.ifft2(np.fft.fft2(sig_padded) * np.fft.fft2(np.fft.ifftshift(kernel_padded))))[pad_length[0]:signal.shape[0]+pad_length[0], pad_length[1]:signal.shape[1]+pad_length[1]].real



#computes the pseudo-inverse for the given kernel assuming convolution
def pseudoinv(kernel, const):
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

	#lazy init to save memory in cases where this is not necessary, although this seems to not quite be the case currently
	def squared(self, linear_conv=False):
		if self.kernel2 is None:
			self.kernel2 = convolve2d(self.kernel, self.adj, linear=linear_conv)

		return self.kernel2


def dwt_single(img, wavelet):
	return pywt.coeffs_to_array(pywt.wavedec2(img, wavelet=wavelet, mode="periodization"))

def idwt_single(coeffs, slice, wavelet):
	return pywt.waverec2(pywt.array_to_coeffs(coeffs, slice, output_format="wavedec2"), wavelet=wavelet, mode="periodization")

#parallelize later after making sure this works correctly
def dwt(img, wavelets, parallel=False):
	coeffs = np.zeros((len(wavelets),) + img.shape)
	slices = [None for x in range(len(wavelets))]

	if parallel:
		results = Parallel(n_jobs=len(wavelets))(delayed(dwt_single)(img, w) for w in wavelets)

		for i, result in enumerate(results):
			coeffs[i,:,:] = result[0]
			slices[i] = result[1]
	else:
		for i, wavelet in enumerate(wavelets):
			curr_wavelet_coeffs = dwt_single(img, wavelet)
			coeffs[i,:,:] = curr_wavelet_coeffs[0]
			slices[i] = curr_wavelet_coeffs[1]

	return coeffs, slices

def idwt(coeffs, slices, wavelets, parallel=False):
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
	return np.max(np.abs(2*dwt(meas_term, wavelets, parallel=parallelize_wavelets)[0])) 

#computes the step size using the lipschitz, which in this case is the largest eigenvalue of the operator
#this is done analytically for orthogonal wavelets, or using power iteration otherwise
def compute_step(psf, wavelets, orthowavelets, wavelet_slices, filters=None, partition=None, power_iterations=50, parallelize_wavelets=False, linear_conv=False):
	#analytic
	if orthowavelets:
		lipschitz = None

		if filters is None:
			lipschitz = np.max(np.abs(np.fft.fft2(psf.squared(linear_conv=linear_conv))))
		else:
			fpsf = np.fft.fft2(psf.kernel)
			ffilt = np.fft.fft2(filters[partition].kernel)
			lipschitz = np.max(np.abs(fpsf*fpsf*ffilt*ffilt))

			for i, filt in enumerate(filters):
				if i == partition:
					continue

				lipschitz += np.max(np.abs(filters[i].squared(linear_conv=linear_conv)))


		return 1 / (2 * len(wavelets) * lipschitz)
	#power iteration
	else:
		operator_kernel = None
		if filters is None:
			operator_kernel = psf.squared(linear_conv=linear_conv)
		else:
			operator_kernel = convolve2d(convolve2d(psf, curr_filter.squared(linear_conv=linear_conv), linear=linear_conv), psf.adj, linear=linear_conv)

		for i, filt in enumerate(filters):
			if i == partition:
				continue

			operator_kernel += filt.squared(linear_conv=linear_conv)

		img_dims = psf.kernel.shape
		alpha = np.random.rand((len(wavelets),) + img_dims)

		for i in range(power_iterations):
			alpha = dwt(convolve2d(idwt(alpha, wavelet_slices, wavelets, parallel=parallelize_wavelets), operator_kernel, linear=linear_conv), wavelets, parallel=parallelize_wavelets)[0]
			alpha = alpha / np.linalg.norm(alpha)

		alpha_new = dwt(convolve2d(idwt(alpha, wavelet_slices, wavelets, parallel=parallelize_wavelets), operator_kernel, linear=linear_conv), wavelets, parallel=parallelize_wavelets)[0]

		return np.dot(alpha.T, alpha_new)

def compute_step_cov(psf, cov_pi, wavelets, orthowavelets, wavelet_slices, power_iterations=50, parallelize_wavelets=False, linear_conv=False):
	#analytic
	if orthowavelets:
		lipschitz = np.max(np.abs(np.fft.fft2(convolve2d(convolve2d(psf.kernel, cov_pi, linear=linear_conv), psf.adj, linear=linear_conv))))

		return 1 / (2 * len(wavelets) * lipschitz)
	#power iteration
	else:
		operator_kernel = convolve2d(convolve2d(psf.kernel, cov_pi, linear=linear_conv), psf.adj, linear=linear_conv)

		img_dims = psf.kernel.shape
		alpha = np.random.rand((len(wavelets),) + img_dims)

		for i in range(power_iterations):
			alpha = dwt(convolve2d(idwt(alpha, wavelet_slices, wavelets, parallel=parallelize_wavelets), operator_kernel, linear=linear_conv), wavelets, parallel=parallelize_wavelets)[0]
			alpha = alpha / np.linalg.norm(alpha)

		alpha_new = dwt(convolve2d(idwt(alpha, wavelet_slices, wavelets, parallel=parallelize_wavelets), operator_kernel, linear=linear_conv), wavelets, parallel=parallelize_wavelets)[0]

		return np.dot(alpha.T, alpha_new)

def soft_thresh(coeffs, step):
	return np.sign(coeffs) * np.maximum(np.abs(coeffs) - step, 0)

#fista for solving multi-partition deconvolution
def fista(psf, dirty, reg_param, wavelets, niter, orthowavelets=None, filters=None, constraint_images=None, partition=None, lambda_max=None, parallelize_wavelets=False, linear_conv=False):
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

	#precomputations for gradient
	meas_term = None
	sol_term = None
	step = None

	if filters is None:
		sol_term = psf.squared(linear_conv=linear_conv)
		meas_term = convolve2d(dirty, psf.adj, linear=linear_conv)
		step = compute_step(psf, wavelets, orthowavelets, slices, parallelize_wavelets=parallelize_wavelets, linear_conv=linear_conv)
	else:
		curr_filter = filters[partition]

		#initialize to local fidelity term
		sol_term = convolve2d(convolve2d(psf.kernel, curr_filter.squared(linear_conv=linear_conv), linear=linear_conv), psf.adj, linear=linear_conv)
		meas_term = convolve2d(convolve2d(dirty, curr_filter.squared(linear_conv=linear_conv), linear=linear_conv), psf.adj, linear=linear_conv)

		#add additional constraint fidelity terms
		for i, constraint in enumerate(constraint_images):
			if i == partition:
				continue

			curr_filter = filters[i]

			sol_term += curr_filter.squared(linear_conv=linear_conv)
			meas_term += convolve2d(constraint_images[i], curr_filter.squared(linear_conv=linear_conv), linear=linear_conv)

		step = compute_step(psf, wavelets, orthowavelets, slices, filters=filters, partition=partition, parallelize_wavelets=parallelize_wavelets, linear_conv=linear_conv)

	if lambda_max is None:
		lambda_max = compute_lambda_max(meas_term, wavelets, parallelize_wavelets=parallelize_wavelets)

	lambd = lambda_max * reg_param

	#iterate
	for i in range(niter):
		curr_img = idwt(alpha, slices, wavelets, parallel=parallelize_wavelets)

		grad_img = convolve2d(curr_img, sol_term, linear=linear_conv) - meas_term

		gradient = 2 * dwt(grad_img, wavelets, parallel=parallelize_wavelets)[0]

		old_beta[:,:,:] = beta
		beta[:,:,:] = soft_thresh(alpha - step * gradient, step * lambd)

		new_t = (1 + np.sqrt(1 + 4 * t**2))/2

		alpha[:,:,:] = beta + (t - 1)/new_t * (beta - old_beta)

		t = new_t

	return idwt(alpha, slices, wavelets, parallel=parallelize_wavelets)


#fista taking into account covariance approximation, just a prototype for now
def fista_cov(psf, dirty, reg_param, wavelets, niter, vis_var, orthowavelets=None, lambda_max=None, parallelize_wavelets=False):
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

	cov_pi = pseudoinv(psf.kernel, ivs_var)
	#cov_pi = np.fft.ifftshift(np.fft.ifft2(np.ones(psf.kernel.shape)))

	#precomputations for gradient
	sol_term = convolve2d(convolve2d(psf.kernel, cov_pi), psf.adj) 
	meas_term = convolve2d(convolve2d(dirty, cov_pi), psf.adj)
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

		images.append(convolve2d(curr_img, sol_term))
		images.append(np.abs(np.fft.ifftshift(np.fft.fft2(convolve2d(curr_img, sol_term))).real)**0.1)

		grad_img = convolve2d(curr_img, sol_term) - meas_term

		images.append(grad_img)

		gradient = 2 * dwt(grad_img, wavelets, parallel=parallelize_wavelets)[0]

		old_beta[:,:,:] = beta
		beta[:,:,:] = soft_thresh(alpha - step * gradient, step * lambd)

		new_t = (1 + np.sqrt(1 + 4 * t**2))/2

		alpha[:,:,:] = beta + (t - 1)/new_t * (beta - old_beta)

		t = new_t

	return idwt(alpha, slices, wavelets, parallel=parallelize_wavelets), images