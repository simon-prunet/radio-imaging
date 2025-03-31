import numpy as np
import pywt


#computes adjoint of a convolution for some given kernel, as with convolve2d, assumes circular convolution
def adj(kernel):
	shift = (1 if kernel.shape[0] % 2 == 0 else 0, 1 if kernel.shape[1] % 2 == 0 else 0)
	return np.roll(np.flip(kernel), shift=shift, axis=(0, 1)) 

#circular, add padding later if this is an issue, step size should still be the same as the max eigenvalue is the same, will need to update some other functions however (e.g. adj) accordingly
def convolve2d(signal, kernel):
	return np.fft.ifft2(np.fft.fft2(signal) * np.fft.fft2(np.fft.ifftshift(kernel))).real

class Filter2D:
	def __init__(self, kernel):
		self.kernel = kernel
		self.adj = adj(kernel)
		self.kernel2 = None

	#lazy init to save memory in cases where this is not necessary, although this seems to not quite be the case currently
	def squared(self):
		if self.kernel2 is None:
			self.kernel2 = convolve2d(self.kernel, self.adj)

		return self.kernel2


#parallelize later after making sure this works correctly
def dwt(img, wavelets):
	coeffs = np.zeros((len(wavelets),) + img.shape)
	slices = []

	for i, wavelet in enumerate(wavelets):
		curr_wavelet_coeffs = pywt.coeffs_to_array(pywt.wavedec2(img, wavelet=wavelet, mode="periodization"))
		coeffs[i,:,:] = curr_wavelet_coeffs[0]
		slices.append(curr_wavelet_coeffs[1])

	return coeffs, slices

def idwt(coeffs, slices, wavelets):
	img = np.zeros(coeffs[0].shape)

	for i, wavelet in enumerate(wavelets):
		img += pywt.waverec2(pywt.array_to_coeffs(coeffs[i], slices[i], output_format="wavedec2"), wavelet=wavelet, mode="periodization")

	return img


def compute_lambda_max(meas_term, wavelets):
	return np.max(np.abs(2*dwt(meas_term, wavelets)[0])) 

#computes the step size using the lipschitz, which in this case is the largest eigenvalue of the operator
#this is done analytically for orthogonal wavelets, or using power iteration otherwise
def compute_step(psf, wavelets, orthowavelets, wavelet_slices, filters=None, partition=None, power_iterations=50):
	#analytic
	if orthowavelets:
		lipschitz = None

		if filters is None:
			lipschitz = np.max(np.abs(np.fft.fft2(psf.squared())))
		else:
			fpsf = np.fft.fft2(psf.kernel)
			ffilt = np.fft.fft2(filters[partition].kernel)
			lipschitz = np.max(np.abs(fpsf*fpsf*ffilt*ffilt))

			for i, filt in enumerate(filters):
				if i == partition:
					continue

				lipschitz += np.max(np.abs(filters[i].squared()))


		return 1 / (2 * len(wavelets) * lipschitz)
	#power iteration
	else:
		operator_kernel = None
		if filters is None:
			operator_kernel = psf.squared()
		else:
			operator_kernel = convolve2d(convolve2d(psf, curr_filter.squared()), psf.adj)

		for i, filt in enumerate(filters):
			if i == partition:
				continue

			operator_kernel += filt.squared()

		img_dims = psf.kernel.shape
		alpha = np.random.rand((len(wavelets),) + img_dims)

		for i in range(power_iterations):
			alpha = dwt(convolve2d(idwt(alpha, wavelet_slices, wavelets), operator_kernel), wavelets)[0]
			alpha = alpha / np.linalg.norm(alpha)

		alpha_new = dwt(convolve2d(idwt(alpha, wavelet_slices, wavelets), operator_kernel), wavelets)[0]

		return np.dot(alpha.T, alpha_new)


def soft_thresh(coeffs, step):
	return np.sign(coeffs) * np.maximum(np.abs(coeffs) - step, 0)

def fista(psf, dirty, reg_param, wavelets, niter, orthowavelets=None, filters=None, constraint_images=None, partition=None, lambda_max=None):
	img_dims = psf.kernel.shape

	if orthowavelets is None:
		orthowavelets = True

		for w in wavelets:
			orthowavelets &= w.orthogonal

	beta = np.zeros((len(wavelets),) + img_dims)
	alpha = np.zeros((len(wavelets),) + img_dims)
	old_beta = np.zeros((len(wavelets),) + img_dims)

	#just compute the slice data which will remain the same throughout as the wavelet dictionary and image dimensions don't change
	_, slices = dwt(beta[0], wavelets)

	t = 1

	#precomputations for gradient
	meas_term = None
	sol_term = None
	step = None

	if filters is None:
		sol_term = psf.squared()
		meas_term = convolve2d(dirty, psf.adj)
		step = compute_step(psf, wavelets, orthowavelets, slices)
	else:
		curr_filter = filters[partition]

		#initialize to local fidelity term
		sol_term = convolve2d(convolve2d(psf, curr_filter.squared()), psf.adj)
		meas_term = convolve2d(convolve2d(dirty, curr_filter.squared()), psf.adj)

		#add additional constraint fidelity terms
		for i, constraint in enumerate(constraint_images):
			if i == partition:
				continue

			curr_filter = filters[partition]

			sol_term += curr_filter.squared()
			meas_term += convolve2d(constraint_images[i], curr_filter.squared())

		step = compute_step(psf, wavelets, orthowavelets, slices, filters=filters, partition=partition)

	if lambda_max is None:
		lambda_max = compute_lambda_max(meas_term, wavelets)

	lambd = lambda_max * reg_param

	#iterate
	for i in range(niter):
		curr_img = idwt(alpha, slices, wavelets)

		grad_img = convolve2d(curr_img, sol_term) - meas_term

		gradient = 2 * dwt(grad_img, wavelets)[0]

		for j in range(len(beta)):
			old_beta[j][:,:] = beta[j]
			beta[j][:,:] = soft_thresh(alpha[j] - step * gradient[j], step * lambd)

		new_t = (1 + np.sqrt(1 + 4 * t**2))/2

		for j in range(len(alpha)):
			alpha[j][:,:] = beta[j] + (t - 1)/new_t * (beta[j] - old_beta[j])

	return idwt(alpha, slices, wavelets)