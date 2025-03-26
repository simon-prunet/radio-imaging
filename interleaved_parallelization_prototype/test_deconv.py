import ip_helpers as iph

residual = iph.fromfits("10ktest/residual_0_1.fits")
deconv_vl = iph.fromfits("10ktest/deconv_vl_0.fits")
deconv_vh = iph.fromfits("10ktest/deconv_vh_0.fits")
constraint = [deconv_vl, deconv_vh]
psf = iph.fromfits("10ktest/psf_0.fits")

deconvolved = iph.deconvolve(0, residual, psf, constraint, 50, 0, 1, 0.01, 2, 40, 1, 5, 1000)

iph.tofits(deconvolved, "10ktest/deconvolved.fits")