using DeconvMultiStep
using FITSIO
using FFTW

root = "."

i_prev_lowres = read(FITS(joinpath(root, ARGS[1]))[1])
i_prev_deconv = read(FITS(joinpath(root, ARGS[2]))[1])

# make LP and HP filters
ℓ = 61
δ = 10
n_pix, _ = size(i_prev_lowres)

G = make_filters(ℓ, δ, n_pix)

# make lowres constraint
i_constraint = i_prev_lowres - real(ifft(fft(i_prev_deconv).*fft(ifftshift(G.LowPass))))

f = FITS(joinpath(root, ARGS[3]), "w")
write(f, i_constraint)
close(f)