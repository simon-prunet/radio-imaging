using DeconvMultiStep
using FITSIO
using FFTW

root = "."

a = 71
b = 0.000000001

G = make_filters(a, b, 512; σ² = 1, η² = 1)

sky = read(FITS(joinpath(root, ARGS[1]))[1])
lowsky = real(ifft(fft(sky).*fft(ifftshift(G.LowPass))))

f = FITS(joinpath(root, ARGS[2]), "w")
write(f, lowsky)
close(f)