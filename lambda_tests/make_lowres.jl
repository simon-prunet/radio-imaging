# reconstruction using only short bases

using DeconvMultiStep
using FITSIO
using FFTW

root = "."

psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])

i_lowres = fista(psf, dirty, parse(Float64, ARGS[1]), 100)

a = 71
b = 0.000000001

G = make_filters(a, b, 512; σ² = 1, η² = 1)

filtered = real(ifft(fft(i_lowres).*fft(ifftshift(G.LowPass))))

f = FITS(joinpath(root, ARGS[4]), "w")
write(f, filtered)
close(f)
