using DeconvMultiStep
using FITSIO
using FFTW

root = "."

a = parse(Float64, ARGS[2])
b = 0.000000001

sky = read(FITS(joinpath(root, ARGS[1]))[1])
G = make_filters(a, b, size(sky, 1); σ² = 1, η² = 1)

lowsky = real(ifft(fft(sky).*fft(ifftshift(G.LowPass))))

f = FITS(joinpath(root, ARGS[3]), "w")
write(f, lowsky)
close(f)