# reconstruction using only short bases

using DeconvMultiStep
using FITSIO
using FFTW

root = "."

psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])
wavelet_dict = parse(Int, ARGS[4])

if wavelet_dict == 0
    i_lowres = fista(psf, dirty, parse(Float64, ARGS[1]), 100)
else
    i_lowres, m = fista_iuwt(psf, dirty, parse(Float64, ARGS[1]), 100)
end

a = parse(Float64, ARGS[5])
b = 0.000000001

G = make_filters(a, b, size(i_lowres, 1); σ² = 1, η² = 1)

filtered = real(ifft(fft(i_lowres).*fft(ifftshift(G.LowPass))))

f = FITS(joinpath(root, ARGS[6]), "w")
write(f, filtered)
close(f)
