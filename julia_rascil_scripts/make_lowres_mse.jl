# reconstruction using only short bases

using DeconvMultiStep
using FITSIO
using FFTW
using CSV
using Tables

root = "."

psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])
gt = read(FITS(joinpath(root, ARGS[4]))[1])
wavelet_dict = parse(Int, ARGS[5])

mse = nothing

if wavelet_dict == 0
    i_lowres, mse = fista(psf, dirty, parse(Float64, ARGS[1]), 100, sky=gt)
else
    i_lowres, mse = fista_iuwt(psf, dirty, parse(Float64, ARGS[1]), 100, sky=gt)
end

CSV.write(ARGS[8], Tables.table(mse))

a = parse(Float64, ARGS[6],)
b = 0.000000001

G = make_filters(a, b, size(i_lowres, 1); σ² = 1, η² = 1)

filtered = real(ifft(fft(i_lowres).*fft(ifftshift(G.LowPass))))

f = FITS(joinpath(root, ARGS[7]), "w")
write(f, filtered)
close(f)
