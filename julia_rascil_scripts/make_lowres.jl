# reconstruction using only short bases

using DeconvMultiStep
using FITSIO
using FFTW

root = "."

lambda = parse(Float64, ARGS[1])
psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])
wavelet_dict = parse(Int, ARGS[4])
num_fista_iter = parse(Int, ARGS[5])
output_filename = joinpath(root, ARGS[6])

if wavelet_dict == 0
    i_lowres = fista(psf, dirty, lambda, num_fista_iter)
else
    i_lowres = fista_iuwt(psf, dirty, lambda, num_fista_iter)
end

f = FITS(output_filename, "w")
write(f, i_lowres)
close(f)
