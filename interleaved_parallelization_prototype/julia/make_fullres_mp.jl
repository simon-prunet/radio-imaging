# reconstruction using all bases, 
# i.e. full telescope resolution

using DeconvMultiStep
using FITSIO

root = "."

lambda = parse(Float64, ARGS[1])
psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])
wavelet_dict = parse(Int, ARGS[4])
num_fista_iter = parse(Int, ARGS[5])
output_filename = joinpath(root, ARGS[6])

i_fullres = fista_multipartition(psf, dirty, lambda, num_fista_iter)

f = FITS(output_filename, "w")
write(f, i_fullres)
close(f)
