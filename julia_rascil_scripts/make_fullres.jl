# reconstruction using all bases, 
# i.e. full telescope resolution

using DeconvMultiStep
using FITSIO

root = "."

psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])
wavelet_dict = parse(Int, ARGS[4])

if wavelet_dict == 0
    i_fullres = fista(psf, dirty, parse(Float64, ARGS[1]), 100)
else
    i_fullres, m = fista_iuwt(psf, dirty, parse(Float64, ARGS[1]), 100)
end

f = FITS(joinpath(root, ARGS[5]), "w")
write(f, i_fullres)
close(f)
