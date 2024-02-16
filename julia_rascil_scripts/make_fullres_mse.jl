# reconstruction using all bases, 
# i.e. full telescope resolution

using DeconvMultiStep
using FITSIO
using CSV
using Tables

root = "."

psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])
gt = read(FITS(joinpath(root, ARGS[4]))[1])
wavelet_dict = parse(Int, ARGS[5])

mse = nothing

if wavelet_dict == 0
    i_fullres, mse = fista(psf, dirty, parse(Float64, ARGS[1]), 100, sky=gt)
else
    i_fullres, m, mse = fista_iuwt(psf, dirty, parse(Float64, ARGS[1]), 100, sky=gt)
end

CSV.write(ARGS[7], Tables.table(mse))

f = FITS(joinpath(root, ARGS[6]), "w")
write(f, i_fullres)
close(f)
