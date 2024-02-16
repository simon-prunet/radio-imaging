using DelimitedFiles
using DeconvMultiStep
using FITSIO

root = "."

psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])

i_fullres, m = fista_iuwt(psf, dirty, parse(Float64, ARGS[1]), 10)
writedlm(joinpath(root, ARGS[4]), transpose(m), ',')