using DeconvMultiStep
using FITSIO
using FFTW

root = "."

i_prev_lowres = read(FITS(joinpath(root, ARGS[1]))[1])
i_prev_deconv = read(FITS(joinpath(root, ARGS[2]))[1])

# make lowres constraint, no filtering as that is done in the actual second step
i_constraint = i_prev_lowres - i_prev_deconv

f = FITS(joinpath(root, ARGS[3]), "w")
write(f, i_constraint)
close(f)