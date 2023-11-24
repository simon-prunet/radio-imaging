using DeconvMultiStep
using FITSIO

root = "."

psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])
i_lowres = read(FITS(joinpath(root, ARGS[4]))[1])

# make LP and HP filters
ℓ = 61
δ = 10
n_pix, _ = size(psf)

G = make_filters(ℓ, δ, n_pix)

# reconstruction
i_multistep = fista(psf, dirty, parse(Float64, ARGS[1]), 100; G=G, ip = i_lowres)

f = FITS(joinpath(root, ARGS[5]), "w")
write(f, i_multistep)
close(f)