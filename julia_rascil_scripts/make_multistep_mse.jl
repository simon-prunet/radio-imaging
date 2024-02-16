using DeconvMultiStep
using FITSIO
using CSV
using Tables

root = "."

psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])
i_lowres = read(FITS(joinpath(root, ARGS[4]))[1])
gt = read(FITS(joinpath(root, ARGS[5]))[1])

wavelet_dict = parse(Int, ARGS[6])

# make LP and HP filters
ℓ = parse(Float64, ARGS[9])
δ = parse(Float64, ARGS[10])
n_pix, _ = size(psf)

recon_noise = parse(Float64, ARGS[7])
vis_noise = parse(Float64, ARGS[8])

len = sqrt(recon_noise * recon_noise + vis_noise * vis_noise)

recon_noise /= len
vis_noise /= len

G = make_filters(ℓ, δ, n_pix, σ² = vis_noise, η² = recon_noise)

mse = Nothing

# reconstruction
if wavelet_dict == 0
    i_multistep, mse = fista(psf, dirty, parse(Float64, ARGS[1]), 100; G=G, ip = i_lowres, sky=gt)
else
    i_multistep, m, mse = fista_iuwt(psf, dirty, parse(Float64, ARGS[1]), 100; G=G, ip = i_lowres, sky=gt)
end

CSV.write(ARGS[12], Tables.table(mse))

f = FITS(joinpath(root, ARGS[11]), "w")
write(f, i_multistep)
close(f)