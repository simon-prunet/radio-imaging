using DeconvMultiStep
using FITSIO
using CSV
using Tables

root = "."

lambda = parse(Float64, ARGS[1])
psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])
i_lowres = read(FITS(joinpath(root, ARGS[4]))[1])
gt = read(FITS(joinpath(root, ARGS[5]))[1])
wavelet_dict = parse(Int, ARGS[6])
num_fista_iter = parse(Int, ARGS[7])
recon_noise = parse(Float64, ARGS[8])
vis_noise = parse(Float64, ARGS[9])

# make LP and HP filters
ℓ = parse(Float64, ARGS[10])
δ = parse(Float64, ARGS[11])
n_pix, _ = size(psf)

output_fits = joinpath(root, ARGS[12])
output_mse = joinpath(root, ARGS[13])
output_coeff_dist = joinpath(root, ARGS[14])
output_costs = joinpath(root, ARGS[15])

len = sqrt(recon_noise * recon_noise + vis_noise * vis_noise)

recon_noise /= len
vis_noise /= len

G = make_filters(ℓ, δ, n_pix, σ² = vis_noise, η² = recon_noise)

mse = Nothing

# reconstruction
if wavelet_dict == 0
    i_multistep, coeff_dist, mse, costs = fista(psf, dirty, lambda, num_fista_iter; G=G, ip = i_lowres, sky=gt)
else
    i_multistep, coeff_dist, mse, costs = fista_iuwt(psf, dirty, lambda, num_fista_iter; G=G, ip = i_lowres, sky=gt)
end

CSV.write(output_mse, Tables.table(mse))
CSV.write(output_coeff_dist, Tables.table(coeff_dist))
CSV.write(output_costs, Tables.table(costs))

f = FITS(output_fits, "w")
write(f, i_multistep)
close(f)