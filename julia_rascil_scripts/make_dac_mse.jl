using DeconvMultiStep
using FITSIO
using CSV
using Tables

root = "."

lambda = parse(Float64, ARGS[1])
psf_high = read(FITS(joinpath(root, ARGS[2]))[1])
dirty_high = read(FITS(joinpath(root, ARGS[3]))[1])
psf_low = read(FITS(joinpath(root, ARGS[4]))[1])
dirty_low = read(FITS(joinpath(root, ARGS[5]))[1])
gt = read(FITS(joinpath(root, ARGS[6]))[1])
wavelet_dict = parse(Int, ARGS[7])
num_fista_iter = parse(Int, ARGS[8])
eta = parse(Float64, ARGS[9])
sigma = parse(Float64, ARGS[10])

# make LP and HP filters
ℓ = parse(Float64, ARGS[11])
δ = parse(Float64, ARGS[12])
n_pix, _ = size(psf_high)

output_fits = joinpath(root, ARGS[13])
output_mse = joinpath(root, ARGS[14])
output_coeff_dist = joinpath(root, ARGS[15])
output_costs = joinpath(root, ARGS[16])

len = sqrt(eta * eta + sigma * sigma)

#eta /= len
#sigma /= len

G = make_filters(ℓ, δ, n_pix, σ² = sigma, η² = eta)

# reconstruction
if wavelet_dict == 0
    i_multistep, coeff_dist, mse, costs = fista(psf_high, dirty_high, lambda, num_fista_iter; G=G, ip = dirty_low, sky=gt, Hl = psf_low)
else
    i_multistep, coeff_dist, mse, costs = fista_iuwt(psf_high, dirty_high, lambda, num_fista_iter; G=G, ip = dirty_low, sky=gt, Hl = psf_low)
end

CSV.write(output_mse, Tables.table(mse))
CSV.write(output_coeff_dist, Tables.table(coeff_dist))
CSV.write(output_costs, Tables.table(costs))

f = FITS(output_fits, "w")
write(f, i_multistep)
close(f)