# reconstruction using all bases, 
# i.e. full telescope resolution

using DeconvMultiStep
using FITSIO
using CSV
using Tables

root = "."

lambda = parse(Float64, ARGS[1])
psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])
gt = read(FITS(joinpath(root, ARGS[4]))[1])
wavelet_dict = parse(Int, ARGS[5])
num_fista_iter = parse(Int, ARGS[6])
output_fits = joinpath(root, ARGS[7])
output_mse = joinpath(root, ARGS[8])
output_coeff_dist = joinpath(root, ARGS[9])
output_costs = joinpath(root, ARGS[10])

mse = nothing
coeff_dist = nothing

if wavelet_dict == 0
    i_fullres, coeff_dist, mse, costs = fista(psf, dirty, lambda, num_fista_iter, sky=gt)
else
    i_fullres, coeff_dist, mse, costs = fista_iuwt(psf, dirty, lambda, num_fista_iter, sky=gt)
end

CSV.write(output_mse, Tables.table(mse))
CSV.write(output_coeff_dist, Tables.table(coeff_dist))
CSV.write(output_costs, Tables.table(costs))

f = FITS(output_fits, "w")
write(f, i_fullres)
close(f)
